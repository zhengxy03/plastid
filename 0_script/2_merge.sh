#!/usr/bin/env bash

BASH_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "${BASH_DIR}"

if [[ "$PWD" == *0_script* ]]; then
    cd ..
fi

#----------------------------#
# Colors in term
#----------------------------#
# http://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
GREEN=
RED=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    NC='\033[0m' # No Color
fi

log_warn () {
    echo >&2 -e "${RED}==> $@ <==${NC}"
}

log_info () {
    echo >&2 -e "${GREEN}==> $@${NC}"
}

log_debug () {
    echo >&2 -e "==> $@"
}

#----------------------------#
# helper functions
#----------------------------#
set +e

# set stacksize to unlimited
if [[ "$OSTYPE" != "darwin"* ]]; then
    ulimit -s unlimited
fi

signaled () {
    log_warn Interrupted
    exit 1
}
trap signaled TERM QUIT INT

# save environment variables
save () {
    printf ". + { %s: \"%s\"}" $1 $(eval "echo -n \"\$$1\"") > jq.filter.txt

    if [ -e env.json ]; then
        cat env.json |
            jq --sort-keys --from-file jq.filter.txt \
            > env.json.new
        rm env.json
    else
        jq --from-file jq.filter.txt --null-input \
            > env.json.new
    fi

    mv env.json.new env.json
    rm jq.filter.txt
}

stat_format () {
    echo $(faops n50 -H -N 50 -S -C $@) |
        perl -nla -MNumber::Format -e '
            printf qq(%d\t%s\t%d\n), $F[0], Number::Format::format_bytes($F[1], base => 1000,), $F[2];
        '
}

byte_format () {
    echo "$@" |
        perl -nl -MNumber::Format -e '
            print Number::Format::format_bytes($_, base => 1000,);
        '
}

time_format () {
    echo "$@" |
        perl -nl -e '
            sub parse_duration {
                use integer;
                sprintf("%d:%02d:%02d", $_[0]/3600, $_[0]/60%60, $_[0]%60);
            }
            print parse_duration($_);
        '
}

readlinkf () {
    perl -MCwd -l -e 'print Cwd::abs_path shift' "$1";
}

#----------------------------#
# run
#----------------------------#
log_warn 2_merge.sh

if [ -e 2_illumina/merge/pe.cor.fa.gz ]; then
    log_debug "2_illumina/merge/pe.cor.fa.gz presents"
    exit;
fi

mkdir -p 2_illumina/merge
cd 2_illumina/merge

START_TIME=$(date +%s)
save START_TIME

NUM_THREADS=8
save NUM_THREADS

# save genome size
ESTIMATED_GENOME_SIZE=580076
save ESTIMATED_GENOME_SIZE

for PREFIX in R S T; do
    if [ ! -e ../trim/${PREFIX}1.fq.gz ]; then
        continue;
    fi

    # rotate ascii characters https://en.wikipedia.org/wiki/ROT13
    PREFIXM=$(echo ${PREFIX} | tr 'A-Z' 'V-ZA-U')   # M N O
    PREFIXU=$(echo ${PREFIX} | tr 'A-Z' 'D-ZA-C')   # U V W

    if [ -e ${PREFIX}1.fq.gz ]; then
        log_debug "2_illumina/merge/${PREFIXM}1.fq.gz presents"
        continue;
    fi

    anchr mergeread \
        ../trim/${PREFIX}1.fq.gz ../trim/${PREFIX}2.fq.gz ../trim/${PREFIX}s.fq.gz \
        --ecphase "1 2 3" \
        --parallel 8 --xmx 12g \
        --prefixm ${PREFIXM} \
        --prefixu ${PREFIXU} \
        -o mergeread.sh
    bash mergeread.sh

    # Create .cor.fa.gz
    # 合并双端未配对序列（PE）
    faops interleave -p "unmerged" \
        ${PREFIXU}1.fq.gz \
        ${PREFIXU}2.fq.gz \
        > ${PREFIXM}.interleave.fa

# 追加单端未配对序列（SE）
    faops interleave -p "single" \
        ${PREFIXU}s.fq.gz \
        >> ${PREFIXM}.interleave.fa

# 追加合并后的双端序列（Merged PE）
    faops interleave -p "merged" \
        ${PREFIXM}1.fq.gz \
        >> ${PREFIXM}.interleave.fa

    # Shuffle interleaved reads.
    log_info Shuffle interleaved reads.
    cat ${PREFIXM}.interleave.fa |
        awk '{
            OFS="\t"; \
            getline seq; \
            getline name2; \
            getline seq2; \
            print $0,seq,name2,seq2}' |
        tsv-sample |
        awk '{OFS="\n"; print $1,$2,$3,$4}' \
        > ${PREFIXM}.cor.fa
    rm ${PREFIXM}.interleave.fa
    pigz -p 8 ${PREFIXM}.cor.fa

    log_info "stats of all .fq.gz files"
    if [ ! -e statMergeReads.tsv ]; then
        printf "%s\t%s\t%s\t%s\n" \
            "Name" "N50" "Sum" "#" \
            > statMergeReads.tsv
    fi

    for NAME in clumped ecco eccc ecct extended merged.raw unmerged.raw unmerged.trim ${PREFIXM}1 ${PREFIXU}1 ${PREFIXU}2 ${PREFIXU}s; do
        if [ ! -e ${NAME}.fq.gz ]; then
            continue
        fi
        printf "%s\t%s\t%s\t%s\n" \
            $(echo ${NAME}; stat_format ${NAME}.fq.gz;) >> statMergeReads.tsv
    done

    printf "%s\t%s\t%s\t%s\n" \
        $(echo ${PREFIXM}.cor; stat_format ${PREFIXM}.cor.fa.gz;) >> statMergeReads.tsv

    log_info "stats of insert sizes"
    if [ ! -e statMergeInsert.tsv ]; then
        printf "%s\t%s\t%s\t%s\t%s\n" \
            "Group" "Mean" "Median" "STDev" "Pairs%" \
            > statMergeInsert.tsv
    fi

    #Mean	339.868
    #Median	312
    #Mode	251
    #STDev	134.676
    #PercentOfPairs	36.247
    for NAME in ${PREFIXM}.ihist.merge1.txt ${PREFIXM}.ihist.merge.txt; do
        printf "| %s " ${NAME} >> statMergeReads.md
        cat ${NAME} |
            GROUP="${NAME}" perl -nla -e '
                BEGIN { our $stat = { }; };

                m{\#(Mean|Median|STDev|PercentOfPairs)} or next;
                $stat->{$1} = $F[1];

                END {
                    printf qq(%s\t%.1f\t%s\t%.1f\t%.2f%%\n),
                        qq($ENV{GROUP}),
                        $stat->{Mean},
                        $stat->{Median},
                        $stat->{STDev},
                        $stat->{PercentOfPairs};
                }
                ' \
            >> statMergeInsert.tsv
    done

    log_info "clear unneeded .fq.gz files"
    for NAME in temp clumped ecco eccc ecct extended merged.raw unmerged.raw unmerged.trim; do
        if [ -e ${NAME}.fq.gz ]; then
            rm ${NAME}.fq.gz
        fi
    done

done

log_debug "Combine .cor.fa.gz files"
if [ -e ../S1.fq.gz ]; then
    gzip -d -c [MNO].cor.fa.gz |
        awk '{
            OFS="\t"; \
            getline seq; \
            getline name2; \
            getline seq2; \
            print $0,seq,name2,seq2}' |
        tsv-sample |
        awk '{OFS="\n"; print $1,$2,$3,$4}' \
        > pe.cor.fa
    pigz -p 8 pe.cor.fa
    rm [MNO].cor.fa.gz
else
    mv M.cor.fa.gz pe.cor.fa.gz
fi

log_debug "Reads stats with hnsm"
SUM_OUT=$( faops n50 -H -N 0 -S pe.cor.fa.gz )
save SUM_OUT

cat statMergeReads.tsv |
    rgr md stdin --right 2-4 \
    > statMergeReads.md

echo -e "\nTable: statMergeReads\n" >> statMergeReads.md

cat statMergeReads.md
mv statMergeReads.md ${BASH_DIR}/../9_markdown

cat statMergeInsert.tsv |
    rgr md stdin --right 2-5 \
    > statMergeInsert.md

echo -e "\nTable: statMergeInsert\n" >> statMergeInsert.md

cat statMergeInsert.md
mv statMergeInsert.md ${BASH_DIR}/../9_markdown

END_TIME=$(date +%s)
save END_TIME
RUNTIME=$((END_TIME-START_TIME))
save RUNTIME

log_info Done.

exit 0


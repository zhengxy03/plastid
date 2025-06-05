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
log_warn 2_trim.sh

mkdir -p 2_illumina/trim
cd 2_illumina/trim

for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ -e ${PREFIX}1.fq.gz ]; then
        log_debug "2_illumina/trim/${PREFIX}1.fq.gz presents"
        continue;
    fi

    anchr trim \
        --dedupe --cutoff 30 --cutk 31 \
        --qual "25 30" \
        --len "60" \
    --filter "adapter artifact" \
    --parallel 8 --xmx 12g \
        ../${PREFIX}1.fq.gz ../${PREFIX}2.fq.gz \
        --prefix ${PREFIX} \
        -o trim.sh
    bash trim.sh

    log_info "stats of all .fq.gz files"

    if [ ! -e statTrimReads.tsv ]; then
        printf "%s\t%s\t%s\t%s\n" \
            "Name" "N50" "Sum" "#" \
            > statTrimReads.tsv
    fi

    for NAME in clumpify filteredbytile highpass sample trim filter ${PREFIX}1 ${PREFIX}2 ${PREFIX}s; do
        if [ ! -e ${NAME}.fq.gz ]; then
            continue;
        fi

        printf "%s\t%s\t%s\t%s\n" \
            $(echo ${NAME}; stat_format ${NAME}.fq.gz;) >> statTrimReads.tsv
    done

    log_info "clear unneeded .fq.gz files"
    for NAME in temp clumpify filteredbytile highpass sample trim filter; do
        if [ -e ${NAME}.fq.gz ]; then
            rm ${NAME}.fq.gz
        fi
    done
done

cat statTrimReads.tsv |
    rgr md stdin --right 2-4 \
    > statTrimReads.md

echo -e "\nTable: statTrimReads\n" >> statTrimReads.md

for PREFIX in R S T; do
    if [ ! -s statTrimReads.md ]; then
        continue;
    fi

    if [ -e ${PREFIX}.trim.stats.txt ]; then
        echo >> statTrimReads.md
        echo '```text' >> statTrimReads.md
        echo "#${PREFIX}.trim" >> statTrimReads.md
        cat ${PREFIX}.trim.stats.txt |
            perl -nla -F"\t" -e '
                /^#(Matched|Name)/ and print and next;
                /^#/ and next;
                $F[2] =~ m{([\d.]+)} and $1 > 0.1 and print;
            ' \
            >> statTrimReads.md
        echo '```' >> statTrimReads.md
    fi

    if [ -e ${PREFIX}.filter.stats.txt ]; then
        echo >> statTrimReads.md
        echo '```text' >> statTrimReads.md
        echo "#${PREFIX}.filter" >> statTrimReads.md
        cat ${PREFIX}.filter.stats.txt |
            perl -nla -F"\t" -e '
                /^#(Matched|Name)/ and print and next;
                /^#/ and next;
                $F[2] =~ m{([\d.]+)} and $1 > 0.01 and print;
            ' \
            >> statTrimReads.md
        echo '```' >> statTrimReads.md
    fi

#    if [ -e ${PREFIX}.peaks.txt ]; then
#        echo >> statTrimReads.md
#        echo '```text' >> statTrimReads.md
#        echo "#${PREFIX}.peaks" >> statTrimReads.md
#        cat ${PREFIX}.peaks.txt |
#            grep "^#" \
#            >> statTrimReads.md
#        echo '```' >> statTrimReads.md
#    fi
done

if [ -s statTrimReads.md ]; then
    cat statTrimReads.md
    mv statTrimReads.md ${BASH_DIR}/../9_markdown
fi

cd ${BASH_DIR}/..
cd 2_illumina

parallel --no-run-if-empty --linebuffer -k -j 2 "
    ln -fs ./trim/Q{1}L{2}/ ./Q{1}L{2}
    " ::: 25 30 ::: 60
ln -fs ./trim ./Q0L0

log_info Done.

exit 0


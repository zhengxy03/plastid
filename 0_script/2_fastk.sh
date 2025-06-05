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
# Run
#----------------------------#
log_warn 2_fastk.sh

mkdir -p 2_illumina/fastk
cd 2_illumina/fastk

for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ -d ${PREFIX}-GeneScope-21 ]; then
        continue;
    fi

    for KMER in 21 51 81; do
        log_info "PREFIX: ${PREFIX}; KMER: ${KMER}"

        log_info "FastK"
        FastK -v -T8 -t1 -k${KMER} \
            ../${PREFIX}1.fq.gz ../${PREFIX}2.fq.gz \
            -NTable-${KMER}

        log_info "GeneScope"
        Histex -G Table-${KMER} |
            Rscript ../../0_script/genescopefk.R -k ${KMER} -p 1 -o ${PREFIX}-GeneScope-${KMER}

        KatGC -T8 -x1.9 -s Table-${KMER} ${PREFIX}-Merqury-KatGC-${KMER}

        Fastrm Table-${KMER}
    done
done

for PREFIX in R S T; do
    for KMER in 21 51 81; do
        if [ ! -e ${PREFIX}-GeneScope-${KMER}/summary.txt ]; then
            continue
        fi

        COV=$(
            cat ${PREFIX}-GeneScope-${KMER}/model.txt |
                grep '^kmercov' |
                tr -s ' ' '\t' |
                cut -f 2 |
                perl -nl -e 'printf qq(%.1f\n), $_'
        )

        cat ${PREFIX}-GeneScope-${KMER}/summary.txt |
            sed '1,6 d' |
            sed "1 s/^/K\t/" |
            sed "2 s/^/${PREFIX}.${KMER}\t/" |
            sed "3,7 s/^/\t/" |
            perl -nlp -e 's/\s{2,}/\t/g; s/\s+$//g;' |
            perl -nla -F'\t' -e '
                @fields = map {/\bNA\b/ ? q() : $_ } @F;        # Remove NA fields
                $fields[2] = q() if $fields[2] eq $fields[3];   # Remove identical fields
                print join qq(\t), @fields
            '

        printf "\tKmer Cov\t\t${COV}\n"
    done
done |
    keep-header -- grep -v 'property' \
    > statFastK.tsv

cat statFastK.tsv |
    rgr md stdin --right 3-4 \
    > statFastK.md

echo -e "\nTable: statFastK\n" >> statFastK.md

cat statFastK.md
mv statFastK.md ${BASH_DIR}/../9_markdown


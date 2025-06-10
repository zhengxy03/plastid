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
log_warn 1_repetitive.sh

mkdir -p 1_genome/repetitive
cd 1_genome/repetitive

if [ ! -s repetitive.fa ]; then
    hnsm size ../genome.fa > chr.sizes

    FastK -v -p -k21 ../genome

    cat chr.sizes |
        number-lines |
        parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
            Profex ../genome {1} |
                sed "1,2 d" |
                perl -nl -e '\''/(\d+).+(\d+)/ and printf qq(%s\t%s\n), $1 + 1, $2'\'' |
                tsv-filter --ge 2:2 |
                cut -f 1 |
                sed "s/^/{2}:/" |
                spanr cover stdin |
                spanr span --op fill -n 2 stdin |
                spanr span --op excise -n 100 stdin |
                spanr span --op fill -n 10 stdin
        ' |
        spanr combine stdin \
        > repetitive.json

    Fastrm ../genome

    spanr convert repetitive.json > region.txt
    hnsm range ../genome.fa -r region.txt |
        hnsm filter -N -d -a 100 stdin \
        > repetitive.fa

    spanr stat chr.sizes repetitive.json |
        tr ',' '\t' \
        > statRepetitive.tsv
fi

cat statRepetitive.tsv |
    rgr md stdin --num \
    > statRepetitive.md

echo -e "\nTable: statRepetitive\n" >> statRepetitive.md

cat statRepetitive.md
mv statRepetitive.md ${BASH_DIR}/../9_markdown


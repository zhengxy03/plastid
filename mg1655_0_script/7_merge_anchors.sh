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
log_warn 7_merge_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [DIR_PREFIX] [DIR_MERGE]"

DIR_PREFIX=${1:-"4_unitigs"}
DIR_MERGE=${2:-"7_merge_anchors"}

if [ -e ${DIR_MERGE}/anchor.merge.fasta ]; then
    echo >&2 "${DIR_MERGE}/anchor.merge.fasta presents"
    exit;
fi

#----------------------------#
# merge anchors
#----------------------------#
log_info "anchor.non-contained"

mkdir -p ${DIR_MERGE}

# reversely sorted files, so that Q30L60X80 will be infile_0
anchr contained \
    $( find . -path "*${DIR_PREFIX}*" -name "anchor.fasta" -or -path "*${DIR_PREFIX}*" -name "anchor.merge.fasta" | sort -r ) \
    --len 1000 --idt 0.9999 --ratio 0.99999 --parallel 24 \
    -o stdout |
    hnsm filter -a 1000 stdin -o ${DIR_MERGE}/anchor.non-contained.fasta

anchr orient \
    ${DIR_MERGE}/anchor.non-contained.fasta \
    --len 1000 --idt 0.999 --parallel 24 \
    -o ${DIR_MERGE}/anchor.intermediate_0.fasta
anchr merge \
    ${DIR_MERGE}/anchor.intermediate_0.fasta \
    --len 1000 --idt 0.9999 --parallel 24 \
    -o ${DIR_MERGE}/anchor.intermediate_1.fasta
anchr contained \
    ${DIR_MERGE}/anchor.intermediate_1.fasta \
    --len 1000 --idt 0.98 --ratio 0.99 --parallel 24 \
    -o stdout |
    hnsm filter -a 1000 stdin -o ${DIR_MERGE}/anchor.merge.fasta


#----------------------------#
# others
#----------------------------#
log_info "others"

cd ${BASH_DIR}/..

anchr contained \
    $( find . -path "*${DIR_PREFIX}*" -name "pe.others.fa" -or -path "*${DIR_PREFIX}*" -name "others.non-contained.fasta" | sort -r ) \
--len 500 --idt 0.9999 --ratio 0.99999 --parallel 24 \
    -o stdout |
    hnsm filter -a 500 stdin -o ${DIR_MERGE}/others.intermediate_0.fasta

anchr contained \
    ${DIR_MERGE}/anchor.merge.fasta \
    ${DIR_MERGE}/others.intermediate_0.fasta \
    --len 500 --idt 0.98 --ratio 0.99999 --parallel 24 \
    -o stdout |
    hnsm filter -a 500 stdin -o ${DIR_MERGE}/others.intermediate_1.fasta

cat ${DIR_MERGE}/others.intermediate_1.fasta |
    grep '>infile_1/' |
    sed 's/>//' \
    > ${DIR_MERGE}/others.txt

hnsm some -l 0 \
    ${DIR_MERGE}/others.intermediate_1.fasta \
    ${DIR_MERGE}/others.txt \
    -o ${DIR_MERGE}/others.non-contained.fasta

find ${DIR_MERGE} -name "anchor.intermediate*" | parallel --no-run-if-empty -j 1 rm
find ${DIR_MERGE} -name "others.intermediate*" | parallel --no-run-if-empty -j 1 rm

log_info Done.

exit 0


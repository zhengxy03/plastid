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
log_warn 7_fill_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 FILE_ANCHOR FILE_LONG GAP_COV"

if [ "$#" -lt 2 ]; then
    echo >&2 "$USAGE"
    exit 1
fi

FILE_ANCHOR=$1
FILE_LONG=$2
GAP_COV=${3:-2}

if [ -e 7_fill_anchors/contig.fasta ]; then
    echo >&2 "7_fill_anchors/contig.fasta presents"
    exit;
fi

#----------------------------#
# fill anchors
#使用长序列填充锚定序列之间的间隙，生成连续的序列（contigs）
#----------------------------#
mkdir -p 7_fill_anchors

log_info "overlap: between anchor-long"

anchr overlap2 \
    --parallel 24 \
    ${FILE_ANCHOR} \
    ${FILE_LONG} \
    -d 7_fill_anchors \
    -b 50 --len 1000 --idt 0.999 --all

cd 7_fill_anchors

CONTIG_COUNT=$(hnsm n50 -H -N 0 -C anchor.fasta)
log_debug "contigs: ${CONTIG_COUNT}"

log_info "group: anchor-long"
rm -fr group
dazz group \
    anchorLong.db \
    anchorLong.ovlp.tsv \
    --parallel 24 \
    --keep \
    --range "1-${CONTIG_COUNT}" --len 1000 --idt 0.999 --max 100 -c ${GAP_COV}

log_info "Processing each groups"
cat group/groups.txt |
    parallel --no-run-if-empty --linebuffer -k -j 12 '
        echo {};
        anchr orient \
            --len 1000 --idt 0.999 \
            group/{}.anchor.fasta \
            group/{}.long.fasta \
            -r group/{}.restrict.tsv \
            -o group/{}.strand.fasta;

        anchr overlap --len 1000 --idt 0.9999 --all \
            group/{}.strand.fasta \
            -o stdout |
            anchr restrict \
                stdin group/{}.restrict.tsv \
                -o group/{}.ovlp.tsv;

        dazz layout \
            group/{}.strand.fasta \
            group/{}.ovlp.tsv \
            group/{}.relation.tsv \
            -o group/{}.contig.fasta
    '

log_info "Build contigs"
cat \
   group/non_grouped.fasta \
   group/*.contig.fasta |
   hnsm filter -a 1000 stdin -o contig.fasta

log_info Done.

exit 0


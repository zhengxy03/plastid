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
log_warn 8_megahit.sh
#使用 Megahit 进行基因组组装，并对组装结果进行非包含序列筛选和锚定序列生成
#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 DIR_READS"

DIR_READS=${1:-"2_illumina/trim"}

# Convert to abs path
DIR_READS="$(cd "$(dirname "$DIR_READS")"; pwd)/$(basename "$DIR_READS")"

if [ -e 8_megahit/anchor/anchor.fasta ]; then
    log_info "8_megahit/anchor/anchor.fasta presents"
    exit;
fi

#----------------------------#
# spades
#----------------------------#
#MEGAHIT 是专门针对宏基因组数据优化的组装工具，能够高效处理海量短读长测序数据
if [ -e 8_megahit/megahit.non-contained.fasta ]; then
    log_info "8_megahit/megahit.non-contained.fasta presents"
else
    log_info "Run megahit"

    megahit \
        -t 8 \
        --k-list 31,41,51,61,71,81 \
        --12 ${DIR_READS}/pe.cor.fa.gz \
        --min-count 3 \
        -o 8_megahit

    anchr contained \
        8_megahit/final.contigs.fa \
        --len 1000 --idt 0.98 --ratio 0.99999 --parallel 16 \
        -o stdout |
        hnsm filter -a 1000 stdin -o 8_megahit/megahit.non-contained.fasta

    log_info "Clear intermediate files"
    find . -type d -path "*8_megahit/*" -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

mkdir -p 8_megahit/anchor
cd 8_megahit/anchor

anchr anchors \
    ../megahit.non-contained.fasta \
    ${DIR_READS}/pe.cor.fa.gz \
    --readl 125 \
    --uscale 2 \
    --lscale 3 \
    -p 8 \
    --ratio 0.98 \
    -o anchors.sh
bash anchors.sh

exit 0;


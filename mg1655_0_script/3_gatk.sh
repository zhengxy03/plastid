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
log_warn 3_gatk.sh

if [ ! -e 1_genome/genome.fa ]; then
    log_info "1_genome/genome.fa does not exist"
    exit;
fi

if [ ! -e 3_bwa/R.sort.bai ]; then
    log_info "3_bwa/R.sort.bai does not exist"
    exit;
fi

mkdir -p 3_gatk
cd 3_gatk

#----------------------------#
# Mutect2
#----------------------------#
#https://github.com/gatk-workflows/gatk4-mitochondria-pipeline/blob/master/tasks/align-and-call.wdl
#基因组变异检测
if [ ! -e R.raw.vcf ]; then
    gatk --java-options "-Xmx80g" \
        Mutect2 \
        --native-pair-hmm-threads 21 \
        -R ../3_bwa/genome.fa \
        -I ../3_bwa/R.sort.bam \
        --max-reads-per-alignment-start 100 \ #限制每个比对起始位置的最大 reads 数
        --max-mnp-distance 0 \  #禁用多核苷酸多态性（MNP）检测，仅检测单核苷酸变异（SNV）
        --annotation StrandBiasBySample \ #添加链偏倚注释，评估变异在正负链上的分布是否异常（可能提示测序误差）
        --mitochondria-mode \
        -O R.raw.vcf
fi

#----------------------------#
# Filter
#----------------------------#
if [ ! -e R.filtered.vcf ]; then
    gatk --java-options "-Xmx80g" \
        FilterMutectCalls \ #过滤低质量变异
        -R ../3_bwa/genome.fa \
        -V R.raw.vcf \
        --max-alt-allele-count 4 \
        --mitochondria-mode \
        -O R.filtered.vcf
fi

log_info Done.

exit 0


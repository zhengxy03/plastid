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
log_warn 3_bwa.sh

if [ ! -e 1_genome/genome.fa ]; then
    log_info "1_genome/genome.fa does not exist"
    exit;
fi

mkdir -p 3_bwa
cd 3_bwa

ln -fs ../1_genome/genome.fa genome.fa

# bwa index
if [ ! -e genome.fa.bwt ]; then
    bwa index genome.fa
fi

# faidx
if [ ! -e genome.fa.fai ]; then
    samtools faidx genome.fa
fi

# dict
if [ ! -e genome.dict ]; then
    picard CreateSequenceDictionary --REFERENCE genome.fa
fi

# chr.sizes
if [ ! -e chr.sizes ]; then
    hnsm size genome.fa > chr.sizes
fi

#----------------------------#
# Mapping
#----------------------------#
SAMPLE=$(readlinkf "${BASH_DIR}/.." | xargs dirname)
if [ ! -e R.sort.bai ]; then
    bwa mem -t 21 \
        -M -K 100000000 -v 3 -Y \
        genome.fa \
        ../2_illumina/Q25L60/R1.fq.gz \
        ../2_illumina/Q25L60/R2.fq.gz \
        2> >(tee R.bwa.log >&2) |
        samtools view -F 4 -u - | # Remove unmapped reads, write uncompressed BAM output
        picard AddOrReplaceReadGroups \
            --INPUT /dev/stdin \
            --OUTPUT /dev/stdout \
            --RGLB ${SAMPLE} \
            --RGPL ILLUMINA \
            --RGPU ${SAMPLE} \
            --RGSM ${SAMPLE} \
            --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 0 |
        picard CleanSam \ #移除无效的 reads
            --INPUT /dev/stdin \
            --OUTPUT /dev/stdout \
            --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 0 |
        picard FixMateInformation \  #修复配对信息
            --INPUT /dev/stdin \
            --OUTPUT R.mate.bam \
            --SORT_ORDER queryname \
            --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 1

    picard MarkDuplicates \
        --INPUT R.mate.bam \
        --OUTPUT /dev/stdout \
        --METRICS_FILE R.dedup.metrics \
        --ASSUME_SORT_ORDER queryname \
        --REMOVE_DUPLICATES true \
        --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 0 |
    picard SortSam \
        --INPUT /dev/stdin \
        --OUTPUT R.sort.bam \
        --SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY LENIENT --COMPRESSION_LEVEL 1

    picard BuildBamIndex \
        --INPUT R.sort.bam \
        --VALIDATION_STRINGENCY LENIENT

    picard CollectWgsMetrics \
        --INPUT R.sort.bam \
        --REFERENCE_SEQUENCE genome.fa \
        --OUTPUT R.wgs.metrics \
        --VALIDATION_STRINGENCY SILENT

    rm -f R.mate.bam
fi

#----------------------------#
# Depth
#----------------------------#
if [ ! -f R.mosdepth.summary.txt ]; then
    if hash mosdepth 2>/dev/null; then
        # depth
        mosdepth R R.sort.bam

        # covered
        gzip -dcf R.per-base.bed.gz |
            perl -nla -F"\t" -e '
                $F[3] == 0 and next;
                $start = $F[1] + 1;
                $end = $F[2];
                if ($start == $F[2]) {
                    print qq($F[0]:$start);
                }
                else {
                    print qq($F[0]:$start-$end);
                }
            ' |
            spanr cover stdin -o covered.json #基因组组装

        spanr stat chr.sizes covered.json -o stdout |
            grep -v "^all" |
            sed 's/^chr/chrom/' |
            sed 's/,size/,covLength/' |
            sed 's/,coverage/,covRate/' |
            sed 's/,/\t/g' \
            > coverage.tsv

        # join
        cat coverage.tsv |
            tsv-join -H --filter-file R.mosdepth.summary.txt \
                --key-fields chrom --append-fields 3-6 \
            > join.tsv

    fi
fi

log_info Done.

exit 0


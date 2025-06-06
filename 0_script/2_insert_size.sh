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
log_warn 2_insert_size.sh

mkdir -p 2_illumina/insert_size
cd 2_illumina/insert_size

for PREFIX in R S T; do
    if [ ! -e ../${PREFIX}1.fq.gz ]; then
        continue;
    fi

    if [ -e ${PREFIX}.ihist.tadpole.txt ]; then
        continue;
    fi

    tadpole.sh \ #BBTools 中的序列组装工具
        in=../${PREFIX}1.fq.gz \
        in2=../${PREFIX}2.fq.gz \
        out=${PREFIX}.tadpole.contig.fasta \
        threads=8 \
        overwrite 

    cat ${PREFIX}.tadpole.contig.fasta |
        anchr dazzname --no-replace --prefix T stdin \
        > ${PREFIX}.tadpole.contig.fa

    bbmap.sh \ #BBTools 中的短序列比对工具
        in=../${PREFIX}1.fq.gz \
        in2=../${PREFIX}2.fq.gz \
        out=${PREFIX}.tadpole.sam.gz \
        ref=${PREFIX}.tadpole.contig.fa \
        threads=8 \
        pairedonly \ #仅保留双端都比对上的 reads
        reads=1000000 \ #随机抽取1000000条reads
        nodisk overwrite

    reformat.sh \ #BBTools 中的格式转换工具
        in=${PREFIX}.tadpole.sam.gz \
        ihist=${PREFIX}.ihist.tadpole.txt \  #生成插入片段长度直方图
        overwrite

    picard SortSam \  #SortSam：按染色体坐标排序 SAM/BAM 文件
        -I ${PREFIX}.tadpole.sam.gz \
        -O ${PREFIX}.tadpole.sort.bam \
        --SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY LENIENT

    picard CollectInsertSizeMetrics \  #计算插入片段的统计信息（均值、标准差等）
        -I ${PREFIX}.tadpole.sort.bam \
        -O ${PREFIX}.insert_size.tadpole.txt \
        --Histogram_FILE ${PREFIX}.insert_size.tadpole.pdf

    if [ -e ../../1_genome/genome.fa ]; then
        bbmap.sh \
            in=../${PREFIX}1.fq.gz \
            in2=../${PREFIX}2.fq.gz \
            out=${PREFIX}.genome.sam.gz \
            ref=../../1_genome/genome.fa \
            threads=8 \
            maxindel=0 strictmaxindel \
            reads=1000000 \
            nodisk overwrite

        reformat.sh \
            in=${PREFIX}.genome.sam.gz \
            ihist=${PREFIX}.ihist.genome.txt \
            overwrite

        picard SortSam \
            -I ${PREFIX}.genome.sam.gz \
            -O ${PREFIX}.genome.sort.bam \
            --SORT_ORDER coordinate

        picard CollectInsertSizeMetrics \
            -I ${PREFIX}.genome.sort.bam \
            -O ${PREFIX}.insert_size.genome.txt \
            --Histogram_FILE ${PREFIX}.insert_size.genome.pdf
    fi

    find . -name "${PREFIX}.*.sam.gz" -or -name "${PREFIX}.*.sort.bam" |
        parallel --no-run-if-empty -j 1 rm
done

printf "%s\t%s\t%s\t%s\t%s\n" \
    "Group" "Mean" "Median" "STDev" "Pairs%/Orientation" \
    > statInsertSize.tsv

# bbtools reformat.sh
#Mean	339.868
#Median	312
#Mode	251
#STDev	134.676
#PercentOfPairs	36.247
for PREFIX in R S T; do
    for G in genome tadpole; do
        if [ ! -e ${PREFIX}.ihist.${G}.txt ]; then
            continue;
        fi

        cat ${PREFIX}.ihist.${G}.txt |
            GROUP="${PREFIX}.${G}" perl -nla -e '
                BEGIN { our $stat = { }; };

                m{\#(Mean|Median|STDev|PercentOfPairs)} or next;
                $stat->{$1} = $F[1];

                END {
                    printf qq(%s\t%.1f\t%s\t%.1f\t%.2f%%\n),
                        qq($ENV{GROUP}.bbtools),
                        $stat->{Mean},
                        $stat->{Median},
                        $stat->{STDev},
                        $stat->{PercentOfPairs};
                }
                '
    done
done \
    >> statInsertSize.tsv

# picard CollectInsertSizeMetrics
#MEDIAN_INSERT_SIZE	MODE_INSERT_SIZE	MEDIAN_ABSOLUTE_DEVIATION	MIN_INSERT_SIZE	MAX_INSERT_SIZE	MEAN_INSERT_SIZE	STANDARD_DEVIATION	READ_PAIRS	PAIR_ORIENTATION	WIDTH_OF_10_PERCENT	WIDTH_OF_20_PERCENT	WIDTH_OF_30_PERCENT	WIDTH_OF_40_PERCENT	WIDTH_OF_50_PERCENT	WIDTH_OF_60_PERCENT	WIDTH_OF_70_PERCENT	WIDTH_OF_80_PERCENT	WIDTH_OF_90_PERCENT	WIDTH_OF_95_PERCENT	WIDTH_OF_99_PERCENT	SAMPLE	LIBRARY	READ_GROUP
#296	287	14	92	501	294.892521	21.587526	1611331	FR	7	11	17	23	29	35	41	49	63	81	145
for PREFIX in R S T; do
    for G in genome tadpole; do
        if [ ! -e ${PREFIX}.insert_size.${G}.txt ]; then
            continue;
        fi

        cat ${PREFIX}.insert_size.${G}.txt |
            GROUP="${PREFIX}.${G}" perl -nla -F"\t" -e '
                next if @F < 9;
                next unless /^\d/;
                printf qq(%s\t%.1f\t%s\t%.1f\t%s\n),
                    qq($ENV{GROUP}.picard),
                    $F[5],
                    $F[0],
                    $F[6],
                    $F[8];
                '
    done
done \
    >> statInsertSize.tsv

cat statInsertSize.tsv |
    rgr md stdin --right 2-5 \
    > statInsertSize.md

echo -e "\nTable: statInsertSize\n" >> statInsertSize.md

cat statInsertSize.md
mv statInsertSize.md ${BASH_DIR}/../9_markdown


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
log_warn 9_busco.sh
ARRAY=()

if [ -e 1_genome/genome.fa ]; then
    ARRAY+=('Genome::1_genome/genome.fa')
fi
if [ -e 1_genome/paralogs.fa ]; then
    ARRAY+=('Paralogs::1_genome/paralogs.fa')
fi
if [ -e 1_genome/repetitive.fa ]; then
    ARRAY+=('repetitive::1_genome/repetitive.fa')
fi

if [ -e 7_merge_unitigs_bcalm/anchor.merge.fasta ]; then
    ARRAY+=('merge_bcalm::7_merge_unitigs_bcalm/anchor.merge.fasta')
fi
if [ -e 7_merge_unitigs_bifrost/anchor.merge.fasta ]; then
    ARRAY+=('merge_bifrost::7_merge_unitigs_bifrost/anchor.merge.fasta')
fi
if [ -e 7_merge_unitigs_superreads/anchor.merge.fasta ]; then
    ARRAY+=('merge_superreads::7_merge_unitigs_superreads/anchor.merge.fasta')
fi

if [ -e 7_merge_mr_unitigs_bcalm/anchor.merge.fasta ]; then
    ARRAY+=('merge_mr_bcalm::7_merge_mr_unitigs_bcalm/anchor.merge.fasta')
fi
if [ -e 7_merge_mr_unitigs_bifrost/anchor.merge.fasta ]; then
    ARRAY+=('merge_mr_bifrost::7_merge_mr_unitigs_bifrost/anchor.merge.fasta')
fi
if [ -e 7_merge_mr_unitigs_superreads/anchor.merge.fasta ]; then
    ARRAY+=('merge_mr_superreads::7_merge_mr_unitigs_superreads/anchor.merge.fasta')
fi

if [ -e 7_merge_anchors/anchor.merge.fasta ]; then
    ARRAY+=('merge_anchors::7_merge_anchors/anchor.merge.fasta')
fi

if [ -e 7_glue_anchors/contig.fasta ]; then
    ARRAY+=('glue_anchors::7_glue_anchors/contig.fasta')
fi
if [ -e 7_fill_anchors/contig.fasta ]; then
    ARRAY+=('fill_anchors::7_fill_anchors/contig.fasta')
fi

if [ -e 8_spades/spades.non-contained.fasta ]; then
    ARRAY+=('spades::8_spades/spades.non-contained.fasta')
fi
if [ -e 8_mr_spades/spades.non-contained.fasta ]; then
    ARRAY+=('mr_spades::8_mr_spades/spades.non-contained.fasta')
fi

if [ -e 8_megahit/megahit.non-contained.fasta ]; then
    ARRAY+=('megahit::8_megahit/megahit.non-contained.fasta')
fi
if [ -e 8_mr_megahit/megahit.non-contained.fasta ]; then
    ARRAY+=('mr_megahit::8_mr_megahit/megahit.non-contained.fasta')
fi

mkdir -p 9_busco

for item in "${ARRAY[@]}" ; do
    NAME="${item%%::*}"
    FILE="${item##*::}"

    echo "==> ${NAME}"

    cat ${FILE} |
        sed "s/\//_/g" > tmp.fasta
    #通过比对一组保守的单拷贝直系同源基因（universal single-copy orthologs）来量化评估组装结果的完整性 
    busco \
        -i tmp.fasta \
        -o ${NAME} \
        --out_path 9_busco \
        --auto-lineage-prok \ #自动选择原核生物谱系数据集
        -m genome --cpu 24

    rm tmp.fasta

    # Clean unused directories
    find 9_busco/${NAME} -maxdepth 1 -mindepth 1 -type d |
        grep -v "run_" |
        xargs rm -fr

    find 9_busco/${NAME} -maxdepth 2 -mindepth 2 -type d |
        xargs rm -fr

done

# Save lineages
find 9_busco/ -maxdepth 2 -mindepth 2 -type d |
    parallel -j 1 'echo {/}' |
    grep "run_" |
    sort -u \
    > 9_busco/lineages.txt

for L in $(cat 9_busco/lineages.txt); do
    echo "Table: statBusco ${L}"
    echo
    echo '| NAME | C | S | D | F | M | Total |'
    echo '|:--|--:|--:|--:|--:|--:|--:|'

    for item in "${ARRAY[@]}" ; do
        NAME="${item%%::*}"

        if [ -f 9_busco/${NAME}/${L}/short_summary.txt ]; then
            cat 9_busco/${NAME}/${L}/short_summary.txt |
                NAME=${NAME} perl -n -MJSON::PP -e '
                    if (m/(\d+)\s*Complete BUSCOs/){
                        $C = $1;
                    }
                    if (m/(\d+)\s*Complete and single-copy BUSCOs/){
                        $S = $1;
                    }
                    if (m/(\d+)\s*Complete and duplicated BUSCOs/){
                        $D = $1;
                    }
                    if (m/(\d+)\s*Fragmented BUSCOs/){
                        $F = $1;
                    }
                    if (m/(\d+)\s*Missing BUSCOs/){
                        $M = $1;
                    }
                    if (m/(\d+)\s*Total BUSCO groups searched/){
                        $T = $1;
                    }
                    END{
                        printf "| %s | %d | %d | %d | %d | %d | %d |\n",
                            $ENV{NAME},
                            $C || 0,
                            $S || 0,
                            $D || 0,
                            $F || 0,
                            $M || 0,
                            $T || 0;
                    }
                '
        fi
    done

    echo
    echo

done \
    > statBusco.md

find . -type f -name "busco*log" | xargs rm


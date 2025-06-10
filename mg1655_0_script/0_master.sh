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
log_warn 0_master.sh

#----------------------------#
# Illumina QC
#----------------------------#
if [ -e 0_script/2_fastqc.sh ]; then
    bash 0_script/2_fastqc.sh
fi

if [ -e 0_script/2_insert_size.sh ]; then
    bash 0_script/2_insert_size.sh
fi

if [ -e 0_script/2_fastk.sh ]; then
    bash 0_script/2_fastk.sh
fi

#----------------------------#
# trim reads
#----------------------------#
if [ -e 0_script/2_trim.sh ]; then
    bash 0_script/2_trim.sh
fi

if [ -e 0_script/9_stat_reads.sh ]; then
    bash 0_script/9_stat_reads.sh
fi

#----------------------------#
# merge reads
#----------------------------#
if [ -e 0_script/2_merge.sh ]; then
    bash 0_script/2_merge.sh
fi

#----------------------------#
# quorum
#----------------------------#
if [ -e 0_script/2_quorum.sh ]; then
    bash 0_script/2_quorum.sh
fi

#----------------------------#
# mapping
#----------------------------#
if [ -e 0_script/3_bwa.sh ]; then
    bash 0_script/3_bwa.sh
fi
if [ -e 0_script/3_gatk.sh ]; then
    bash 0_script/3_gatk.sh
fi

#----------------------------#
# down sampling trimmed reads; build unitigs and anchors
#----------------------------#
if [ -e 0_script/4_down_sampling.sh ]; then
    bash 0_script/4_down_sampling.sh
fi

if [ -e 0_script/4_unitigs_bcalm.sh ]; then
    bash 0_script/4_unitigs_bcalm.sh
fi
if [ -e 0_script/4_anchors.sh ]; then
    bash 0_script/4_anchors.sh 4_unitigs_bcalm
fi
if [ -e 0_script/9_stat_anchors.sh ]; then
    bash 0_script/9_stat_anchors.sh 4_unitigs_bcalm statUnitigsBcalm.md
fi
if [ -e 0_script/4_unitigs_bifrost.sh ]; then
    bash 0_script/4_unitigs_bifrost.sh
fi
if [ -e 0_script/4_anchors.sh ]; then
    bash 0_script/4_anchors.sh 4_unitigs_bifrost
fi
if [ -e 0_script/9_stat_anchors.sh ]; then
    bash 0_script/9_stat_anchors.sh 4_unitigs_bifrost statUnitigsBifrost.md
fi
if [ -e 0_script/4_unitigs_superreads.sh ]; then
    bash 0_script/4_unitigs_superreads.sh
fi
if [ -e 0_script/4_anchors.sh ]; then
    bash 0_script/4_anchors.sh 4_unitigs_superreads
fi
if [ -e 0_script/9_stat_anchors.sh ]; then
    bash 0_script/9_stat_anchors.sh 4_unitigs_superreads statUnitigsSuperreads.md
fi

#----------------------------#
# down sampling merged reads
#----------------------------#
if [ -e 0_script/6_down_sampling.sh ]; then
    bash 0_script/6_down_sampling.sh
fi

if [ -e 0_script/6_unitigs_bcalm.sh ]; then
    bash 0_script/6_unitigs_bcalm.sh
fi
if [ -e 0_script/6_anchors.sh ]; then
    bash 0_script/6_anchors.sh 6_unitigs_bcalm
fi
if [ -e 0_script/9_stat_anchors.sh ]; then
    bash 0_script/9_stat_mr_anchors.sh 6_unitigs_bcalm statMRUnitigsBcalm.md
fi
if [ -e 0_script/6_unitigs_bifrost.sh ]; then
    bash 0_script/6_unitigs_bifrost.sh
fi
if [ -e 0_script/6_anchors.sh ]; then
    bash 0_script/6_anchors.sh 6_unitigs_bifrost
fi
if [ -e 0_script/9_stat_anchors.sh ]; then
    bash 0_script/9_stat_mr_anchors.sh 6_unitigs_bifrost statMRUnitigsBifrost.md
fi
if [ -e 0_script/6_unitigs_superreads.sh ]; then
    bash 0_script/6_unitigs_superreads.sh
fi
if [ -e 0_script/6_anchors.sh ]; then
    bash 0_script/6_anchors.sh 6_unitigs_superreads
fi
if [ -e 0_script/9_stat_anchors.sh ]; then
    bash 0_script/9_stat_mr_anchors.sh 6_unitigs_superreads statMRUnitigsSuperreads.md
fi

#----------------------------#
# merge anchors
#----------------------------#
if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 4_unitigs_bcalm 7_merge_unitigs_bcalm
fi
if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 4_unitigs_bifrost 7_merge_unitigs_bifrost
fi
if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 4_unitigs_superreads 7_merge_unitigs_superreads
fi

if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 6_unitigs_bcalm 7_merge_mr_unitigs_bcalm
fi
if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 6_unitigs_bifrost 7_merge_mr_unitigs_bifrost
fi
if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 6_unitigs_superreads 7_merge_mr_unitigs_superreads
fi

if [ -e 0_script/7_merge_anchors.sh ]; then
    bash 0_script/7_merge_anchors.sh 7_merge 7_merge_anchors
fi

if [ -e 0_script/9_stat_merge_anchors.sh ]; then
    bash 0_script/9_stat_merge_anchors.sh
fi

#----------------------------#
# spades, megahit
#----------------------------#
if [ -e 0_script/8_spades.sh ]; then
    bash 0_script/8_spades.sh
fi
if [ -e 0_script/8_mr_spades.sh ]; then
    bash 0_script/8_mr_spades.sh
fi
if [ -e 0_script/8_megahit.sh ]; then
    bash 0_script/8_megahit.sh
fi
if [ -e 0_script/8_mr_megahit.sh ]; then
    bash 0_script/8_mr_megahit.sh
fi

if [ -e 0_script/9_stat_other_anchors.sh ]; then
    bash 0_script/9_stat_other_anchors.sh
fi

#----------------------------#
# extend anchors
#----------------------------#
rm -fr 7_extend_anchors
mkdir -p 7_extend_anchors
cat \
    8_spades/spades.non-contained.fasta \
    8_megahit/megahit.non-contained.fasta \
8_mr_spades/spades.non-contained.fasta \
    8_mr_megahit/megahit.non-contained.fasta \
| anchr dazzname --no-replace stdin \
    | hnsm filter -a 1000 stdin -o 7_extend_anchors/contigs.2GS.fasta

if [ -e 0_script/7_glue_anchors.sh ]; then
    bash 0_script/7_glue_anchors.sh 7_merge_anchors/anchor.merge.fasta 7_extend_anchors/contigs.2GS.fasta 3
fi
if [ -e 0_script/7_fill_anchors.sh ]; then
    bash 0_script/7_fill_anchors.sh 7_glue_anchors/contig.fasta 7_extend_anchors/contigs.2GS.fasta 3
fi

#----------------------------#
# final stats
#----------------------------#
if [ -e 0_script/9_stat_final.sh ]; then
    bash 0_script/9_stat_final.sh
fi
if [ -e 0_script/9_quast.sh ]; then
    bash 0_script/9_quast.sh
fi


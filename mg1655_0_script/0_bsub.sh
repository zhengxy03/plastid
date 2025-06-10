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
log_warn 0_bsub.sh

#----------------------------#
# Illumina QC
#----------------------------#
if [ -e 0_script/2_fastqc.sh ]; then
    bsub -q mpi -n 24 -J "${BASE_NAME}-2_fastqc" \
        "bash 0_script/2_fastqc.sh"
fi

if [ -e 0_script/2_insert_size.sh ]; then
    bsub -q mpi -n 24 -J "${BASE_NAME}-2_insert_size" \
        "bash 0_script/2_insert_size.sh"
fi

if [ -e 0_script/2_fastk.sh ]; then
    bsub -q mpi -n 24 -J "${BASE_NAME}-2_fastk" \
        "bash 0_script/2_fastk.sh"
fi

#----------------------------#
# trim reads
#----------------------------#
bsub -q mpi -n 24 -J "${BASE_NAME}-2_trim" \
    "bash 0_script/2_trim.sh"

bsub -w "ended(${BASE_NAME}-2_trim)" \
    -q mpi -n 24 -J "${BASE_NAME}-9_stat_reads" \
    "bash 0_script/9_stat_reads.sh"

if [ -e 0_script/3_bwa.sh ]; then
    bsub  -w "ended(${BASE_NAME}-2_trim)" \
        -q mpi -n 24 -J "${BASE_NAME}-3_bwa" \
        "bash 0_script/3_bwa.sh"
fi

if [ -e 0_script/3_gatk.sh ]; then
    bsub  -w "ended(${BASE_NAME}-3_bwa)" \
        -q mpi -n 24 -J "${BASE_NAME}-3_gatk" \
        "bash 0_script/3_gatk.sh"
fi

#----------------------------#
# merge reads
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_trim)" \
    -q mpi -n 24 -J "${BASE_NAME}-2_merge" \
    "bash 0_script/2_merge.sh"

#----------------------------#
# quorum
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_trim)" \
    -q mpi -n 24 -J "${BASE_NAME}-2_quorum" \
    "bash 0_script/2_quorum.sh"

#----------------------------#
# down sampling trimmed reads; build unitigs and anchors
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_quorum)" \
    -q mpi -n 24 -J "${BASE_NAME}-4_down_sampling" \
    "bash 0_script/4_down_sampling.sh"

bsub -w "ended(${BASE_NAME}-4_down_sampling)" \
    -q mpi -n 24 -J "${BASE_NAME}-4_unitigs_bcalm" \
    "
    bash 0_script/4_unitigs_bcalm.sh
    bash 0_script/4_anchors.sh 4_unitigs_bcalm
    bash 0_script/9_stat_anchors.sh 4_unitigs_bcalm statUnitigsBcalm.md
    "

bsub -w "ended(${BASE_NAME}-4_down_sampling)" \
    -q mpi -n 24 -J "${BASE_NAME}-4_unitigs_bifrost" \
    "
    bash 0_script/4_unitigs_bifrost.sh
    bash 0_script/4_anchors.sh 4_unitigs_bifrost
    bash 0_script/9_stat_anchors.sh 4_unitigs_bifrost statUnitigsBifrost.md
    "

bsub -w "ended(${BASE_NAME}-4_down_sampling)" \
    -q mpi -n 24 -J "${BASE_NAME}-4_unitigs_superreads" \
    "
    bash 0_script/4_unitigs_superreads.sh
    bash 0_script/4_anchors.sh 4_unitigs_superreads
    bash 0_script/9_stat_anchors.sh 4_unitigs_superreads statUnitigsSuperreads.md
    "


#----------------------------#
# down sampling merged reads
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_merge)" \
    -q mpi -n 24 -J "${BASE_NAME}-6_down_sampling" \
    "bash 0_script/6_down_sampling.sh"

bsub -w "ended(${BASE_NAME}-6_down_sampling)" \
    -q mpi -n 24 -J "${BASE_NAME}-6_unitigs_bcalm" \
    "
    bash 0_script/6_unitigs_bcalm.sh
    bash 0_script/6_anchors.sh 6_unitigs_bcalm
    bash 0_script/9_stat_mr_anchors.sh 6_unitigs_bcalm statMRUnitigsBcalm.md
    "

bsub -w "ended(${BASE_NAME}-6_down_sampling)" \
    -q mpi -n 24 -J "${BASE_NAME}-6_unitigs_bifrost" \
    "
    bash 0_script/6_unitigs_bifrost.sh
    bash 0_script/6_anchors.sh 6_unitigs_bifrost
    bash 0_script/9_stat_mr_anchors.sh 6_unitigs_bifrost statMRUnitigsBifrost.md
    "

bsub -w "ended(${BASE_NAME}-6_down_sampling)" \
    -q mpi -n 24 -J "${BASE_NAME}-6_unitigs_superreads" \
    "
    bash 0_script/6_unitigs_superreads.sh
    bash 0_script/6_anchors.sh 6_unitigs_superreads
    bash 0_script/9_stat_mr_anchors.sh 6_unitigs_superreads statMRUnitigsSuperreads.md
    "


#----------------------------#
# merge anchors
#----------------------------#
bsub -w "ended(${BASE_NAME}-4_unitigs_bcalm)" \
    -q mpi -n 24 -J "${BASE_NAME}-7_merge_anchors_4_unitigs_bcalm" \
    "bash 0_script/7_merge_anchors.sh 4_unitigs_bcalm 7_merge_unitigs_bcalm"
bsub -w "ended(${BASE_NAME}-4_unitigs_bifrost)" \
    -q mpi -n 24 -J "${BASE_NAME}-7_merge_anchors_4_unitigs_bifrost" \
    "bash 0_script/7_merge_anchors.sh 4_unitigs_bifrost 7_merge_unitigs_bifrost"
bsub -w "ended(${BASE_NAME}-4_unitigs_superreads)" \
    -q mpi -n 24 -J "${BASE_NAME}-7_merge_anchors_4_unitigs_superreads" \
    "bash 0_script/7_merge_anchors.sh 4_unitigs_superreads 7_merge_unitigs_superreads"

bsub -w "ended(${BASE_NAME}-6_unitigs_bcalm)" \
    -q mpi -n 24 -J "${BASE_NAME}-7_merge_anchors_6_unitigs_bcalm" \
    "bash 0_script/7_merge_anchors.sh 6_unitigs_bcalm 7_merge_mr_unitigs_bcalm"
bsub -w "ended(${BASE_NAME}-6_unitigs_bifrost)" \
    -q mpi -n 24 -J "${BASE_NAME}-7_merge_anchors_6_unitigs_bifrost" \
    "bash 0_script/7_merge_anchors.sh 6_unitigs_bifrost 7_merge_mr_unitigs_bifrost"
bsub -w "ended(${BASE_NAME}-6_unitigs_superreads)" \
    -q mpi -n 24 -J "${BASE_NAME}-7_merge_anchors_6_unitigs_superreads" \
    "bash 0_script/7_merge_anchors.sh 6_unitigs_superreads 7_merge_mr_unitigs_superreads"

bsub -w "ended(${BASE_NAME}-2_quorum) && ended(${BASE_NAME}-7_merge_anchors_4_unitigs_bcalm)&& ended(${BASE_NAME}-7_merge_anchors_4_unitigs_bifrost)&& ended(${BASE_NAME}-7_merge_anchors_4_unitigs_superreads) && ended(${BASE_NAME}-7_merge_anchors_4_unitigs_bcalm)&& ended(${BASE_NAME}-7_merge_anchors_4_unitigs_bifrost)&& ended(${BASE_NAME}-7_merge_anchors_4_unitigs_superreads)" \
    -q mpi -n 24 -J "${BASE_NAME}-7_merge_anchors" \
    "bash 0_script/7_merge_anchors.sh 7_merge 7_merge_anchors"
bsub -w "ended(${BASE_NAME}-7_merge_anchors)" \
    -q mpi -n 24 -J "${BASE_NAME}-9_stat_merge_anchors" \
    "bash 0_script/9_stat_merge_anchors.sh"

#----------------------------#
# spades, megahit
#----------------------------#
bsub -w "ended(${BASE_NAME}-2_quorum)" \
    -q mpi -n 24 -J "${BASE_NAME}-8_spades" \
    "bash 0_script/8_spades.sh"

bsub -w "ended(${BASE_NAME}-2_quorum)" \
    -q mpi -n 24 -J "${BASE_NAME}-8_megahit" \
    "bash 0_script/8_megahit.sh"

bsub -w "ended(${BASE_NAME}-2_merge)" \
    -q mpi -n 24 -J "${BASE_NAME}-8_mr_spades" \
    "bash 0_script/8_mr_spades.sh"
bsub -w "ended(${BASE_NAME}-2_merge)" \
    -q mpi -n 24 -J "${BASE_NAME}-8_mr_megahit" \
    "bash 0_script/8_mr_megahit.sh"

bsub -w "ended(${BASE_NAME}-8_spades) && ended(${BASE_NAME}-8_megahit) && ended(${BASE_NAME}-8_mr_spades) && ended(${BASE_NAME}-8_mr_megahit)" \
    -q mpi -n 24 -J "${BASE_NAME}-9_stat_other_anchors" \
    "bash 0_script/9_stat_other_anchors.sh"

#----------------------------#
# extend anchors
#----------------------------#
bsub -w "ended(${BASE_NAME}-8_spades) && ended(${BASE_NAME}-8_megahit) && ended(${BASE_NAME}-8_mr_spades)&& ended(${BASE_NAME}-8_mr_megahit)" \
    -q mpi -n 24 -J "${BASE_NAME}-contigs_2GS" \
    '
    rm -fr 7_extend_anchors
    mkdir -p 7_extend_anchors
    cat \
        8_spades/spades.non-contained.fasta \
        8_megahit/megahit.non-contained.fasta \
8_mr_spades/spades.non-contained.fasta \
        8_mr_megahit/megahit.non-contained.fasta \
| anchr dazzname --no-replace stdin \
        | hnsm filter -a 1000 stdin -o 7_extend_anchors/contigs.2GS.fasta
    '

bsub -w "ended(${BASE_NAME}-7_merge_anchors) && ended(${BASE_NAME}-contigs_2GS)" \
    -q mpi -n 24 -J "${BASE_NAME}-7_glue_anchors" \
    "bash 0_script/7_glue_anchors.sh 7_merge_anchors/anchor.merge.fasta 7_extend_anchors/contigs.2GS.fasta 3"

bsub -w "ended(${BASE_NAME}-7_glue_anchors)" \
    -q mpi -n 24 -J "${BASE_NAME}-7_fill_anchors" \
    "bash 0_script/7_fill_anchors.sh 7_glue_anchors/contig.fasta 7_extend_anchors/contigs.2GS.fasta 3"

#----------------------------#
# final stats
#----------------------------#
bsub -w "ended(${BASE_NAME}-7_merge_anchors) && ended(${BASE_NAME}-8_spades) && ended(${BASE_NAME}-7_fill_anchors)" \
    -q mpi -n 24 -J "${BASE_NAME}-9_stat_final" \
    "bash 0_script/9_stat_final.sh"

bsub -w "ended(${BASE_NAME}-7_merge_anchors) && ended(${BASE_NAME}-8_spades) && ended(${BASE_NAME}-7_fill_anchors)" \
    -q mpi -n 24 -J "${BASE_NAME}-9_quast" \
    "bash 0_script/9_quast.sh"

bsub -w "ended(${BASE_NAME}-9_stat_final) && ended(${BASE_NAME}-9_quast)" \
    -q mpi -n 24 -J "${BASE_NAME}-0_cleanup" \
    "bash 0_script/0_cleanup.sh"


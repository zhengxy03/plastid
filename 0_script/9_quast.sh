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
log_warn 9_quast.sh
QUAST_TARGET=
QUAST_LABEL=

if [ -e 1_genome/genome.fa ]; then
    QUAST_TARGET+=" -R 1_genome/genome.fa "
fi

if [ -e 7_merge_unitigs_bcalm/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_unitigs_bcalm/anchor.merge.fasta "
    QUAST_LABEL+="merge_bcalm,"
fi
if [ -e 7_merge_unitigs_bifrost/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_unitigs_bifrost/anchor.merge.fasta "
    QUAST_LABEL+="merge_bifrost,"
fi
if [ -e 7_merge_unitigs_superreads/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_unitigs_superreads/anchor.merge.fasta "
    QUAST_LABEL+="merge_superreads,"
fi
if [ -e 7_merge_unitigs_tadpole/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_unitigs_tadpole/anchor.merge.fasta "
    QUAST_LABEL+="merge_tadpole,"
fi

if [ -e 7_merge_mr_unitigs_bcalm/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_mr_unitigs_bcalm/anchor.merge.fasta "
    QUAST_LABEL+="merge_mr_bcalm,"
fi
if [ -e 7_merge_mr_unitigs_bifrost/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_mr_unitigs_bifrost/anchor.merge.fasta "
    QUAST_LABEL+="merge_mr_bifrost,"
fi
if [ -e 7_merge_mr_unitigs_superreads/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_mr_unitigs_superreads/anchor.merge.fasta "
    QUAST_LABEL+="merge_mr_superreads,"
fi
if [ -e 7_merge_mr_unitigs_tadpole/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_mr_unitigs_tadpole/anchor.merge.fasta "
    QUAST_LABEL+="merge_mr_tadpole,"
fi

if [ -e 7_merge_anchors/anchor.merge.fasta ]; then
    QUAST_TARGET+=" 7_merge_anchors/anchor.merge.fasta "
    QUAST_LABEL+="merge_anchors,"
fi

if [ -e 7_glue_anchors/contig.fasta ]; then
    QUAST_TARGET+=" 7_glue_anchors/contig.fasta "
    QUAST_LABEL+="glue_anchors,"
fi
if [ -e 7_fill_anchors/contig.fasta ]; then
    QUAST_TARGET+=" 7_fill_anchors/contig.fasta "
    QUAST_LABEL+="fill_anchors,"
fi

if [ -e 8_spades/spades.non-contained.fasta ]; then
    QUAST_TARGET+=" 8_spades/spades.non-contained.fasta "
    QUAST_LABEL+="spades,"
fi
if [ -e 8_mr_spades/spades.non-contained.fasta ]; then
    QUAST_TARGET+=" 8_mr_spades/spades.non-contained.fasta "
    QUAST_LABEL+="mr_spades,"
fi

if [ -e 8_megahit/megahit.non-contained.fasta ]; then
    QUAST_TARGET+=" 8_megahit/megahit.non-contained.fasta "
    QUAST_LABEL+="megahit,"
fi
if [ -e 8_mr_megahit/megahit.non-contained.fasta ]; then
    QUAST_TARGET+=" 8_mr_megahit/megahit.non-contained.fasta "
    QUAST_LABEL+="mr_megahit,"
fi

if [ -e 1_genome/paralogs.fa ]; then
    QUAST_TARGET+=" 1_genome/paralogs.fa "
    QUAST_LABEL+="paralogs,"
fi

if [ -e 1_genome/repetitive/repetitive.fa ]; then
    QUAST_TARGET+=" 1_genome/repetitive/repetitive.fa "
    QUAST_LABEL+="repetitive,"
fi

QUAST_LABEL=$( echo "${QUAST_LABEL}" | sed 's/,$//' )

rm -fr 9_quast
quast.py --threads 8 --min-contig 100 \
    ${QUAST_TARGET} \
    --label ${QUAST_LABEL} \
    -o 9_quast


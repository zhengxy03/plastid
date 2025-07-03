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
log_warn 9_statFinal.sh

tempfile=$(mktemp /tmp/stat_final_XXXXXXXX)
trap 'rm -f "$tempfile"' EXIT

printf "%s\t" \
    "Name" "N50" "Sum" "#" |
    sed 's/\t$/\n/' \
    > ${tempfile}

# genome
if [ -e 1_genome/genome.fa ]; then
    printf "%s\t" "Genome" $( stat_format 1_genome/genome.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 1_genome/paralogs.fa ]; then
    printf "%s\t" "Paralogs" $( stat_format 1_genome/paralogs.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 1_genome/repetitive/repetitive.fa ]; then
    printf "%s\t" "repetitive" $( stat_format 1_genome/repetitive/repetitive.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# anchors
for D in 7_merge_anchors; do
    if [ -e ${D}/anchor.merge.fasta ]; then
        printf "%s\t" "${D}.anchors" $( stat_format ${D}/anchor.merge.fasta ) |
            sed 's/\t$/\n/' \
            >> ${tempfile}
    fi
    if [ -e ${D}/others.non-contained.fasta ]; then
        printf "%s\t" "${D}.others" $( stat_format ${D}/others.non-contained.fasta ) |
            sed 's/\t$/\n/' \
            >> ${tempfile}
    fi
done

# extended anchors
if [ -e 7_glue_anchors/contig.fasta ]; then
    printf "%s\t" "glue_anchors" $( stat_format 7_glue_anchors/contig.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 7_fill_anchors/contig.fasta ]; then
    printf "%s\t" "fill_anchors" $( stat_format 7_fill_anchors/contig.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# spades
if [ -e 8_spades/contigs.fasta ]; then
    printf "%s\t" "spades.contig" $( stat_format 8_spades/contigs.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_spades/scaffolds.fasta ]; then
    printf "%s\t" "spades.scaffold" $( stat_format 8_spades/scaffolds.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_spades/spades.non-contained.fasta ]; then
    printf "%s\t" "spades.non-contained" $( stat_format 8_spades/spades.non-contained.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# mr_spades
if [ -e 8_mr_spades/contigs.fasta ]; then
    printf "%s\t" "mr_spades.contig" $( stat_format 8_mr_spades/contigs.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_mr_spades/scaffolds.fasta ]; then
    printf "%s\t" "mr_spades.scaffold" $( stat_format 8_mr_spades/scaffolds.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_mr_spades/spades.non-contained.fasta ]; then
    printf "%s\t" "mr_spades.non-contained" $( stat_format 8_mr_spades/spades.non-contained.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# megahit
if [ -e 8_megahit/final.contigs.fa ]; then
    printf "%s\t" "megahit.contig" $( stat_format 8_megahit/final.contigs.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_megahit/megahit.non-contained.fasta ]; then
    printf "%s\t" "megahit.non-contained" $( stat_format 8_megahit/megahit.non-contained.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

# mr_megahit
if [ -e 8_mr_megahit/final.contigs.fa ]; then
    printf "%s\t" "mr_megahit.contig" $( stat_format 8_mr_megahit/final.contigs.fa ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi
if [ -e 8_mr_megahit/megahit.non-contained.fasta ]; then
    printf "%s\t" "mr_megahit.non-contained" $( stat_format 8_mr_megahit/megahit.non-contained.fasta ) |
        sed 's/\t$/\n/' \
        >> ${tempfile}
fi

rgr md ${tempfile} --right 2-4 -o statFinal.md
echo -e "\nTable: statFinal\n" >> statFinal.md

cat statFinal.md
mv statFinal.md ${BASH_DIR}/../9_markdown


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
log_warn 8_mr_spades.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0"

if [ -e 8_mr_spades/anchor/anchor.fasta ]; then
    log_info "8_mr_spades/anchor/anchor.fasta presents"
    exit;
fi

#----------------------------#
# spades
#----------------------------#
if [ -e 8_mr_spades/spades.non-contained.fasta ]; then
    log_info "8_mr_spades/spades.non-contained.fasta presents"
else
    log_info "Run spades"

    mkdir -p 8_mr_spades
    cd 8_mr_spades

    mkdir -p re-pair
    faops filter -l 0 -a 60 ${BASH_DIR}/../2_illumina/merge/pe.cor.fa.gz stdout |
        repair.sh \
            in=stdin.fa \
            out=re-pair/R1.fa \
            out2=re-pair/R2.fa \
            outs=re-pair/Rs.fa \
            threads=24 \
            fint overwrite

    # spades seems ignore non-properly paired reads
    spades.py \
        -t 24 \
        --only-assembler \
        -k 25,55,95,125 \
        -1 re-pair/R1.fa \
        -2 re-pair/R2.fa \
        -s re-pair/Rs.fa \
        -o .

    anchr contained \
        contigs.fasta \
        --len 1000 --idt 0.98 --ratio 0.99999 --parallel 24 \
        -o stdout |
        faops filter -a 1000 -l 0 stdin spades.non-contained.fasta

    log_info "Clear intermediate files"
    find . -type d -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr
fi

#----------------------------#
# anchor
#----------------------------#
log_info "Create anchors"

cd ${BASH_DIR}/..
mkdir -p 8_mr_spades/anchor
cd 8_mr_spades/anchor

anchr anchors \
    ../spades.non-contained.fasta \
    ${BASH_DIR}/../2_illumina/merge/pe.cor.fa.gz \
    --readl 151 \
    --uscale 2 \
    --lscale 3 \
    -p 24 \
    --ratio 0.98 \
    -o anchors.sh
bash anchors.sh

exit 0;


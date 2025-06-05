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
log_warn 9_stat_reads.sh

cd 2_illumina

if [ -e statReads.tsv ]; then
    log_debug "statReads.tsv presents"
else

printf "%s\t%s\t%s\t%s\n" \
    "Name" "N50" "Sum" "#" \
    > statReads.tsv

for NAME in genome paralogs; do
    if [ -e ../1_genome/${NAME}.fa ]; then
        printf "%s\t%s\t%s\t%s\n" \
            $(echo "${NAME}"; stat_format ../1_genome/${NAME}.fa;)
    fi
done \
    >> statReads.tsv

if [ -e ../1_genome/repetitive/repetitive.fa ]; then
    printf "%s\t%s\t%s\t%s\n" \
        $(echo "repetitive"; stat_format ../1_genome/repetitive/repetitive.fa;)
fi \
    >> statReads.tsv

for PREFIX in R S T; do
    if [ -e ${PREFIX}1.fq.gz ]; then
        printf "%s\t%s\t%s\t%s\n" \
            $(echo "Illumina.${PREFIX}"; stat_format ${PREFIX}1.fq.gz ${PREFIX}2.fq.gz;)
    fi
    if [ -e trim/${PREFIX}1.fq.gz ]; then
        printf "%s\t%s\t%s\t%s\n" \
            $(echo "trim.${PREFIX}"; stat_format trim/${PREFIX}1.fq.gz trim/${PREFIX}2.fq.gz trim/${PREFIX}s.fq.gz;)
    fi
done \
    >> statReads.tsv

for PREFIX in R S T; do
    for Q in 0 25 30; do
        for L in 0 60; do
            if [ ! -e Q${Q}L${L}/${PREFIX}1.fq.gz ]; then
                continue
            fi

            printf "%s\t%s\t%s\t%s\n" \
                $(
                    echo Q${Q}L${L};

                    stat_format \
                        Q${Q}L${L}/${PREFIX}1.fq.gz \
                        Q${Q}L${L}/${PREFIX}2.fq.gz \
                        Q${Q}L${L}/${PREFIX}s.fq.gz;

                )
        done
    done
done \
    >> statReads.tsv

fi # end of statReads

rgr md statReads.tsv --right 2-4 -o statReads.md
echo -e "\nTable: statReads\n" >> statReads.md

cat statReads.md
mv statReads.md ${BASH_DIR}/../9_markdown


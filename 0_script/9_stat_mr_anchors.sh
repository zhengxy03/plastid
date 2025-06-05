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
log_warn 9_stat_mr_anchors.sh

#----------------------------#
# set parameters
#----------------------------#
USAGE="Usage: $0 [DIR_PREFIX] [FILENAME_MD]"

DIR_PREFIX=${1:-"6_unitigs"}
FILENAME_MD=${2:-"statMRAnchors.md"}

tempfile=$(mktemp /tmp/stat_mr_anchor_XXXXXXXX)
trap 'rm -f "$tempfile"' EXIT

printf "%s\t" \
    "Name" "CovCor" "Mapped" \
    "N50Anchor" "Sum" "#" \
    "SumOthers" \
    "median" "MAD" "lower" "upper" |
    sed 's/\t$/\n/' \
    > ${tempfile}

for X in 40 80; do
	for P in $(printf "%03d " {0..2}); do
		if [ ! -e ${DIR_PREFIX}/MRX${X}P${P}/anchor/anchor.fasta ]; then
			continue;
		fi

		pushd ${DIR_PREFIX}/MRX${X}P${P}/ > /dev/null

		SUM_COR=$( cat env.json | jq '.SUM_COR | tonumber' )

        printf "%s\t" \
			"MRX${X}P${P}" \
			$( perl -e "printf qq(%.1f), ${SUM_COR} / 580076;" ) \
            $( cat anchor/env.json | jq '.MAPPED_RATIO | tonumber | (. * 1000 | round) / 1000' ) \
			$( stat_format anchor/anchor.fasta ) \
            $( stat_format anchor/pe.others.fa | cut -f 2 ) \
            $( cat anchor/env.json | jq '.median | tonumber | (. * 10 | round) / 10' ) \
            $( cat anchor/env.json | jq '.MAD    | tonumber | (. * 10 | round) / 10' ) \
            $( cat anchor/env.json | jq '.lower  | tonumber | (. * 10 | round) / 10' ) \
            $( cat anchor/env.json | jq '.upper  | tonumber | (. * 10 | round) / 10' ) |
            sed 's/\t$/\n/'

		popd > /dev/null
	done
done \
>> ${tempfile}

rgr md ${tempfile} --right 2-11 -o ${FILENAME_MD}
echo -e "\nTable: ${FILENAME_MD}\n" >> ${FILENAME_MD}

cat ${FILENAME_MD}
mv ${FILENAME_MD} ${BASH_DIR}/../9_markdown


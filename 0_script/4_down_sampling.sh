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
log_warn 4_down_sampling.sh

parallel --no-run-if-empty --linebuffer -k -j 2 "
    if [ ! -e 2_illumina/Q{1}L{2}/pe.cor.fa.gz ]; then
        exit;
    fi
    echo '==> Q{1}L{2}X{3}'

    if [ -d 4_down_sampling/Q{1}L{2}X{3} ]; then
        echo '    Skip'
        exit
    fi

    # shortcut if cov == all
    if [[ {3} == 'all' ]]; then
        mkdir -p 4_down_sampling/Q{1}L{2}XallP000
        cd 4_down_sampling/Q{1}L{2}XallP000
        gzip -dcf ../../2_illumina/Q{1}L{2}/pe.cor.fa.gz > pe.cor.fa
        cp ../../2_illumina/Q{1}L{2}/env.json .
        exit;
    fi

    # actual sampling
    mkdir -p 4_down_sampling/Q{1}L{2}X{3}
    faops split-about -e \
        2_illumina/Q{1}L{2}/pe.cor.fa.gz \
        \$(( 580076 * {3} )) \
        4_down_sampling/Q{1}L{2}X{3}

    MAX_SERIAL=\$(
        cat 2_illumina/Q{1}L{2}/env.json |
            jq '.SUM_OUT | tonumber | . / 580076 / {3} | floor | . - 1'
    )
    MAX_SERIAL=\$(( \${MAX_SERIAL} < 20 ? \${MAX_SERIAL} : 20 ))

    for i in \$( seq 0 1 \${MAX_SERIAL} ); do
        P=\$( printf '%03d' \${i})
        printf \"  * Part: %s\n\" \${P}

        mkdir -p \"4_down_sampling/Q{1}L{2}X{3}P\${P}\"

        mv  \"4_down_sampling/Q{1}L{2}X{3}/\${P}.fa\" \
            \"4_down_sampling/Q{1}L{2}X{3}P\${P}/pe.cor.fa\"
        cp 2_illumina/Q{1}L{2}/env.json \"4_down_sampling/Q{1}L{2}X{3}P\${P}\"
    done

    " ::: 0 25 30 ::: 0 60 ::: 40 80


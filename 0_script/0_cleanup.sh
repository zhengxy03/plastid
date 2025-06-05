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
    echo $(hnsm n50 -H -N 50 -S -C $@) |
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
log_warn 0_cleanup.sh

# Illumina
parallel --no-run-if-empty --linebuffer -k -j 1 "
    if [ -e 2_illumina/{}.fq.gz ]; then
        rm 2_illumina/{}.fq.gz;
        touch 2_illumina/{}.fq.gz;
    fi
    " ::: clumpify filteredbytile sample trim filter

# insertSize
rm -f 2_illumina/insert_size/*tadpole.contig.fa*

# bwa
if [ -d 3_bwa ]; then
    find 3_bwa -type f -name "genome.fa*"        | parallel --no-run-if-empty -j 1 rm
    find 3_bwa -type f -name "*mate.ba[mi]"      | parallel --no-run-if-empty -j 1 rm
    find 3_bwa -type f -name "*.per-base.bed.gz" | parallel --no-run-if-empty -j 1 rm
fi

# quorum
if [ -d 2_illumina ]; then
    find 2_illumina -type f -name "quorum_mer_db.jf" | parallel --no-run-if-empty -j 1 rm
    find 2_illumina -type f -name "k_u_hash_0"       | parallel --no-run-if-empty -j 1 rm
    find 2_illumina -type f -name "*.tmp"            | parallel --no-run-if-empty -j 1 rm
    find 2_illumina -type f -name "pe.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
    find 2_illumina -type f -name "se.renamed.fastq" | parallel --no-run-if-empty -j 1 rm
    find 2_illumina -type f -name "pe.cor.sub.fa"    | parallel --no-run-if-empty -j 1 rm
    find 2_illumina -type f -name "pe.cor.log"       | parallel --no-run-if-empty -j 1 rm
fi

# down sampling
find . -type f -path "*4_unitigs_bcalm/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_bcalm/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_bcalm/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm

find . -type f -path "*6_unitigs_bcalm/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_bcalm/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_bcalm/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_bifrost/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_bifrost/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_bifrost/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm

find . -type f -path "*6_unitigs_bifrost/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_bifrost/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_bifrost/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_superreads/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_superreads/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_superreads/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm

find . -type f -path "*6_unitigs_superreads/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_superreads/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_superreads/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_tadpole/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_tadpole/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*4_unitigs_tadpole/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm

find . -type f -path "*6_unitigs_tadpole/*" -name "unitigs_K*.fasta"  | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_tadpole/*/anchor*" -name "basecov.txt" | parallel --no-run-if-empty -j 1 rm
find . -type f -path "*6_unitigs_tadpole/*/anchor*" -name "*.sam"       | parallel --no-run-if-empty -j 1 rm

# tempdir
find . -type d -name "\?" | xargs rm -fr

# anchorLong and anchorFill
find . -type d -name "group"         -path "*7_anchor*" | parallel --no-run-if-empty -j 1 rm -fr
find . -type f -name "long.fasta"    -path "*7_anchor*" | parallel --no-run-if-empty -j 1 rm
find . -type f -name ".anchorLong.*" -path "*7_anchor*" | parallel --no-run-if-empty -j 1 rm

# spades
find . -type d -path "*8_spades/*" -not -name "anchor" | parallel --no-run-if-empty -j 1 rm -fr

# quast
find . -type d -name "nucmer_output" | parallel --no-run-if-empty -j 1 rm -fr
find . -type f -path "*contigs_reports/*" -name "*.stdout*" -or -name "*.stderr*" | parallel --no-run-if-empty -j 1 rm

# LSF outputs and dumps
find . -type f -name "output.*" | parallel --no-run-if-empty -j 1 rm
find . -type f -name "core.*"   | parallel --no-run-if-empty -j 1 rm

log_info cat all .md
for NAME in \
    statInsertSize \
    statKAT \
    statFastK \
    statReads \
    statTrimReads \
    statMergeReads \
    statMergeInsert \
    statQuorum \
    ; do
    if [ -e 9_markdown/${NAME}.md ]; then
        echo;
        cat 9_markdown/${NAME}.md;
        echo;
    fi

done

if [ -e statAnchors.md ]; then
    echo;
    cat statAnchors.md;
    echo;
fi
if [ -e 9_markdown/statUnitigsBcalm.md ]; then
    echo;
    cat 9_markdown/statUnitigsBcalm.md
    echo;
fi
if [ -e 9_markdown/statMRUnitigsBcalm.md ]; then
    echo;
    cat 9_markdown/statMRUnitigsBcalm.md;
    echo;
fi
if [ -e 9_markdown/statUnitigsBifrost.md ]; then
    echo;
    cat 9_markdown/statUnitigsBifrost.md
    echo;
fi
if [ -e 9_markdown/statMRUnitigsBifrost.md ]; then
    echo;
    cat 9_markdown/statMRUnitigsBifrost.md;
    echo;
fi
if [ -e 9_markdown/statUnitigsSuperreads.md ]; then
    echo;
    cat 9_markdown/statUnitigsSuperreads.md
    echo;
fi
if [ -e 9_markdown/statMRUnitigsSuperreads.md ]; then
    echo;
    cat 9_markdown/statMRUnitigsSuperreads.md;
    echo;
fi
if [ -e 9_markdown/statUnitigsTadpole.md ]; then
    echo;
    cat 9_markdown/statUnitigsTadpole.md
    echo;
fi
if [ -e 9_markdown/statMRUnitigsTadpole.md ]; then
    echo;
    cat 9_markdown/statMRUnitigsTadpole.md;
    echo;
fi
if [ -e 9_markdown/statMergeAnchors.md ]; then
    echo;
    cat 9_markdown/statMergeAnchors.md;
    echo;
fi
if [ -e 9_markdown/statOtherAnchors.md ]; then
    echo;
    cat 9_markdown/statOtherAnchors.md;
    echo;
fi
if [ -e 9_markdown/statFinal.md ]; then
    echo;
    cat 9_markdown/statFinal.md;
    echo;
fi


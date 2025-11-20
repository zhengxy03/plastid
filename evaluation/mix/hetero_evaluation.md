## preparation
```
#ref_genome
mkdir -p eva/genome
cd eva/genome

for ACCESSION in "NC_000932" "Y08501"; do
    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
    curl $URL -o ${ACCESSION}.fa
done

TAB=$'\t'
cat <<EOF > replace.tsv
NC_000932${TAB}Pt
Y08501${TAB}Mt
EOF

cat NC_000932.fa Y08501.fa |
    faops filter -s stdin stdout |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(echo Pt; echo Mt) genome.fa

cd ..
mkdir -p Sample_Col_G/1_genome
mkdir -p Sample_Ler_XL_4/1_genome
cp genome/genome.fa Sample_Col_G/1_genome/genome.fa
cp genome/genome.fa Sample_Ler_XL_4/1_genome/genome.fa



(head -1 ../Atha_cross/ena/ena_info.tsv && grep -E "(Sample_Col_G|Sample_Ler_XL_4)" ../Atha_cross/ena/ena_info.tsv) > ena_info.tsv

export FOLD=100
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/col0/chr.sizes |
        tsv-summarize --sum 2
)

cat ena_info.tsv |
    tsv-select -H -f name,srr,bases |
    grep -E "Sample_Col_G" |
    perl -nla -F'\t' -e '
        BEGIN { our %seen }
        /^name/ and next;
        $seen{$F[0]} and next;
        my $bases = $F[2];
        $bases =~ s/G$//;
        my $cutoff = $bases * 1000 * 1000 * 1000 / $ENV{GENOME_SIZE} * $ENV{FOLD};
        $cutoff = int $cutoff;
        print join qq(\t), ($F[0], $F[1], $cutoff, $ENV{FOLD});
        $seen{$F[0]}++;
    ' \
    > opts.tsv


export GENOME_SIZE=$(
    cat ~/data/plastid/genome/col0/chr.sizes |
        tsv-summarize --sum 2
)
for FOLD in 1 2 4 8 16 32 64; do
    export FOLD=$FOLD
    echo "计算 FOLD=$FOLD..."
    
    cat ena_info.tsv |
        tsv-select -H -f name,srr,bases |
        grep -E "Sample_Ler_XL_4" |
        perl -nla -F'\t' -e '
            BEGIN { our %seen }
            /^name/ and next;
            $seen{$F[0]} and next;
            my $bases = $F[2];
            $bases =~ s/G$//;
            my $cutoff = $bases * 1000 * 1000 * 1000 / $ENV{GENOME_SIZE} * $ENV{FOLD};
            $cutoff = int $cutoff;
            print join qq(\t), ($F[0], $F[1], $cutoff, $ENV{FOLD});
            $seen{$F[0]}++;
        ' \
        >> opts.tsv
done

#syn
cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ ! -f ena/{2}_1.fastq.gz ]; then
            exit;
        fi
        if [ ! -f ena/{2}_2.fastq.gz ]; then
            exit;
        fi

        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        mkdir -p {1}_{4}/1_genome
        pushd {1}_{4}/1_genome

        cp ../../genome/genome.fa genome.fa
        popd > /dev/null

        mkdir -p {1}_{4}/2_illumina
        pushd {1}_{4}/2_illumina

        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd > /dev/null
    '
```
## run anchr 
```
cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        if [ ! -d {1}_{4} ]; then
            exit;
        fi

        if [ ! -e {1}_{4}/2_illumina/R1.fq.gz ]; then
            exit;
        fi
        if [ ! -e {1}_{4}/2_illumina/R2.fq.gz ]; then
            exit;
        fi

        if bjobs -w | tr -s " " | cut -d " " -f 7 | grep -w "^{1}$"; then
            echo Job {1}_{4} exists
            exit;
        fi

        cd {1}_{4}

        echo {1}_{4}

        rm *.sh
        anchr template \
            --genome 1000000 \
            --parallel 24 \
            --xmx 80g \
            \
            --fastqc \
            --insertsize \
            --fastk \
            \
            --trim "--dedupe --cutoff {3} --cutk 31" \
            --qual "25" \
            --len "60" \
            --filter "adapter artifact" \
            \
            --bwa Q25L60 \
            --gatk

        bsub -q mpi -n 24 -J "{1}_{4}" "
            bash 0_script/2_fastqc.sh
            bash 0_script/2_insert_size.sh
            bash 0_script/2_fastk.sh
            bash 0_script/2_trim.sh
            bash 0_script/9_stat_reads.sh
            bash 0_script/3_bwa.sh
            bash 0_script/3_gatk.sh
        "
    '
```
## mix

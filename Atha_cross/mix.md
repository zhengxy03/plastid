# preperation
```
mkdir -p mix/genome
cp genome/col0/* mix/genome
mkdir -p mix/ena
cp ena/SRR61696* mix/ena
mkdir mix/SRR616966
cd mix/SRR616966
mkdir 1_genome
cd 1_genome
cp ../../genome/genome.fa genome.fa
cp ../../genome/chr.sizes chr.sizes
cd ..
mkdir 2_illumina
cd 2_illumina
cp ../../ena/SRR616966_1.fastq.gz R1.fq.gz
cp ../../ena/SRR616966_2.fastq.gz R2.fq.gz
cd ../../ena
(head -1 ../../ena/ena_info.tsv && grep -E "(SRR616966|SRR616965)" ../../ena/ena_info.tsv) > ena_info.tsv

export FOLD=100
export GENOME_SIZE=$(
    cat genome/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.tsv |
    tsv-select -H -f name,srr,bases |
    grep -E "SRR616966" |
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
    cat genome/chr.sizes |
        tsv-summarize --sum 2
)
for FOLD in 0.5 1 2 4 8 16 32 64; do
    export FOLD=$FOLD
    echo "计算 FOLD=$FOLD..."
    
    cat ena/ena_info.tsv |
        tsv-select -H -f name,srr,bases |
        grep -E "SRR616965" |
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

        mkdir -p {2}_{4}/1_genome
        pushd {2}_{4}/1_genome

        cp ../../genome/genome.fa genome.fa
        cp ../../genome/chr.sizes chr.sizes
        popd > /dev/null

        mkdir -p {2}_{4}/2_illumina
        pushd {2}_{4}/2_illumina

        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd > /dev/null
    '
```
## run anchr 
```
cat opts.tsv | grep -E "0.5" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        if [ ! -d {2}_{4} ]; then
            exit;
        fi

        if [ ! -e {2}_{4}/2_illumina/R1.fq.gz ]; then
            exit;
        fi
        if [ ! -e {2}_{4}/2_illumina/R2.fq.gz ]; then
            exit;
        fi

        if bjobs -w | tr -s " " | cut -d " " -f 7 | grep -w "^{2}$"; then
            echo Job {2}_{4} exists
            exit;
        fi

        cd {2}_{4}

        echo {2}_{4}

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

        bsub -q mpi -n 24 -J "{2}_{4}" "
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

# preperation
```
mkdir -p mix2/genome
cp genome/col0/* mix2/genome
mkdir -p mix2/ena
cp ena/SRR61696* mix2/ena
mkdir mix2/SRR616966
cd mix2/SRR616966
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
cd ..
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
cat opts.tsv |
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
            bash 0_script/0_master.sh
        "
    '
```
# vcf
```
cd ~/zxy/plastid/mix/SRR616966_100/3_gatk
bcftools view --apply-filters PASS --max-alleles 2 --targets Pt --include "AF>0.01" -Oz R.filtered.vcf > Col_Pt.vcf.gz
gunzip Col_Pt.vcf.gz

cd ~/zxy/plastid/mix/SRR616965_1/3_gatk
bcftools view --apply-filters PASS --max-alleles 2 --targets Pt --include "AF>0.01" -Oz R.filtered.vcf > Ler_Pt.vcf.gz
gunzip Ler_Pt.vcf.gz
# 提取Col 100×和Ler 1×的变异位点（仅保留位置）
cd ~/zxy/plastid/mix/mix_results
awk '!/#/ {print $1 "\t" $2}' ~/zxy/plastid/mix/SRR616966_100/3_gatk/Col_Pt.vcf > Col_sites.txt
awk '!/#/ {print $1 "\t" $2}' ~/zxy/plastid/mix/SRR616965_1/3_gatk/Ler_Pt.vcf > Ler_sites.txt

# 筛选Ler特有位点（Ler有而Col没有的位点）
comm -23 <(sort Ler_sites.txt) <(sort Col_sites.txt) > Ler_unique_sites.txt

echo "Col 100×的变异位点数量：$(wc -l < Col_sites.txt)"  # 之前是19个
echo "Ler 1×的变异位点数量：$(wc -l < Ler_sites.txt)"    # 之前是60个
echo "Ler特有的位点数量（Ler有、Col没有）：$(wc -l < Ler_unique_sites.txt)"
```
## mix
```

mkdir -p ~/zxy/plastid/mix/mix_results
cd ~/zxy/plastid/mix/mix_results

# 混合Col 100×和Ler 1×的原始sorted BAM（保留所有区域的reads）
samtools merge -c mix_Col100_Ler1.bam \
  ~/zxy/plastid/mix/SRR616966_100/3_bwa/R.sort.bam \
  ~/zxy/plastid/mix/SRR616965_1/3_bwa/R.sort.bam
samtools index /share/home/wangq/zxy/plastid/mix/mix_results/mix_Col100_Ler1.bam
cd ../genome
gatk CreateSequenceDictionary \
  -R genome.fa \
  -O genome.dict 
cd -
gatk --java-options "-Xmx40g" \
    Mutect2 \
    --native-pair-hmm-threads 16 \
    -R ../genome/genome.fa \
    -I mix_Col100_Ler1.bam \
    -L Pt \
    --mitochondria-mode \
    --max-reads-per-alignment-start 100 \
    --max-mnp-distance 0 \
    -O mix_Pt.raw.vcf

gatk --java-options "-Xmx80g" \
    FilterMutectCalls \
    -R ../genome/genome.fa \
    -V mix_Pt.raw.vcf \
    --max-alt-allele-count 4 \
    --mitochondria-mode \
    -O mix.filtered.vcf
bcftools view --apply-filters PASS --max-alleles 2 --targets Pt --include "AF>0.01" -Oz mix.filtered.vcf > mix_Pt.vcf.gz
gunzip mix_Pt.vcf.gz
```
# detective rate
```
gunzip mix_Pt.vcf.gz
mix_vcf="mix_Pt.vcf"
ler_unique="Ler_unique_sites.txt"
col_sites="Col_sites.txt"


## 混合样本中检测到的Ler特有位点
detected_all=$(grep -Ff $ler_unique $mix_vcf | wc -l)
grep -Ff $ler_unique $mix_vcf | awk '!/#/ {print $1 "\t" $2}' | sort > detected_ler_sites.tmp
detected_true=$(comm -23 detected_ler_sites.tmp <(sort $col_sites) | wc -l)
total_ler=$(wc -l < $ler_unique)

detection_rate_all=$(echo "scale=2; $detected_all / $total_ler * 100" | bc)
detection_rate_true=$(echo "scale=2; $detected_true / $total_ler * 100" | bc)

echo -e "Ler特有位点总数：$total_ler"
echo -e "混合样本检测到的Ler位点（含Col假阳性）：$detected_all → 检测率：$detection_rate_all%"
echo -e "排除Col假阳性后，真实Ler检测位点：$detected_true → 真实检测率：$detection_rate_true%"

comm -23 <(sort Ler_unique_sites.txt) <(sort Col_sites.txt) > true_ler_sites.txt
bgzip mix_Pt.vcf
bcftools query -R true_ler_sites.txt -f '%CHROM:%POS\t%REF\t%ALT\t%AF\n' mix_Pt.vcf.gz | grep -v "^$" > ler_af_final.txt
```

cd ~/zxy/plastid/mix/mix_results
cp ~/zxy/plastid/mix/SRR616965_1/3_bwa/R.sort.bam ./ler1x.sort.bam
cp ~/zxy/plastid/mix/SRR616965_1/3_bwa/R.sort.bai ./ler1x.sort.bai
chmod 755 ./ler1x.sort.bam ./ler1x.sort.bai

samtools cat \
  -o mix_col100_ler1x.bam \
  ~/zxy/plastid/mix/SRR616966_100/3_bwa/R.sort.bam \
  ./ler1x.sort.bam

samtools sort -o mix_col100_ler1x.sorted.bam mix_col100_ler1x.bam
samtools index mix_col100_ler1x.sorted.bam


samtools view -b ./mix_col100_ler1x.sorted.bam Pt > ./mix_Pt.bam

samtools sort -o ./mix_Pt_sorted.bam ./mix_Pt.bam
samtools index ./mix_Pt_sorted.bam
samtools view -b mix_Pt_sorted.bam Pt > mix_Pt_sorted.Pt.bam
samtools index mix_Pt_sorted.Pt.bam
        

```
## heterogenity analysis
```

faops some ~/zxy/plastid/mix/genome/genome.fa <(echo Pt) col_Pt.fa
scp ler_Pt.fa wangq@202.119.37.251:/share/home/wangq/zxy/plastid/mix/mix_results/

#gatk
gatk --java-options "-Xmx80g" Mutect2 \
    -R col_Pt.fa \
    -I mix_Pt_sorted.Pt.bam \
    -L Pt \
    --max-reads-per-alignment-start 100 \
    --max-mnp-distance 0 \
    --annotation StrandBiasBySample \
    -O mix_Pt.raw.vcf

gatk --java-options "-Xmx80g" FilterMutectCalls \
    -R col_Pt.fa \
    -V mix_Pt.raw.vcf \
    --max-alt-allele-count 2 \
    --mitochondria-mode \
    -O mix_Pt.filtered.vcf

cp ~/zxy/plastid/mix/SRR616966_100/3_gatk/R.filtered.vcf ./col100x.filtered.vcf
cp ~/zxy/plastid/mix/SRR616965_1/3_gatk/R.filtered.vcf ./ler1x.filtered.vcf 
bcftools mpileup -f col_Pt.fa ./mix_Pt_sorted.bam | bcftools call -mv -Ov -o ./Pt_variants.vcf







blastn -query ./ler_Pt.fa -subject col_Pt.fa -outfmt 6 -out ./Pt_col_ler_blast.txt
awk '$2 == "Pt" && $9 != $10 {print "Pt\t" $4}' ./Pt_col_ler_blast.txt | sort -u > ./Pt_true_sites.txt
awk '!/#/ && $1 == "Pt" {print "Pt\t" $2}' ./Pt_variants.vcf > ./Pt_detected_sites.txt

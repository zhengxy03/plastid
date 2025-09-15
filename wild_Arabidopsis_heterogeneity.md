# *Arabidopsis thaliana* 1001 Genomes Project
## Basic info
* Genome: GCF_000001735.3, TAIR10, 119.668 Mb<br>
* Chloroplast: NC_000932, Columbia, 154478 bp<br>
* Chloroplast: KX551970, Landsberg erecta, 154515 bp<br> Ler-0，另一种常用野生型，与 Col-0 存在遗传差异
* Mitochondrion: Y08501, 366924 bp<br>
# Project
* PRJNA273563
# Download
* reference
```
mkdir -p Atha_1001/genome
cd Atha_1001/genome
cp ../../Atha_cross/genome/genome.fa .
```
* illumina
    * Download Metadata from NCBI SRA Run Selector via a web browser
        * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA273563
    * Download Accessions info
        * http://1001genomes.org/accessions.html
```
mkdir -p Atha_1001/ena
cd Atha_1001/ena

cat SraRunTable.csv |
    mlr --icsv --otsv cat \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    tsv-filter -H \
        --str-in-fld tissue:leaf \
        --ge AvgSpotLen:90

curl -o accessions.csv \
    https://tools.1001genomes.org/api/accessions.csv?query=SELECT%20*%20FROM%20tg_accessions%20ORDER%20BY%20id

echo -e 'Ecotype\tName\tCS Number\tCountry\tLat\tLong\tCollector\tAdmixture Group\tSequenced by' \
    > accessions.tsv

cat accessions.csv |
    mlr --icsv --implicit-csv-header --otsv cat |
    sed 1d |
    tsv-select -f 1,3,10,4,6,7,8,11,2 \
    >> accessions.tsv

cat accessions.tsv |
    tsv-join -H --key-fields "Ecotype" \
        -f SraRunTable.tsv \
        --append-fields tissue,AvgSpotLen,Bases,Experiment |
    tsv-filter -H \
        --str-in-fld tissue:leaf \
        --ge Bases:1000000000 \
        --le Bases:10000000000 \
        --ge AvgSpotLen:90 \
        --regex Name:'^[\w\d-]+$' |
    tsv-select -H -f Experiment,Name,Country |
    mlr --itsv --ocsv cat |
    sed 1d \
    > source.csv



anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

rgr md ena_info.tsv --fmt | head -n 20

# sort by SRR
cat ena_info.ascp.sh |
    sort > tmp.sh &&
    mv tmp.sh ena_info.ascp.sh

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"

md5sum --check ena_info.md5.txt

cat ena_info.tsv |
    tsv-select -H -f name,srr,bases |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -f {2}_1.fastq.gz ]; then
            echo 1>&2 {1}
            echo {1}
            exit;
        fi
        if [ ! -f {2}_2.fastq.gz ]; then
            echo 1>&2 {1}
            echo {1}
            exit;
        fi

        LENGTH=$(gzip -dcf {2}_1.fastq.gz |
            head -n 100 |
            faops n50 -H stdin
        )
        [[ $LENGTH -le 90 ]] && ( echo 1>&2 {1}; echo {1}; exit; )

        LENGTH=$(gzip -dcf {2}_2.fastq.gz |
            head -n 100 |
            faops n50 -H stdin
        )
        [[ $LENGTH -le 90 ]] && ( echo 1>&2 {1}; echo {1}; exit;  )
    ' |
    sort -r |
    uniq \
    > skip.lst

rsync -avP \
    ~/project/plastid/Atha_1001/ \
    wangq@202.119.37.251:zxy/plastid/Atha_1001/ena
```
| name        | srx       | platform | layout | ilength | srr        |      spots | bases |
|-------------|-----------|----------|--------|---------|------------|-----------:|-------|
| 11C1        | SRX972788 | ILLUMINA | PAIRED |         | SRR1946105 | 20,357,230 | 3.79G |
| 328PNA062   | SRX972701 | ILLUMINA | PAIRED |         | SRR1946018 | 73,235,208 | 6.62G |
| 627ME-13Y1  | SRX972634 | ILLUMINA | PAIRED |         | SRR1945951 | 44,902,951 | 8.28G |
| 627ME-1MI1  | SRX972635 | ILLUMINA | PAIRED |         | SRR1945952 | 48,675,298 | 8.98G |
| 627RMX-1MN4 | SRX972632 | ILLUMINA | PAIRED |         | SRR1945949 | 37,102,228 | 6.82G |
| 627RMX-1MN5 | SRX972633 | ILLUMINA | PAIRED |         | SRR1945950 | 44,085,175 | 8.11G |
| ANH-1       | SRX972422 | ILLUMINA | PAIRED |         | SRR1945739 | 10,783,643 | 2.03G |
| ARGE-1-15   | SRX973146 | ILLUMINA | PAIRED |         | SRR1946463 | 25,951,374 | 4.83G |
| ARR-17      | SRX973157 | ILLUMINA | PAIRED |         | SRR1946474 | 14,893,935 | 2.77G |
| Aa-0        | SRX972490 | ILLUMINA | PAIRED |         | SRR1945807 |  9,051,264 | 1.7G  |
| Abd-0       | SRX972484 | ILLUMINA | PAIRED |         | SRR1945801 | 13,961,306 | 2.63G |
| Adam-1      | SRX972881 | ILLUMINA | PAIRED |         | SRR1946198 | 12,010,826 | 2.24G |
| Ag-0        | SRX972432 | ILLUMINA | PAIRED |         | SRR1945749 | 10,141,217 | 1.91G |
| Aiell-1     | SRX972915 | ILLUMINA | PAIRED |         | SRR1946232 | 12,869,053 | 2.4G  |
| Aitba-1     | SRX972878 | ILLUMINA | PAIRED |         | SRR1946195 | 11,730,031 | 2.18G |
| Ak-1        | SRX972485 | ILLUMINA | PAIRED |         | SRR1945802 |  8,524,299 | 1.6G  |
| Alst-1      | SRX972486 | ILLUMINA | PAIRED |         | SRR1945803 |  9,347,634 | 1.76G |
| Alt-1       | SRX973020 | ILLUMINA | PAIRED |         | SRR1946337 | 24,312,469 | 4.53G |
# Symlink
```
export FOLD=8
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/col0/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.tsv |
    tsv-select -H -f name,srr,bases |
    tsv-join -H --filter-file ena/skip.lst --key-fields name --exclude |
    perl -nla -F'\t' -e '
        BEGIN { our %seen }
        /^name/ and next;
        $seen{$F[0]} and next;
        my $bases = $F[2];
        $bases =~ s/G$//;
        my $cutoff = $bases * 1000 * 1000 * 1000 / $ENV{GENOME_SIZE} * $ENV{FOLD};
        $cutoff = int $cutoff;
        print join qq(\t), ($F[0], $F[1], $cutoff);
        $seen{$F[0]}++;
    ' \
    > opts.tsv

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

        mkdir -p {1}/1_genome
        pushd {1}/1_genome

        cp ../../genome/genome.fa genome.fa
        popd > /dev/null

        mkdir -p {1}/2_illumina
        pushd {1}/2_illumina

        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd > /dev/null
    '
```
# Run
```
cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        if [ ! -d {1} ]; then
            exit;
        fi

        if [ ! -e {1}/2_illumina/R1.fq.gz ]; then
            exit;
        fi
        if [ ! -e {1}/2_illumina/R2.fq.gz ]; then
            exit;
        fi

        if bjobs -w | tr -s " " | cut -d " " -f 7 | grep -w "^{1}$"; then
            echo Job {1} exists
            exit;
        fi

        cd {1}

        echo {1}

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

        bsub -q mpi -n 24 -J "{1}" "
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
# heterogeneity
```
mkdir -p vcf
cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -d "{1}" ] || [ ! -f "{1}/3_gatk/R.filtered.vcf" ]; then
            echo "样本 {1} 文件缺失，跳过"
            continue
        fi

        echo "处理样本: {1}"
        
        # 处理VCF：重命名样本、筛选Pt、过滤PASS、双等位基因、最小AF
        bcftools reheader --samples <(echo {1}) "{1}/3_gatk/R.filtered.vcf" |
            bcftools view \
                --apply-filters PASS \
                --max-alleles 2 \
                --targets Pt \
                -Oz |
            bcftools view --include "AF>0.01" -Oz -o vcf/{1}.vcf.gz

        # 索引VCF
        bcftools index -f vcf/{1}.vcf.gz
    '
bcftools merge --merge all -l <(
        cat opts.tsv |
            cut -f 1 |
            parallel -k -j 1 ' [ -f vcf/{}.vcf.gz ] && echo "vcf/{}.vcf.gz" '
    ) \
    > Atha_1001.vcf

bcftools stats Atha_1001.vcf > Atha_1001.stats

cat opts.tsv | cut -f 1 > sample_list.txt
GENOME_SIZE=154478

echo -e "样本名称\t异质位点数\t比例(%)" > sample_heterozygosity_summary.tsv

while read sample; do
    echo "[INFO] 处理样本: $sample"

    bcftools view -s $sample Atha_1001.vcf -Ov -o tmp_${sample}.vcf

    variant_count=$(bcftools query -f '[%GT\n]' tmp_${sample}.vcf | grep -v -E '0/0|\.\/\.' | wc -l)

    ratio=$(echo "scale=6; $variant_count / $GENOME_SIZE * 100" | bc)
    ratio=$(printf "%0.6f" $ratio)

    echo -e "$sample\t$variant_count\t$ratio" >> sample_heterozygosity_summary.tsv

    rm tmp_${sample}.vcf
done < sample_list.txt

echo -e "| 样本名称 | 异质位点数 | 比例(%) |" > sample_heterozygosity_summary.md
echo -e "| --- | --- | --- |" >> sample_heterozygosity_summary.md

tail -n +2 sample_heterozygosity_summary.tsv | awk -F'\t' '{printf("| %s | %s | %s |\n", $1, $2, $3)}' >> sample_heterozygosity_summary.md
tail -n +2 sample_heterozygosity_summary.tsv \
| perl -F'\t' -ane 'printf("| %s | %s | %s |\n", $F[0], $F[1], $F[2])' \
>> sample_heterozygosity_summary_perl.md
```
| 样本名称 | 异质位点数 | 比例(%) |
| --- | --- | --- |
| 11C1 | 182 | 0.117800 |
| 627ME-13Y1 | 81 | 0.052400 |
| 627ME-1MI1 | 70 | 0.045300 |
| 627RMX-1MN4 | 127 | 0.082200 |
| 627RMX-1MN5 | 125 | 0.080900 |
| ANH-1 | 133 | 0.086000 |
| ARGE-1-15 | 192 | 0.124200 |
| ARR-17 | 153 | 0.099000 |
| Aa-0 | 73 | 0.047200 |
| Abd-0 | 143 | 0.092500 |
| Adam-1 | 187 | 0.121000 |
| Ag-0 | 75 | 0.048500 |
| Aiell-1 | 107 | 0.069200 |
| Aitba-1 | 233 | 0.150800 |
| Ak-1 | 139 | 0.089900 |
| Alst-1 | 117 | 0.075700 |
| Alt-1 | 128 | 0.082800 |
```
library(ggplot2)
library(readr)
library(dplyr)


df <- read_tsv("sample_heterozygosity_summary.tsv")
df <- df %>% arrange(`比例(%)`)
df$样本名称 <- factor(df$样本名称, levels = df$样本名称)

p <- ggplot(df, aes(x = `比例(%)`, y = 样本名称)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "heterozygosity ratio(%)", y = "Sample") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_blank()
  )

ggsave("heterozygosity_barplot.png", plot = p, width = 8, height = 6, dpi = 300)
```
![heterozygosity](./pic/heterozygosity_barplot.png)
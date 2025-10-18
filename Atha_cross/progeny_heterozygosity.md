# Progeny Chloroplast Heterozygosity
Use Atha_cross(*Arabidopsis thaliana* Columbia X Ler) as example
## Basic info
* Genome: GCF_000001735.3, TAIR10, 119.668 Mb<br>
* Chloroplast: NC_000932, Columbia, 154478 bp<br>
* Chloroplast: KX551970, Landsberg erecta, 154515 bp<br> Ler-0，另一种常用野生型，与 Col-0 存在遗传差异
* Mitochondrion: Y08501, 366924 bp<br>
## Project
* https://www.pnas.org/content/109/51/20992.long<br>
* PRJNA178613<br>
## Download
* Reference
```
mkdir -p ~/data/plastid/Atha_cross/genome
cd ~/data/plastid/Atha_cross/genome

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
```
* illumina
    * Download Metadata from NCBI SRA Run Selector via a web browser
        * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA178613
```
cat SraRunTable.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f Experiment,"Sample\ Name",Bases \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    sed '1 s/^/#/' |
    keep-header -- tsv-sort -k2,2 -k3,3nr |
    tsv-uniq -H -f "Sample\ Name" --max 1 |
    mlr --itsv --ocsv cat \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

rgr md ena_info.tsv --fmt

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"
```
## Symlink
```
cd ~/data/plastid/Atha_cross/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/col0/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.tsv |
    tsv-select -H -f name,srr,bases |
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
## Run anchr
```
cd Atha_cross
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
## calculate heterozygosity frequency
* frequency
```
PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross_fold0"
CHLORO_NAME="Pt"
opts_file="$PROJECT_ROOT/opts.tsv"

# 提取样本名（去除前后空格和空行）
samples=$(awk -F'\t' 'NR>0 && $1!="" {gsub(/^[ \t\r]+|[ \t\r]+$/, ""); print $1}' "$opts_file")

# 循环处理样本
for sample in $samples; do
    echo -e "\n===== 处理样本: $sample ====="
    
    bam_file="$PROJECT_ROOT/$sample/3_bwa/R.sort.bam"
    bam_index="$PROJECT_ROOT/$sample/3_bwa/R.sort.bai"
    gatk_vcf="$PROJECT_ROOT/$sample/3_gatk/R.filtered.vcf"
    gatk_vcf_gz="$PROJECT_ROOT/$sample/3_gatk/R.filtered.vcf.gz"
    result_file="$PROJECT_ROOT/$sample/chloroplast_heteroplasmy_freq.txt"
    chloro_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_from_gatk.vcf.gz"
    chloro_vcf_tbi="$chloro_vcf.tbi"
    temp_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_temp.vcf"
    filtered_temp_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_filtered_temp.vcf"
    log_file="$PROJECT_ROOT/$sample/chloro_extraction.log"

    if [ ! -f "$bam_file" ] || [ ! -f "$bam_index" ] || [ ! -f "$gatk_vcf" ]; then
        echo "警告: 缺少必要文件（BAM/索引/GATK VCF），跳过样本 $sample"
        continue
    fi

    echo "开始处理样本 $sample" | tee -a "$log_file"

    # 提取叶绿体区域reads
    samtools view -b -q 10 -o "$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam" "$bam_file" "$CHLORO_NAME"

    samtools sort -o "$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam" "$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam"
    samtools index "$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam"

    avg_depth=$(samtools depth "$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam" | perl calc_coverage.pl)
    echo "叶绿体平均覆盖度: ${avg_depth}×" | tee -a "$log_file"

    # 压缩并索引 GATK VCF
    bgzip -c "$gatk_vcf" > "$gatk_vcf_gz"
    tabix -p vcf "$gatk_vcf_gz"

    # 提取叶绿体变异
    echo '提取叶绿体变异' | tee -a "$log_file"
    bcftools view --apply-filters PASS --targets "$CHLORO_NAME" -Oz "$gatk_vcf_gz" > "$filtered_temp_vcf.gz"
    tabix -p vcf "$filtered_temp_vcf.gz"

    # 解压供下游 Perl 处理
    bcftools view "$filtered_temp_vcf.gz" > "$temp_vcf"

    var_count=$(grep -v "^#" "$temp_vcf" | wc -l)
    echo "提取到的叶绿体变异总数: $var_count" | tee -a "$log_file"

    if [ "$var_count" -eq 0 ]; then
        echo '未检测到任何变异' | tee -a "$log_file"
        touch "$result_file"
        continue
    fi

    # 压缩最终叶绿体VCF
    bgzip -c "$temp_vcf" > "$chloro_vcf"
    tabix -p vcf "$chloro_vcf"
    echo '生成最终叶绿体VCF完成' | tee -a "$log_file"

    # 计算异质性频率，并过滤 <1%
    echo '计算异质性频率并过滤小于1%位点...' | tee -a "$log_file"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t%QUAL\n' "$chloro_vcf" | \
        perl calc_frequency.pl > "$result_file"

    result_count=$(wc -l "$result_file" | cut -d' ' -f1)
    echo "保留异质性变异数 (≥1%): $result_count" | tee -a "$log_file"

    # 清理临时文件
    rm -f "$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam"
    rm -f "$gatk_vcf_gz" "$gatk_vcf_gz.tbi"
    rm -f "$temp_vcf" "$filtered_temp_vcf.gz" "$filtered_temp_vcf.gz.tbi"

    echo "样本 $sample 处理完成" | tee -a "$log_file"
done
```
* 汇总
```
PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross"
OPTS_FILE="$PROJECT_ROOT/opts.tsv"
OUTPUT_PREFIX="$PROJECT_ROOT/chloro_multi_stats"

CSV_OUTPUT="${OUTPUT_PREFIX}.csv"
MD_OUTPUT="${OUTPUT_PREFIX}.md"


EXCLUDED_SAMPLES=("Sample_c1c2" "Sample_l2c2" "Sample_l2l3" "Sample_l4c1" "Sample_l4l3" "")
samples=($(perl -ne '
    chomp;
    next if /^\s*$/;                # 跳过空行
    my @f = split(/\t/);
    $f[0] =~ s/^\s+|\s+$//g;        # 去掉前后空格
    next if $f[0] =~ /^(Sample_c1c2|Sample_l2c2|Sample_l2l3|Sample_l4c1|Sample_l4l3)$/;  # 排除指定样本
    print "$f[0]\n";
' "$OPTS_FILE"))

echo "找到 ${#samples[@]} 个有效样本："
printf "→ %s\n" "${samples[@]}"


echo "样本名,有效变异数,高频(≥10%),中高频(5%-10%),低频(1%-5%),平均异质性频率(%),最高异质性频率(%),平均总reads数(×),最高频率位点(Pt:POS)" > "$CSV_OUTPUT"
for sample in "${samples[@]}"; do
    result_file="$PROJECT_ROOT/$sample/chloroplast_heteroplasmy_freq.txt"
    
    total_vars=$(wc -l < "$result_file")
    freq_1_5=$(awk '$5>=0.01 && $5<0.05 {c++} END{print c+0}' "$result_file")
    freq_5_10=$(awk '$5>=0.05 && $5<0.1 {c++} END{print c+0}' "$result_file")
    freq_10_plus=$(awk '$5>=0.1 {c++} END{print c+0}' "$result_file")
    avg_freq=$(awk '{sum+=$5} END{printf "%.2f", sum/NR*100}' "$result_file")
    max_freq=$(awk '{max=$5>max?$5:max} END{printf "%.2f", max*100}' "$result_file")
    avg_reads=$(awk '{sum+=$6} END{printf "%.1f", sum/NR}' "$result_file")
    max_pos=$(awk -v max=0 '$5>max{max=$5; pos=$2} END{print pos}' "$result_file")

    echo "${sample},${total_vars},${freq_10_plus},${freq_5_10},${freq_1_5},${avg_freq},${max_freq},${avg_reads},${max_pos}" >> "$CSV_OUTPUT"
done

cat > "$MD_OUTPUT" << EOF
# 叶绿体异质性变异多样本统计报告
**生成时间**：$(date "+%Y-%m-%d %H:%M:%S")  
**项目目录**：$PROJECT_ROOT  
**opts.tsv路径**：$OPTS_FILE  
**排除样本**：${EXCLUDED_SAMPLES[*]}  
**有效样本数量**：${#samples[@]}

| 样本名 | 有效变异数 | 高频(≥10%) | 中高频(5%-10%) | 低频(1%-5%) | 平均频率(%) | 最高频率(%) | 平均覆盖(×) | 最高频率位点(Pt:POS) |
|---------|------------|------------|----------------|-------------|-------------|-------------|--------------|-----------------------|
EOF

for sample in "${samples[@]}"; do
    result_file="$PROJECT_ROOT/$sample/chloroplast_heteroplasmy_freq.txt"
    [ ! -f "$result_file" ] && continue

    total_vars=$(wc -l < "$result_file")
    freq_1_5=$(awk '$5>=0.01 && $5<0.05 {c++} END{print c+0}' "$result_file")
    freq_5_10=$(awk '$5>=0.05 && $5<0.1 {c++} END{print c+0}' "$result_file")
    freq_10_plus=$(awk '$5>=0.1 {c++} END{print c+0}' "$result_file")
    avg_freq=$(awk '{sum+=$5} END{printf "%.2f", sum/NR*100}' "$result_file")
    max_freq=$(awk '{max=$5>max?$5:max} END{printf "%.2f", max*100}' "$result_file")
    avg_reads=$(awk '{sum+=$6} END{printf "%.1f", sum/NR}' "$result_file")
    max_pos=$(awk -v max=0 '$5>max{max=$5; pos=$2} END{print pos}' "$result_file")

    echo "| ${sample} | ${total_vars} | ${freq_10_plus} | ${freq_5_10} | ${freq_1_5} | ${avg_freq} | ${max_freq} | ${avg_reads} | Pt:${max_pos} |" >> "$MD_OUTPUT"
done
```

| 样本名       | 有效变异数 | 高频变异数（≥10%） | 中高频变异数（5%-10%） | 低频变异数（1%-5%） | 平均异质性频率（%） | 最高异质性频率（%） | 平均总reads数（×） | 最高频率位点（Pt:POS） |
|--------------|------------|--------------------|------------------------|--------------------|---------------------|---------------------|--------------------|------------------------|
| Sample_14 | 104 | 6 | 21 | 77 | 4.62 | 84.69 | 1513.5 | Pt:28672 |
| Sample_18 | 111 | 8 | 24 | 79 | 4.69 | 84.67 | 1437.8 | Pt:28672 |
| Sample_19 | 115 | 10 | 21 | 84 | 7.98 | 100.00 | 859.0 | Pt:136734 |
| Sample_20 | 110 | 6 | 26 | 78 | 4.53 | 84.19 | 1222.8 | Pt:28672 |
| Sample_21 | 106 | 4 | 27 | 75 | 4.63 | 84.63 | 1051.8 | Pt:28672 |
| Sample_4 | 100 | 6 | 21 | 73 | 5.51 | 100.00 | 954.7 | Pt:36321 |
| Sample_5 | 102 | 7 | 26 | 69 | 8.35 | 100.00 | 822.1 | Pt:136734 |
| Sample_6 | 193 | 8 | 39 | 146 | 4.00 | 81.46 | 2126.7 | Pt:28672 |
| Sample_7 | 112 | 10 | 17 | 85 | 8.30 | 100.00 | 1867.2 | Pt:136734 |
| Sample_8 | 96 | 9 | 14 | 73 | 8.62 | 100.00 | 670.9 | Pt:136734 |
| Sample_Col_G | 106 | 9 | 19 | 78 | 4.69 | 84.05 | 2619.9 | Pt:28672 |
| Sample_Ler_XL_4 | 266 | 91 | 36 | 139 | 30.49 | 100.00 | 5239.3 | Pt:4870 |
| Sample_c41 | 104 | 6 | 21 | 77 | 4.49 | 83.85 | 1717.3 | Pt:28672 |
| Sample_c42 | 102 | 9 | 20 | 73 | 4.77 | 81.92 | 1409.8 | Pt:28672 |
| Sample_c45 | 100 | 8 | 19 | 73 | 4.95 | 83.30 | 2163.1 | Pt:28672 |
| Sample_c47 | 92 | 6 | 16 | 70 | 4.69 | 85.96 | 773.7 | Pt:28672 |
| Sample_c48 | 99 | 8 | 18 | 73 | 4.77 | 85.14 | 1594.1 | Pt:28672 |
| Sample_c51 | 104 | 9 | 16 | 79 | 4.89 | 84.39 | 1625.3 | Pt:28672 |
| Sample_c52 | 102 | 7 | 13 | 82 | 4.49 | 87.40 | 1564.1 | Pt:28672 |
| Sample_c54 | 100 | 7 | 15 | 78 | 4.52 | 82.81 | 1276.6 | Pt:28672 |
| Sample_c57 | 97 | 6 | 16 | 75 | 4.50 | 83.96 | 1745.2 | Pt:28672 |
| Sample_c61 | 95 | 8 | 16 | 71 | 4.78 | 84.27 | 1606.9 | Pt:28672 |
| Sample_c62 | 106 | 6 | 19 | 81 | 4.51 | 84.35 | 1560.3 | Pt:28672 |
| Sample_c63 | 99 | 8 | 15 | 76 | 4.58 | 84.34 | 1445.2 | Pt:28672 |
| Sample_c64 | 96 | 7 | 15 | 74 | 4.62 | 84.32 | 1207.2 | Pt:28672 |
| Sample_c65 | 198 | 10 | 30 | 158 | 4.05 | 78.42 | 2000.1 | Pt:28672 |
| Sample_c66 | 100 | 9 | 17 | 74 | 5.17 | 84.58 | 2230.4 | Pt:28672 |
| Sample_c73 | 94 | 4 | 20 | 70 | 4.54 | 85.79 | 1374.9 | Pt:28672 |
| Sample_c81 | 107 | 8 | 21 | 78 | 4.65 | 83.00 | 1413.9 | Pt:28672 |
| Sample_c82 | 91 | 6 | 18 | 67 | 5.34 | 83.82 | 2006.1 | Pt:28672 |
| Sample_c83 | 112 | 9 | 17 | 86 | 4.63 | 82.83 | 2545.9 | Pt:28672 |
| Sample_c84 | 114 | 7 | 21 | 86 | 4.59 | 85.05 | 2067.3 | Pt:28672 |
| Sample_c85 | 109 | 7 | 18 | 84 | 4.77 | 84.99 | 2316.5 | Pt:28672 |
| Sample_c87 | 107 | 8 | 18 | 81 | 4.85 | 82.09 | 2039.6 | Pt:28672 |
| Sample_c88 | 112 | 10 | 17 | 85 | 4.96 | 83.66 | 2336.7 | Pt:28672 |
| Sample_c89 | 107 | 8 | 19 | 80 | 4.76 | 83.43 | 1743.7 | Pt:28672 |
| Sample_c90 | 110 | 9 | 16 | 85 | 4.78 | 85.34 | 2493.4 | Pt:28672 |
| Sample_c91 | 109 | 10 | 20 | 79 | 5.12 | 81.59 | 2433.9 | Pt:28672 |
| Sample_c92 | 110 | 8 | 17 | 85 | 4.65 | 82.75 | 1914.9 | Pt:28672 |
| Sample_c93 | 111 | 8 | 19 | 84 | 4.58 | 82.15 | 1719.8 | Pt:28672 |
| Sample_c94 | 103 | 9 | 13 | 81 | 4.55 | 85.02 | 3791.4 | Pt:28672 |
| Sample_c95 | 106 | 9 | 18 | 79 | 4.63 | 83.70 | 3651.7 | Pt:28672 |

## analysis heterogeneity origin
```
WORKDIR="/share/home/wangq/zxy/plastid/Atha_cross"
MOTHER_FREQ="${WORKDIR}/Sample_Col_G/chloroplast_heteroplasmy_freq.txt"
FATHER_FREQ="${WORKDIR}/Sample_Ler_XL_4/chloroplast_heteroplasmy_freq.txt"

declare -A mom_map
declare -A dad_map

while read -r chr pos ref alt freq depth; do
    key="${chr}_${pos}_${alt}"
    mom_map["$key"]="$freq"
done < <(tail -n +1 "$MOTHER_FREQ")

while read -r chr pos ref alt freq depth; do
    key="${chr}_${pos}_${alt}"
    dad_map["$key"]="$freq"
done < <(tail -n +1 "$FATHER_FREQ")

# 获取子代样本
SAMPLES=($(find "$WORKDIR" -maxdepth 1 -type d -name "Sample_*" \
    | grep -vE "Sample_Col_G|Sample_Ler_XL_4|Sample_l2c2|Sample_l2l3|Sample_l4c1|Sample_l4l3" \
    | sort))

echo "有效子代样本数: ${#SAMPLES[@]}"

for SAMPLE_DIR in "${SAMPLES[@]}"; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "----- 处理 ${SAMPLE} -----"

    FREQ_FILE="${SAMPLE_DIR}/chloroplast_heteroplasmy_freq.txt"
    OUT_FILE="${SAMPLE_DIR}/chloroplast_hetero_parent_ratio.txt"

    # 新增“变异reads数”列标题
    echo -e "染色体\t位点\t参考碱基\t变异碱基\t异质性频率\t总reads数\t变异reads数\t母本频率\t父本频率\t来源判断" > "$OUT_FILE"

    tail -n +1 "$FREQ_FILE" | while read -r chr pos ref alt hetfreq depth; do
        key="${chr}_${pos}_${alt}"
        mom_f=${mom_map["$key"]:-0}
        dad_f=${dad_map["$key"]:-0}

        # 计算变异reads数：总reads数 × 异质性频率，并用printf四舍五入取整数
        # 先通过bc计算乘积，再用printf "%.0f"四舍五入
        var_reads=$(echo "$depth * $hetfreq" | bc -l | xargs printf "%.0f")

        # 判断来源
        if (( $(echo "$mom_f>0" | bc -l) )) && (( $(echo "$dad_f==0" | bc -l) )); then
            origin="母本"
        elif (( $(echo "$dad_f>0" | bc -l) )) && (( $(echo "$mom_f==0" | bc -l) )); then
            origin="父本"
        elif (( $(echo "$mom_f>0" | bc -l) )) && (( $(echo "$dad_f>0" | bc -l) )); then
            origin="混合（父母均含）"
        else
            origin="自发变异"
        fi

        # 输出时增加“变异reads数”列
        echo -e "${chr}\t${pos}\t${ref}\t${alt}\t${hetfreq}\t${depth}\t${var_reads}\t${mom_f}\t${dad_f}\t${origin}" >> "$OUT_FILE"
    done

    echo "已完成 ${SAMPLE} -> $(basename "$OUT_FILE")"
done
```
可视化
```R
library(dplyr)
library(ggplot2)
library(readr)
library(scales)

workdir <- "/share/home/wangq/zxy/plastid/Atha_cross"

exclude_samples <- c(
  "Sample_Col_G", "Sample_Ler_XL_4",
  "Sample_l2c2", "Sample_l2l3",
  "Sample_l4c1", "Sample_l4l3", "Sample_c1c2"
)

all_samples <- list.dirs(workdir, full.names = FALSE, recursive = FALSE)
samples <- setdiff(all_samples, exclude_samples)

cat("Valid offspring samples:", length(samples), "\n")

all_stats <- data.frame()

for (sample in samples) {
  infile <- file.path(workdir, sample, "chloroplast_hetero_parent_ratio.txt")

  if (!file.exists(infile)) {
    cat("⚠️ File missing, skipping:", infile, "\n")
    next
  }

  df <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  df$Sample <- sample
  all_stats <- bind_rows(all_stats, df)
}

all_stats$来源判断 <- gsub("（.*", "", all_stats$来源判断)  # 去掉括号说明
all_stats$来源判断 <- trimws(all_stats$来源判断)            # 去掉空格

all_stats$Origin <- recode(all_stats$来源判断,
                           "母本" = "Mother",
                           "父本" = "Father",
                           "混合" = "Mixed",
                           "自发变异" = "Spontaneous")

summary_stats <- all_stats %>%
  group_by(Sample, Origin) %>%
  summarise(total_var_reads = sum(`变异reads数`, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(sum_var_reads = sum(total_var_reads, na.rm = TRUE)) %>%
  mutate(proportion = ifelse(sum_var_reads > 0, total_var_reads / sum_var_reads, 0)) %>%
  ungroup()

summary_stats$proportion <- pmin(summary_stats$proportion, 1)
summary_stats$proportion[is.na(summary_stats$proportion)] <- 0

report <- summary_stats %>%
  group_by(Sample) %>%
  summarise(
    Mother = round(sum(proportion[Origin == "Mother"], na.rm = TRUE), 3),
    Father = round(sum(proportion[Origin == "Father"], na.rm = TRUE), 3),
    Mixed = round(sum(proportion[Origin == "Mixed"], na.rm = TRUE), 3),
    Spontaneous = round(sum(proportion[Origin == "Spontaneous"], na.rm = TRUE), 3)
  ) %>%
  mutate(Total = Mother + Father + Mixed + Spontaneous)

out_table <- file.path(workdir, "chloroplast_hetero_summary.csv")
write.csv(report, out_table, row.names = FALSE)

summary_stats$Origin <- factor(
  summary_stats$Origin,
  levels = c("Mother", "Father", "Mixed", "Spontaneous")
)

fill_colors <- c("Mother" = "black", "Father" = "grey50", "Mixed" = "grey80", "Spontaneous" = "grey30")

p <- ggplot(summary_stats, aes(x = Sample, y = proportion, fill = Origin)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    x = "Sample",
    y = "Variant Read Proportion",
    fill = "Origin Type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

out_plot <- file.path(workdir, "chloroplast_hetero_summary_plot.png")
ggsave(out_plot, p, width = 12, height = 6, dpi = 300)
```
![origin](./results/chloroplast_hetero_summary_plot.png)

```

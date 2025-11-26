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
mkdir -p ~/zxy/plastid/Atha_cross_new/genome
cd ~/zxy/plastid/Atha_cross_new/genome

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
cd ~/data/plastid/Atha_cross_new/

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
cd Atha_cross_new
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
PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross_new"
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

    # 初始化日志文件
    echo "开始处理样本 $sample" > "$log_file"

    if [ ! -f "$bam_file" ] || [ ! -f "$bam_index" ] || [ ! -f "$gatk_vcf" ]; then
        echo "警告: 缺少必要文件（BAM/索引/GATK VCF），跳过样本 $sample" | tee -a "$log_file"
        # 创建空的结果文件
        echo -e "Chr\tPos\tRef\tAlt\tAF\tDP\tQUAL" > "$result_file"
        continue
    fi

    echo "开始处理样本 $sample" | tee -a "$log_file"

    # 提取叶绿体区域reads
    samtools view -b -q 10 -o "$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam" "$bam_file" "$CHLORO_NAME" 2>/dev/null

    samtools sort -o "$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam" "$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam" 2>/dev/null
    samtools index "$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam" 2>/dev/null

    # 计算平均深度（如果calc_coverage.pl不存在，使用替代方法）
    if [ -f "calc_coverage.pl" ]; then
        avg_depth=$(samtools depth "$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam" 2>/dev/null | perl calc_coverage.pl 2>/dev/null || echo "N/A")
    else
        # 使用awk计算平均深度
        avg_depth=$(samtools depth "$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam" 2>/dev/null | awk '{sum+=$3} END {if(NR>0) printf "%.2f", sum/NR; else print "N/A"}' 2>/dev/null || echo "N/A")
    fi
    echo "叶绿体平均覆盖度: ${avg_depth}×" | tee -a "$log_file"

    # 压缩并索引 GATK VCF
    bgzip -c "$gatk_vcf" > "$gatk_vcf_gz" 2>/dev/null
    tabix -p vcf "$gatk_vcf_gz" 2>/dev/null

    # 提取叶绿体变异 - 使用更严格的过滤标准
    echo '提取叶绿体变异（使用严格过滤）' | tee -a "$log_file"
    bcftools view \
        --apply-filters PASS \
        --types snps \
        --max-alleles 2 \
        --targets "$CHLORO_NAME" \
        --include "AF>0.01" \
        -Oz "$gatk_vcf_gz" > "$filtered_temp_vcf.gz" 2>/dev/null
    
    if [ -f "$filtered_temp_vcf.gz" ]; then
        tabix -p vcf "$filtered_temp_vcf.gz" 2>/dev/null
    fi

    # 检查过滤后的变异数量
    if [ -f "$filtered_temp_vcf.gz" ]; then
        var_count=$(bcftools view "$filtered_temp_vcf.gz" 2>/dev/null | grep -v "^#" | wc -l 2>/dev/null || echo 0)
    else
        var_count=0
    fi
    echo "过滤后的叶绿体变异总数 (PASS, 双等位基因, AF>1%): $var_count" | tee -a "$log_file"

    # 创建结果文件头（无论是否有变异都创建）
    echo -e "Chr\tPos\tRef\tAlt\tAF\tDP\tQUAL" > "$result_file"

    if [ "$var_count" -eq 0 ]; then
        echo '未检测到任何符合条件的变异' | tee -a "$log_file"
        # 即使没有变异，也保留空的结果文件（只有标题行）
        echo "创建空的结果文件: $result_file" | tee -a "$log_file"
    else
        # 首先检查VCF文件中AF字段的位置
        echo "检查VCF文件结构..." | tee -a "$log_file"
        bcftools view -h "$filtered_temp_vcf.gz" 2>/dev/null | grep -E "(INFO|FORMAT).*AF" | tee -a "$log_file"

        # 直接提取AF频率信息 - 修正格式
        echo '提取AF频率信息...' | tee -a "$log_file"
        
        # 尝试不同的AF字段提取方式
        # 方式1：如果AF在FORMAT区域
        if bcftools view -h "$filtered_temp_vcf.gz" 2>/dev/null | grep -q "FORMAT.*AF"; then
            echo "AF字段在FORMAT区域，使用[ %AF]格式提取" | tee -a "$log_file"
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\t%DP\t%QUAL\n' "$filtered_temp_vcf.gz" 2>/dev/null >> "$result_file"
        
        # 方式2：如果AF在INFO区域  
        elif bcftools view -h "$filtered_temp_vcf.gz" 2>/dev/null | grep -q "INFO.*AF"; then
            echo "AF字段在INFO区域，使用%AF格式提取" | tee -a "$log_file"
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%DP\t%QUAL\n' "$filtered_temp_vcf.gz" 2>/dev/null >> "$result_file"
        
        # 方式3：如果没有AF字段，从AD字段计算
        else
            echo "未找到AF字段，从AD字段计算频率" | tee -a "$log_file"
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t%DP\t%QUAL\n' "$filtered_temp_vcf.gz" 2>/dev/null | \
            awk 'BEGIN{OFS="\t"} {
                split($5, ad, ",");
                total = ad[1] + ad[2];
                af = (total > 0) ? ad[2]/total : 0;
                print $1, $2, $3, $4, af, $6, $7
            }' >> "$result_file" 2>/dev/null
        fi

        result_count=$(wc -l < "$result_file")
        result_count=$((result_count - 1))  # 减去标题行
        echo "保留异质性变异数 (AF>1%): $result_count" | tee -a "$log_file"

        # 显示前几个变异作为示例
        if [ "$result_count" -gt 0 ]; then
            echo "前5个变异位点示例:" | tee -a "$log_file"
            head -6 "$result_file" | tee -a "$log_file"
        fi

        # 可选：保存过滤后的VCF文件供后续使用
        cp "$filtered_temp_vcf.gz" "$chloro_vcf" 2>/dev/null
        if [ -f "$filtered_temp_vcf.gz.tbi" ]; then
            cp "$filtered_temp_vcf.gz.tbi" "$chloro_vcf.tbi" 2>/dev/null
        fi
        echo "最终过滤的VCF已保存: $chloro_vcf" | tee -a "$log_file"
    fi

    # 清理临时文件
    rm -f "$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam"
    rm -f "$gatk_vcf_gz" "$gatk_vcf_gz.tbi" 2>/dev/null
    rm -f "$temp_vcf" 2>/dev/null
    rm -f "$filtered_temp_vcf.gz" "$filtered_temp_vcf.gz.tbi" 2>/dev/null

    echo "样本 $sample 处理完成" | tee -a "$log_file"
done
```
* 汇总
```
PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross_new"

# 需要排除的样本列表
EXCLUDED_SAMPLES=("Sample_c1c2" "Sample_l2c2" "Sample_l2l3" "Sample_l4c1" "Sample_l4l3")

# 提取样本名（去除前后空格和空行）
samples=$(find "$PROJECT_ROOT" -name "chloroplast_heteroplasmy_freq.txt" -type f | xargs -I {} dirname {} | xargs -I {} basename {} | sort -u)

# 创建汇总文件
summary_file="$PROJECT_ROOT/chloroplast_heteroplasmy_summary.txt"
echo -e "Sample\tChr\tPos\tRef\tAlt\tFrequency\tDepth" > "$summary_file"

# 循环处理已有结果的样本
for sample in $samples; do
    # 检查是否在排除列表中
    if [[ " ${EXCLUDED_SAMPLES[@]} " =~ " ${sample} " ]]; then
        echo "跳过排除样本: $sample"
        continue
    fi
    
    result_file="$PROJECT_ROOT/$sample/chloroplast_heteroplasmy_freq.txt"
    
    if [ ! -f "$result_file" ]; then
        echo "警告: 结果文件不存在，跳过样本 $sample"
        continue
    fi
    
    # 检查文件是否为空或只有标题
    file_size=$(wc -l < "$result_file")
    if [ "$file_size" -le 1 ]; then
        echo "样本 $sample 无变异数据"
        continue
    fi
    
    echo "处理样本: $sample, 变异数: $((file_size - 1))"
    
    # 提取数据行（跳过标题行），并添加样本名
    tail -n +2 "$result_file" | \
    awk -v sample="$sample" 'BEGIN{OFS="\t"} {
        # 重命名列以匹配输出格式
        # 输入格式: Chr Pos Ref Alt AF DP QUAL
        # 输出格式: Sample Chr Pos Ref Alt Frequency Depth
        print sample, $1, $2, $3, $4, $5, $6
    }' >> "$summary_file"
done




PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross_new"


# 父母本样本
MOTHER="Sample_Col_G"
FATHER="Sample_Ler_XL_4"

# 需要排除的样本列表（包括父母本）
EXCLUDED_SAMPLES=("Sample_c1c2" "Sample_l2c2" "Sample_l2l3" "Sample_l4c1" "Sample_l4l3" "$MOTHER" "$FATHER")

# 叶绿体基因组总大小（拟南芥叶绿体约154,478 bp）
CHLORO_GENOME_SIZE=154478

# 汇总文件
inheritance_analysis_file="$PROJECT_ROOT/chloroplast_inheritance_analysis_genome.txt"
detailed_inheritance_file="$PROJECT_ROOT/chloroplast_detailed_inheritance_genome.txt"

# 步骤1：读取父母本的变异位点
echo "读取父母本变异位点..."
declare -A mother_variants
declare -A father_variants

# 读取母本变异
mother_file="$PROJECT_ROOT/$MOTHER/chloroplast_heteroplasmy_freq.txt"
if [ -f "$mother_file" ]; then
    echo "读取母本文件: $mother_file"
    while IFS=$'\t' read -r chr pos ref alt freq depth qual; do
        if [[ "$chr" != "Chr" && -n "$chr" ]]; then
            key="${chr}:${pos}:${ref}:${alt}"
            mother_variants["$key"]=1
        fi
    done < "$mother_file"
    echo "母本 $MOTHER 变异数: ${#mother_variants[@]}"
else
    echo "错误: 找不到母本文件 $mother_file"
    exit 1
fi

# 读取父本变异
father_file="$PROJECT_ROOT/$FATHER/chloroplast_heteroplasmy_freq.txt"
if [ -f "$father_file" ]; then
    echo "读取父本文件: $father_file"
    while IFS=$'\t' read -r chr pos ref alt freq depth qual; do
        if [[ "$chr" != "Chr" && -n "$chr" ]]; then
            key="${chr}:${pos}:${ref}:${alt}"
            father_variants["$key"]=1
        fi
    done < "$father_file"
    echo "父本 $FATHER 变异数: ${#father_variants[@]}"
else
    echo "错误: 找不到父本文件 $father_file"
    exit 1
fi

# 步骤2：创建详细遗传分析文件
echo -e "Progeny\tChr\tPos\tRef\tAlt\tAF\tSource" > "$detailed_inheritance_file"
echo -e "Progeny\tTotal_Sites\tFather_Sites\tNovel_Sites\tMother_Sites\tMother_Rate\tFather_Rate\tNovel_Rate" > "$inheritance_analysis_file"

# 步骤3：直接从VCF文件分析每个子代样本
for progeny_dir in "$PROJECT_ROOT"/Sample_*/; do
    progeny=$(basename "$progeny_dir")
    
    # 跳过排除样本和父母本
    if [[ " ${EXCLUDED_SAMPLES[@]} " =~ " ${progeny} " ]] || [[ "$progeny" == "$MOTHER" ]] || [[ "$progeny" == "$FATHER" ]]; then
        echo "跳过样本: $progeny"
        continue
    fi
    
    vcf_file="$progeny_dir/3_bwa/chloroplast_variants_from_gatk.vcf.gz"
    
    # 初始化位点计数器
    total_sites=0
    father_sites=0
    novel_sites=0
    
    if [ ! -f "$vcf_file" ]; then
        echo "警告: 找不到VCF文件，样本 $progeny 无变异数据"
        # 没有VCF文件，说明没有检测到变异，100%来自母本
        total_sites=0
        father_sites=0
        novel_sites=0
    else
        echo "分析子代: $progeny"
        echo "VCF文件: $vcf_file"
        
        # 检查VCF文件是否有内容
        variant_count=$(bcftools view "$vcf_file" 2>/dev/null | grep -v "^#" | wc -l)
        echo "  VCF中变异总数: $variant_count"
        
        # 使用临时文件避免子shell问题
        temp_file=$(mktemp)
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' "$vcf_file" 2>/dev/null > "$temp_file"
        
        # 从临时文件读取，避免子shell问题
        while IFS=$'\t' read -r chr pos ref alt af; do
            if [[ -n "$chr" && -n "$af" ]]; then
                # 检查AF是否大于0.01
                if (( $(echo "$af > 0.01" | bc -l 2>/dev/null) )); then
                    total_sites=$((total_sites + 1))
                    
                    key="${chr}:${pos}:${ref}:${alt}"
                    source=""
                    
                    # 检查变异来源
                    if [[ -n "${father_variants[$key]}" ]]; then
                        father_sites=$((father_sites + 1))
                        source="Father"
                    else
                        novel_sites=$((novel_sites + 1))
                        source="Novel"
                    fi
                    
                    # 记录详细信息
                    echo -e "$progeny\t$chr\t$pos\t$ref\t$alt\t$af\t$source" >> "$detailed_inheritance_file"
                fi
            fi
        done < "$temp_file"
        
        # 清理临时文件
        rm -f "$temp_file"
    fi
    
    # 计算占基因组比例
    if [ $CHLORO_GENOME_SIZE -gt 0 ]; then
        if [ $total_sites -eq 0 ]; then
            # 没有检测到变异，100%来自母本
            mother_rate=100.0000
            father_rate=0.0000
            novel_rate=0.0000
            mother_sites=$CHLORO_GENOME_SIZE
        else
            # 有变异的情况
            father_rate=$(echo "scale=6; $father_sites / $CHLORO_GENOME_SIZE * 100" | bc)
            novel_rate=$(echo "scale=6; $novel_sites / $CHLORO_GENOME_SIZE * 100" | bc)
            # 母本比例 = 100% - 父本比例 - 新突变比例
            mother_rate=$(echo "scale=6; 100 - $father_rate - $novel_rate" | bc)
            # 母本位点数（理论值）
            mother_sites=$(echo "scale=0; $CHLORO_GENOME_SIZE - $father_sites - $novel_sites" | bc)
        fi
    else
        mother_rate=0
        father_rate=0
        novel_rate=0
        mother_sites=0
    fi
    
    # 输出汇总结果
    echo -e "$progeny\t$total_sites\t$father_sites\t$novel_sites\t$mother_sites\t${mother_rate}%\t${father_rate}%\t${novel_rate}%" >> "$inheritance_analysis_file"
    
    if [ $total_sites -eq 0 ]; then
        echo "  - 无变异检测，100%来自母本"
    else
        echo "  - 变异位点数: $total_sites/$CHLORO_GENOME_SIZE"
        echo "  - 来源: 母本 ${mother_rate}%, 父本 ${father_rate}%, 新突变 ${novel_rate}%"
    fi
done

# 步骤4：生成总体统计
echo -e "\n===== 叶绿体遗传分析完成 (占基因组比例) ====="
echo "叶绿体基因组大小: $CHLORO_GENOME_SIZE bp"
echo "汇总文件: $inheritance_analysis_file"
echo "详细文件: $detailed_inheritance_file"

# 统计无变异的样本数
no_variant_count=$(awk -F'\t' 'NR>1 && $2 == 0 {count++} END {print count+0}' "$inheritance_analysis_file")
total_progeny_count=$(awk -F'\t' 'NR>1 {count++} END {print count+0}' "$inheritance_analysis_file")

echo -e "\n样本统计:"
echo "总子代样本数: $total_progeny_count"
echo "无变异样本数: $no_variant_count (100%母本)"
echo "有变异样本数: $((total_progeny_count - no_variant_count))"

# 显示总体统计
echo -e "\n总体遗传模式统计 (占基因组比例):"
awk -F'\t' -v genome_size=$CHLORO_GENOME_SIZE 'NR>1 {
    total_s += $2; father_s += $3; novel_s += $4; mother_s += $5
} END {
    if(genome_size > 0) {
        total_p = total_s/(NR-1)/genome_size*100
        mother_p = mother_s/(NR-1)/genome_size*100
        father_p = father_s/(NR-1)/genome_size*100  
        novel_p = novel_s/(NR-1)/genome_size*100
        printf "平均每个样本:\n"
        printf "总变异位点: %.2f (%.4f%%)\n", total_s/(NR-1), total_p
        printf "来自母本: %.2f (%.4f%%)\n", mother_s/(NR-1), mother_p
        printf "来自父本: %.2f (%.4f%%)\n", father_s/(NR-1), father_p
        printf "新突变: %.2f (%.4f%%)\n", novel_s/(NR-1), novel_p
    } else {
        printf "基因组大小未定义\n"
    }
}' "$inheritance_analysis_file"

# 显示前几个子代的结果
echo -e "\n前10个子代样本结果:"
head -11 "$inheritance_analysis_file" | column -t


PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross_new"

echo "=== 将TXT汇总文件转换为CSV格式 ==="

# 转换遗传分析汇总文件
if [ -f "$PROJECT_ROOT/chloroplast_inheritance_analysis_genome.txt" ]; then
    echo "转换 chloroplast_inheritance_analysis_genome.txt -> .csv"
    # 将制表符替换为逗号，并转换百分比格式
    sed 's/\t/,/g; s/ %/%/g' "$PROJECT_ROOT/chloroplast_inheritance_analysis_genome.txt" > "$PROJECT_ROOT/chloroplast_inheritance_analysis_genome.csv"
    echo "完成: chloroplast_inheritance_analysis_genome.csv"
else
    echo "警告: 找不到 chloroplast_inheritance_analysis_genome.txt"
fi

# 转换详细遗传分析文件
if [ -f "$PROJECT_ROOT/chloroplast_detailed_inheritance_genome.txt" ]; then
    echo "转换 chloroplast_detailed_inheritance_genome.txt -> .csv"
    sed 's/\t/,/g' "$PROJECT_ROOT/chloroplast_detailed_inheritance_genome.txt" > "$PROJECT_ROOT/chloroplast_detailed_inheritance_genome.csv"
    echo "完成: chloroplast_detailed_inheritance_genome.csv"
else
    echo "警告: 找不到 chloroplast_detailed_inheritance_genome.txt"
fi


# 转换主汇总文件
if [ -f "$PROJECT_ROOT/chloroplast_heteroplasmy_summary.txt" ]; then
    echo "转换 chloroplast_heteroplasmy_summary.txt -> .csv"
    sed 's/\t/,/g' "$PROJECT_ROOT/chloroplast_heteroplasmy_summary.txt" > "$PROJECT_ROOT/chloroplast_heteroplasmy_summary.csv"
    echo "完成: chloroplast_heteroplasmy_summary.csv"
fi
```
# compare fold0
```
awk -F'\t' '
NR==FNR && FNR>1 {
    key = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5
    fold0[key] = $0
    next
}
NR>FNR && FNR>1 {
    key = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5
    new[key] = $0
}
END {
    for (key in fold0) {
        if (!(key in new)) {
            print fold0[key]
        }
    }
}' Atha_cross_fold0_new/chloroplast_heteroplasmy_summary.txt Atha_cross_new/chloroplast_heteroplasmy_summary.txt > fold0_unique_full.tsv
```
## plot
```
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(scales)

# 读取数据
df <- read_csv("chloroplast_inheritance_analysis_genome.csv")

# 转换百分比列为数值型
df <- df %>%
  mutate(
    Mother_Rate = as.numeric(gsub("%", "", Mother_Rate)),
    Father_Rate = as.numeric(gsub("%", "", Father_Rate)),
    Novel_Rate = as.numeric(gsub("%", "", Novel_Rate))
  )

# 转换数据为长格式
df_long <- df %>%
  select(Progeny, Mother_Rate, Father_Rate, Novel_Rate) %>%
  pivot_longer(
    cols = c(Mother_Rate, Father_Rate, Novel_Rate), 
    names_to = "Source", 
    values_to = "Percentage"
  ) %>%
  filter(!is.na(Percentage))

# 清理来源名称并设置顺序
df_long$Source <- recode(df_long$Source,
                        Mother_Rate = "Mother",
                        Father_Rate = "Father", 
                        Novel_Rate = "Novel")

# 设置因子水平以控制顺序
df_long$Source <- factor(df_long$Source, levels = c("Mother", "Father", "Novel"))

# 对数坐标图 - 英文标签和正确顺序
p_log <- ggplot(df_long, aes(x = Source, y = Percentage + 0.0001, color = Source)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2, color = "red") +
  scale_y_log10(
    labels = function(x) paste0(x, "%"),
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100)
  ) +
  scale_color_manual(values = c("Mother" = "#1f77b4", "Father" = "#ff7f0e", "Novel" = "#2ca02c")) +
  labs(x = "Source", y = "Genome Proportion (%)") +  # 英文标签
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )

# 保存图片
ggsave("chloroplast_inheritance_dotplot_log.png", p_log, width = 6, height = 4, dpi = 300)

# 显示图片
print(p_log)

# 父本和新突变图 - 英文标签
df_non_mother <- df_long %>% filter(Source != "Mother")

# 设置因子水平以控制顺序（只保留Father和Novel）
df_non_mother$Source <- factor(df_non_mother$Source, levels = c("Father", "Novel"))

p_non_mother <- ggplot(df_non_mother, aes(x = Source, y = Percentage, color = Source)) +
  geom_point(size = 1, alpha = 0.6, position = position_jitter(width = 0.1, height = 0)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2, color = "red") +
  scale_color_manual(values = c("Father" = "#ff7f0e", "Novel" = "#2ca02c")) +
  labs(x = "Source", y = "Genome Proportion (%)") +  # 英文标签
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )

ggsave("chloroplast_father_novel_only.png", p_non_mother, width = 6, height = 4, dpi = 300)

# 显示图片
print(p_non_mother)
```
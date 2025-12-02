## pure col reads
```
cd ~/zxy/plastid
mkdir -p artificial_heteroplasmy
cd artificial_heteroplasmy

# 提取Col叶绿体reads
samtools view -b -F 0x904 -q 1 ../mix3/SRR611086_2/3_bwa/R.sort.bam Pt > col_pt_clean.bam
samtools index col_pt_clean.bam

# 统计信息
col_reads=$(samtools view -c col_pt_clean.bam)
echo "Col叶绿体有效reads数: $col_reads"
#4760821

# 转为FASTQ
samtools fastq col_pt_clean.bam > col_pt.fq

# 统计reads数
total_reads=$(wc -l < col_pt.fq | awk '{print $1/4}')
echo "总reads数: $total_reads"
```

# Pt:1000 1% mut
```
cd ~/zxy/plastid/artificial_heteroplasmy

echo "检查Pt:1000位点的覆盖情况..."
samtools depth -r Pt:1000-1000 col_pt_clean.bam

# 或者查看整个区域的覆盖
echo "查看Pt:990-1010区域的覆盖:"
samtools depth -r Pt:990-1010 col_pt_clean.bam

samtools view col_pt_clean.bam Pt:1000-1000 | awk '{print $1}' > covering_reads.txt
covering_count=$(wc -l < covering_reads.txt)
echo "覆盖Pt:1000的reads数: $covering_count"
#4854

mutate_count=$(echo "scale=0; $covering_count * 0.0025 / 1" | bc)
echo "需要突变的reads数: $mutate_count"
#485

# 将BAM转为SAM进行修改
echo "将BAM转为SAM进行修改..."
samtools view -h col_pt_clean.bam > col_pt_clean.sam

# 随机选择要突变的reads
samtools view col_pt_clean.bam Pt:1000-1000 | awk '{print $1}' | shuf -n $mutate_count > to_mutate.txt

# 创建简单的awk脚本来修改SAM
cat > sam_format.awk << 'EOF'
BEGIN {
    # 读取要突变的reads
    while ((getline < "to_mutate.txt") > 0) {
        mutate[$1] = 1
    }
    target_pos = 1000
    OFS = "\t"  # 设置输出字段分隔符为制表符
}

/^@/ {
    print
    next
}

{
    # 确保使用制表符分隔
    if (NF < 11) {
        # 如果字段数不够，可能是分隔符问题，重新分割
        split($0, fields, "\t")
        $0 = ""
        for (i=1; i<=length(fields); i++) {
            $i = fields[i]
        }
    }
    
    # 现在应该有正确的字段
    if (NF >= 11) {
        # 修改突变reads
        if ($1 in mutate) {
            start_pos = $4 + 0
            seq = $10
            read_pos = target_pos - start_pos + 1
            
            if (read_pos >= 1 && read_pos <= length(seq)) {
                current_base = substr(seq, read_pos, 1)
                if (current_base == "c" || current_base == "C") {
                    $10 = substr(seq, 1, read_pos-1) "T" substr(seq, read_pos+1)
                }
            }
        }
        
        # 输出，确保用制表符分隔
        printf "%s", $1
        for (i=2; i<=NF; i++) {
            printf "\t%s", $i
        }
        printf "\n"
    }
}
EOF

# 生成正确的SAM文件
samtools view -H col_pt_clean.bam > header_fixed.sam
grep -v "^@" col_pt_clean.sam | awk -f sam_format.awk >> header_fixed.sam

# 验证格式
echo "验证SAM文件格式..."
head -28935 header_fixed.sam | tail -5 | cat -A  # 显示制表符和换行符

samtools view -b header_fixed.sam > artificial_mut_fixed.bam
samtools sort artificial_mut_fixed.bam -o artificial_mut_fixed_sorted.bam
samtools index artificial_mut_fixed_sorted.bam

# 检查BAM文件
echo "检查BAM文件状态..."
samtools view -c artificial_mut_fixed_sorted.bam Pt:1000-1000


samtools mpileup -r Pt:1000-1000 artificial_mut_fixed_sorted.bam | \
    awk '{
        bases = $5
        total = length(bases)
        print "总碱基数: " total
        
        # 分别统计
        a_count = gsub(/[Aa]/, "", bases)
        c_count = gsub(/[Cc]/, "", bases) 
        g_count = gsub(/[Gg]/, "", bases)
        t_count = gsub(/[Tt]/, "", bases)
        other_count = gsub(/[^ACGTacgt]/, "", bases)
        
        print "A: " a_count
        print "C: " c_count
        print "G: " g_count
        print "T: " t_count
        print "其他字符: " other_count
        
        # 计算比例
        if (total > 0) {
            print "T碱基比例: " (t_count/total*100) "%"
            print "预期比例: 10%"
        }
    }'

gatk --java-options "-Xmx40g" \
    Mutect2 \
    --native-pair-hmm-threads 16 \
    -R ../../mix2/genome/genome.fa \
    -I artificial_mut_fixed_sorted.bam \
    -L Pt \
    --mitochondria-mode \
    --max-reads-per-alignment-start 100 \
    --max-mnp-distance 0 \
    -O mix_Pt.raw.vcf



gatk --java-options "-Xmx80g" \
    FilterMutectCalls \
    -R ../../mix2/genome/genome.fa \
    -V mix_Pt.raw.vcf \
    --max-alt-allele-count 4 \
    --mitochondria-mode \
    -O mix.filtered.vcf

# 过滤出PASS的SNP
bcftools view --apply-filters PASS --types snps --max-alleles 2 --targets Pt --include "AF>=0.01" -Oz mix.filtered.vcf > mix_Pt.vcf.gz
gunzip mix_Pt.vcf.gz
```
# summary
```
echo -e "突变频率\t突变位点\t总reads数\t突变reads数\t是否检测到\t检测AF\t过滤状态\tTLOD\t深度" > all_results.tsv

# 逐个处理每个文件夹
for folder in 0.1percent 0.2percent 0.25percent 0.3percent 0.5percent 1percent 5percent 10percent; do
    echo "处理 $folder..."
    
    # 计算突变率
    case $folder in
        "0.1percent") rate=0.001 ;;
        "0.2percent") rate=0.002 ;;
        "0.25percent") rate=0.0025;;
        "0.3percent") rate=0.003 ;;
        "0.5percent") rate=0.005 ;;
        "1percent") rate=0.01 ;;
        "5percent") rate=0.05 ;;
        "10percent") rate=0.10 ;;
    esac
    
    total_reads=4854
    mutated_reads=$(echo "$total_reads * $rate / 1" | bc)
    
    if [ -f "$folder/mix_Pt.vcf" ]; then
        # 检查Pt:1000位置
        vcf_info=$(bcftools view --targets Pt:1000-1000 "$folder/mix_Pt.vcf" 2>/dev/null | grep -v "^#")
        
        if [ -n "$vcf_info" ]; then
            # 解析信息
            chrom=$(echo "$vcf_info" | awk '{print $1}')
            pos=$(echo "$vcf_info" | awk '{print $2}')
            filter=$(echo "$vcf_info" | awk '{print $7}')
            
            # 提取INFO中的DP
            info_dp=$(echo "$vcf_info" | grep -o 'DP=[0-9]*' | head -1 | cut -d= -f2)
            
            # 提取样本中的DP（第4个冒号分隔的值）
            sample_dp=$(echo "$vcf_info" | awk '{split($10,a,":"); print a[4]}')
            
            # 提取TLOD
            tlod=$(echo "$vcf_info" | grep -o 'TLOD=[0-9.]*' | head -1 | cut -d= -f2)
            
            # 提取AF
            af=$(echo "$vcf_info" | awk '{split($10,a,":"); print a[3]}')
            
            if [ "$filter" = "PASS" ]; then
                detected="YES"
            else
                detected="FILTERED"
            fi
            
            # 使用样本DP
            echo -e "$rate\tPt:1000\t$total_reads\t$mutated_reads\t$detected\t$af\t$filter\t$tlod\t$sample_dp" >> all_results.tsv
        else
            echo -e "$rate\tPt:1000\t$total_reads\t$mutated_reads\tNO\t0\tNA\tNA\tNA" >> all_results.tsv
        fi
    else
        echo -e "$rate\tPt:1000\t$total_reads\t$mutated_reads\tNO\t0\tNA\tNA\tNA" >> all_results.tsv
    fi
done
```

#!/bin/bash

# ==== 基础配置 ====
PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross"
OPTS_FILE="$PROJECT_ROOT/opts.tsv"
OUTPUT_PREFIX="$PROJECT_ROOT/chloro_multi_stats"

CSV_OUTPUT="${OUTPUT_PREFIX}.csv"
MD_OUTPUT="${OUTPUT_PREFIX}.md"

# 要排除的样本列表
EXCLUDED_SAMPLES=("Sample_c1c2" "Sample_l2c2" "Sample_l2l3" "Sample_l4c1" "Sample_l4l3" "")

# ==== 提取有效样本 ====
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

# ==== 生成 CSV ====
echo "样本名,有效变异数,高频(≥10%),中高频(5%-10%),低频(1%-5%),平均异质性频率(%),最高异质性频率(%),平均总reads数(×),最高频率位点(Pt:POS)" > "$CSV_OUTPUT"

for sample in "${samples[@]}"; do
    result_file="$PROJECT_ROOT/$sample/chloroplast_heteroplasmy_freq.txt"
    
    if [ ! -f "$result_file" ] || [ ! -s "$result_file" ]; then
        echo "样本 $sample 无效或文件为空，跳过"
        continue
    fi

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

echo "CSV 文件已生成：$CSV_OUTPUT"

# ==== 生成 Markdown 报告 ====
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

echo "Markdown 报告已生成：$MD_OUTPUT"

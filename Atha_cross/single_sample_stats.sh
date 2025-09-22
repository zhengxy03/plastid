#!/bin/bash
set -euo pipefail

# ======================================
# 单个样本叶绿体异质性变异统计脚本
# 使用方法：./single_sample_stats.sh 样本结果文件路径
# 示例：./single_sample_stats.sh /share/home/wangq/zxy/plastid/Atha_cross/Sample_14/chloroplast_heteroplasmy_freq.txt
# ======================================

# 1. 检查输入参数
if [ $# -ne 1 ]; then
    echo "错误：请传入单个样本的结果文件路径！"
    echo "使用方法：./single_sample_stats.sh 样本结果文件路径"
    exit 1
fi

result_file="$1"
sample_name=$(basename $(dirname "$result_file"))  # 从路径中提取样本名（如 Sample_14）
stats_report="${result_file%.txt}_stats.txt"       # 统计报告输出路径（如 ..._stats.txt）

# 2. 检查结果文件是否存在
if [ ! -f "$result_file" ]; then
    echo "错误：结果文件 $result_file 不存在！"
    exit 1
fi

# 3. 检查结果文件是否为空
if [ ! -s "$result_file" ]; then
    echo "警告：样本 $sample_name 的结果文件为空（未检测到有效变异）！"
    exit 0
fi

# 4. 核心统计计算（基于结果文件的6列：CHROM POS REF ALT FREQ TOTAL）
echo "正在统计样本 $sample_name 的变异特征..."

# 4.1 基础指标
total_valid_vars=$(wc -l < "$result_file")  # 有效变异总数（已过滤alt_cnt≥1）
avg_total_reads=$(awk '{sum+=$6} END{printf "%.1f", sum/NR}' "$result_file")  # 平均总reads数
max_freq=$(awk '{max=$5>max?$5:max} END{printf "%.6f", max}' "$result_file")  # 最高异质性频率
max_freq_pos=$(awk -v max=0 '$5>max{max=$5; pos=$2} END{print pos}' "$result_file")  # 最高频率位点

# 4.2 频率区间分类（按1%、5%、10%划分）
freq_0_1=$(awk '$5>0 && $5<0.01 {count++} END{print count+0}' "$result_file")  # 0%-1%（极低频）
freq_1_5=$(awk '$5>=0.01 && $5<0.05 {count++} END{print count+0}' "$result_file")  # 1%-5%（低频）
freq_5_10=$(awk '$5>=0.05 && $5<0.1 {count++} END{print count+0}' "$result_file")  # 5%-10%（中高频）
freq_10_plus=$(awk '$5>=0.1 {count++} END{print count+0}' "$result_file")  # ≥10%（高频）

# 4.3 平均异质性频率
avg_freq=$(awk '{sum+=$5} END{printf "%.6f", sum/NR}' "$result_file")

# 5. 生成统计报告
cat > "$stats_report" << EOF
# 样本 $sample_name 叶绿体异质性变异统计报告
# 生成时间：$(date "+%Y-%m-%d %H:%M:%S")
# 原始结果文件：$result_file

==================================== 核心指标 ====================================
1. 有效变异总数（alt_cnt≥1）：$total_valid_vars 个
2. 平均总reads数（覆盖度）：$avg_total_reads ×
3. 平均异质性频率：$avg_freq （约 $(printf "%.2f%%" $(echo "$avg_freq*100" | bc))）
4. 最高异质性频率：$max_freq （约 $(printf "%.2f%%" $(echo "$max_freq*100" | bc))）
5. 最高频率位点：Pt:$max_freq_pos

==================================== 频率分布 ====================================
- 极低频变异（0% < 频率 < 1%）：$freq_0_1 个
- 低频变异（1% ≤ 频率 < 5%）：$freq_1_5 个
- 中高频变异（5% ≤ 频率 < 10%）：$freq_5_10 个
- 高频变异（频率 ≥ 10%）：$freq_10_plus 个

==================================== 数据可信度 ==================================
- 所有变异均满足：替代碱基支持数 ≥ 1
- 平均覆盖度 $avg_total_reads ×（覆盖度≥30×时，频率计算可信度高）
EOF

# 6. 输出结果
echo "统计完成！"
echo "样本 $sample_name 核心统计："
echo "→ 有效变异数：$total_valid_vars 个"
echo "→ 平均频率：$(printf "%.2f%%" $(echo "$avg_freq*100" | bc))"
echo "→ 高频变异数（≥10%）：$freq_10_plus 个"
echo "详细报告已保存至：$stats_report"
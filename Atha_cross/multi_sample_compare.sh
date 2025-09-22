#!/bin/bash
set -euo pipefail

# ======================================
# 多样本叶绿体异质性统计脚本（双格式输出）
# 功能：1. 从./opts.tsv读取样本（排除指定样本）；2. 生成CSV+MD双格式结果
# 排除样本：Sample_c1c2、Sample_l2c2、Sample_l2l3、Sample_l4c1、Sample_l4l3
# 使用方法：./multi_sample_dual_format.sh 项目根目录 输出前缀
# 示例：./multi_sample_dual_format.sh /share/home/wangq/zxy/plastid/Atha_cross ./chloro_multi_stats
# 输出文件：./chloro_multi_stats.csv + ./chloro_multi_stats.md
# ======================================

# 1. 基础配置
EXCLUDED_SAMPLES=("Sample_c1c2" "Sample_l2c2" "Sample_l2l3" "Sample_l4c1" "Sample_l4l3")
OPTS_FILE="./opts.tsv"  # opts.tsv路径（当前目录）
project_root="$1"       # 项目根目录（含所有Sample文件夹）
output_prefix="$2"      # 输出文件前缀（如 ./chloro_multi_stats）
csv_output="${output_prefix}.csv"  # CSV输出路径
md_output="${output_prefix}.md"    # MD输出路径
valid_samples=()        # 存储过滤后的有效样本（格式：样本名:结果文件路径）

# 2. 输入参数检查
if [ $# -ne 2 ]; then
    echo "错误：参数数量不足！"
    echo "使用方法：./multi_sample_dual_format.sh 项目根目录 输出前缀"
    echo "示例：./multi_sample_dual_format.sh /share/home/wangq/zxy/plastid/Atha_cross ./chloro_multi_stats"
    exit 1
fi

# 3. 依赖文件检查
# 3.1 检查opts.tsv是否存在
if [ ! -f "$OPTS_FILE" ]; then
    echo "错误：./opts.tsv 文件不存在！请确保opts.tsv在当前目录。"
    exit 1
fi
# 3.2 检查项目根目录是否存在
if [ ! -d "$project_root" ]; then
    echo "错误：项目根目录 $project_root 不存在！"
    exit 1
fi

# 4. 从opts.tsv筛选有效样本
echo "=== 从 $OPTS_FILE 筛选有效样本 ==="
# 读取opts.tsv（第一列为样本名，跳过空行和注释行）
while IFS= read -r line; do
    # 跳过空行或#开头的注释行
    if [ -z "$line" ] || [[ "$line" =~ ^# ]]; then
        continue
    fi
    # 提取第一列作为样本名（去除前后空格）
    sample_name=$(echo "$line" | awk -F'\t' '{print $1}' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//')
    # 检查样本是否在排除列表
    is_excluded=0
    for excluded in "${EXCLUDED_SAMPLES[@]}"; do
        if [ "$sample_name" = "$excluded" ]; then
            is_excluded=1
            echo "→ 排除样本：$sample_name（在排除列表中）"
            break
        fi
    done
    if [ $is_excluded -eq 1 ]; then
        continue
    fi
    # 检查样本文件夹和结果文件是否存在
    sample_dir="${project_root}/${sample_name}"
    result_file="${sample_dir}/chloroplast_heteroplasmy_freq.txt"
    if [ ! -d "$sample_dir" ]; then
        echo "→ 跳过样本：$sample_name（样本文件夹 $sample_dir 不存在）"
        continue
    fi
    if [ ! -f "$result_file" ]; then
        echo "→ 跳过样本：$sample_name（结果文件 $result_file 不存在）"
        continue
    fi
    if [ ! -s "$result_file" ]; then
        echo "→ 跳过样本：$sample_name（结果文件为空）"
        continue
    fi
    # 加入有效样本列表
    valid_samples+=("${sample_name}:${result_file}")
    echo "→ 有效样本：$sample_name（结果文件存在且非空）"
done < "$OPTS_FILE"

# 5. 有效样本数量检查
if [ ${#valid_samples[@]} -eq 0 ]; then
    echo "错误：未找到有效样本！请检查opts.tsv和样本结果文件。"
    exit 1
fi
echo "=== 筛选完成：共 ${#valid_samples[@]} 个有效样本 ==="

# 6. 生成CSV文件（便于数据处理）
echo -e "\n=== 生成 CSV 对比表 ==="
# CSV表头
echo "样本名,有效变异数,高频变异数(≥10%),中高频变异数(5%-10%),低频变异数(1%-5%),平均异质性频率(%),最高异质性频率(%),平均总reads数(×),最高频率位点(Pt:POS)" > "$csv_output"

# 批量计算每个样本指标并写入CSV
for sample_info in "${valid_samples[@]}"; do
    # 拆分样本名和结果文件路径
    sample_name=$(echo "$sample_info" | cut -d':' -f1)
    result_file=$(echo "$sample_info" | cut -d':' -f2)
    echo "→ 处理样本：$sample_name"

    # 核心指标计算
    total_valid_vars=$(wc -l < "$result_file")  # 有效变异总数
    # 频率区间计数
    freq_1_5=$(awk '$5>=0.01 && $5<0.05 {c++} END{print c+0}' "$result_file")  # 1%-5%
    freq_5_10=$(awk '$5>=0.05 && $5<0.1 {c++} END{print c+0}' "$result_file")   # 5%-10%
    freq_10_plus=$(awk '$5>=0.1 {c++} END{print c+0}' "$result_file")           # ≥10%
    # 平均/最高频率（转换为百分比）
    avg_freq=$(awk '{sum+=$5} END{printf "%.2f", sum/NR*100}' "$result_file")   # 平均频率（%）
    max_freq=$(awk '{max=$5>max?$5:max} END{printf "%.2f", max*100}' "$result_file")  # 最高频率（%）
    # 平均覆盖度
    avg_total_reads=$(awk '{sum+=$6} END{printf "%.1f", sum/NR}' "$result_file")  # 平均总reads数
    # 最高频率位点
    max_freq_pos=$(awk -v max=0 '$5>max{max=$5; pos=$2} END{print pos}' "$result_file")

    # 写入CSV行
    echo "${sample_name},${total_valid_vars},${freq_10_plus},${freq_5_10},${freq_1_5},${avg_freq},${max_freq},${avg_total_reads},${max_freq_pos}" >> "$csv_output"
done
echo "→ CSV文件已保存至：$csv_output"

# 7. 生成MD文件（便于报告展示）
echo -e "\n=== 生成 MD 报告 ==="
cat > "$md_output" << EOF
# 叶绿体异质性变异多样本统计报告
**生成时间**：$(date "+%Y-%m-%d %H:%M:%S")  
**项目根目录**：${project_root}  
**opts.tsv路径**：${OPTS_FILE}  
**排除样本列表**：${EXCLUDED_SAMPLES[*]}  
**有效样本数量**：${#valid_samples[@]} 个  
**关联文件**：[CSV格式对比表](${csv_output#./})（可用于Excel进一步分析）

---

## 一、统计说明
1. **有效变异定义**：替代碱基支持数 ≥ 1（已过滤无意义假阳性）；
2. **频率区间划分**：
   - 低频变异：1% ≤ 频率 < 5%  
   - 中高频变异：5% ≤ 频率 < 10%  
   - 高频变异：频率 ≥ 10%  
3. **覆盖度评估**：平均总reads数反映测序深度，≥30×时频率计算可信度高。

---

## 二、多样本核心指标对比表
| 样本名       | 有效变异数 | 高频变异数（≥10%） | 中高频变异数（5%-10%） | 低频变异数（1%-5%） | 平均异质性频率（%） | 最高异质性频率（%） | 平均总reads数（×） | 最高频率位点（Pt:POS） |
|--------------|------------|--------------------|------------------------|--------------------|---------------------|---------------------|--------------------|------------------------|
EOF

# 批量写入MD表格内容
for sample_info in "${valid_samples[@]}"; do
    sample_name=$(echo "$sample_info" | cut -d':' -f1)
    result_file=$(echo "$sample_info" | cut -d':' -f2)

    # 重新计算指标（与CSV一致，避免变量复用问题）
    total_valid_vars=$(wc -l < "$result_file")
    freq_1_5=$(awk '$5>=0.01 && $5<0.05 {c++} END{print c+0}' "$result_file")
    freq_5_10=$(awk '$5>=0.05 && $5<0.1 {c++} END{print c+0}' "$result_file")
    freq_10_plus=$(awk '$5>=0.1 {c++} END{print c+0}' "$result_file")
    avg_freq=$(awk '{sum+=$5} END{printf "%.2f", sum/NR*100}' "$result_file")
    max_freq=$(awk '{max=$5>max?$5:max} END{printf "%.2f", max*100}' "$result_file")
    avg_total_reads=$(awk '{sum+=$6} END{printf "%.1f", sum/NR}' "$result_file")
    max_freq_pos=$(awk -v max=0 '$5>max{max=$5; pos=$2} END{print pos}' "$result_file")

    # 写入MD表格行
    echo "| ${sample_name} | ${total_valid_vars} | ${freq_10_plus} | ${freq_5_10} | ${freq_1_5} | ${avg_freq} | ${max_freq} | ${avg_total_reads} | Pt:${max_freq_pos} |" >> "$md_output"
done

# MD报告结尾补充说明
cat >> "$md_output" << EOF

---

## 三、关键结论提示
1. **高频变异样本**：可重点关注「高频变异数≥3」的样本（可能存在显著异质性）；
2. **数据可信度**：「平均总reads数<30×」的样本建议重新评估测序质量；
3. **后续分析**：可基于CSV文件筛选「最高异质性频率≥15%」的位点，进行Sanger验证。

---

*报告生成脚本：multi_sample_dual_format.sh*
EOF

echo "→ MD报告已保存至：$md_output"
echo -e "\n=== 所有文件生成完成 ==="
echo "1. CSV对比表：$csv_output"
echo "2. MD报告：$md_output"
#!/bin/bash
set -euo pipefail

# ====================== 核心配置（无需修改，已适配你的项目） ======================
PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross"  # 项目根目录
CHLORO_NAME="Pt"                                         # 叶绿体染色体名称
OPTS_FILE="$PROJECT_ROOT/opts.tsv"                       # 样本列表文件
EXCLUDE_MIXED=("Sample_l2c2" "Sample_l2l3" "Sample_l4c1" "Sample_l4l3")  # 排除的4个混样
TARGET_BAM="chloroplast_Pt_sorted.bam"                  # 目标叶绿体BAM文件名
TARGET_BAM_BAI="${TARGET_BAM}.bai"                      # 目标BAM索引文件名
# ==================================================================================

# 函数1：检查文件是否存在
check_file() {
    local file_path=$1
    local sample_name=$2
    local file_desc=$3
    if [ ! -f "$file_path" ]; then
        echo -e "\n❌ 错误：样本 $sample_name 的【$file_desc】不存在 - $file_path"
        echo "提示：请先确保该样本已完成前期BWA比对（生成R.sort.bam）"
        exit 1
    fi
}

# 函数2：判断样本是否为混样（是则返回1，否则返回0）
is_mixed_sample() {
    local sample=$1
    for mixed in "${EXCLUDE_MIXED[@]}"; do
        if [ "$sample" = "$mixed" ]; then
            return 1  # 是混样，排除
        fi
    done
    return 0  # 非混样，保留
}

# 函数3：检查目标BAM是否已存在（是则返回1，否则返回0）
is_bam_exist() {
    local sample=$1
    local bam_path="$PROJECT_ROOT/$sample/3_bwa/$TARGET_BAM"
    local bai_path="$PROJECT_ROOT/$sample/3_bwa/$TARGET_BAM_BAI"
    if [ -f "$bam_path" ] && [ -f "$bai_path" ] && [ -s "$bam_path" ]; then
        return 1  # BAM和索引均存在且非空，跳过
    fi
    return 0  # BAM缺失或不完整，需要生成
}

# ====================== 主流程 ======================
echo -e "=================================================="
echo "===== 批量生成所有样本的叶绿体BAM及索引 ====="
echo "=================================================="
echo "项目根目录：$PROJECT_ROOT"
echo "排除的混样：${EXCLUDE_MIXED[*]}"
echo "目标BAM文件名：$TARGET_BAM"
echo -e "==================================================\n"

# 步骤1：验证基础文件
if [ ! -f "$OPTS_FILE" ]; then
    echo "❌ 错误：样本列表文件 opts.tsv 不存在 - $OPTS_FILE"
    exit 1
fi
if [ ! -d "$PROJECT_ROOT" ]; then
    echo "❌ 错误：项目根目录不存在 - $PROJECT_ROOT"
    exit 1
fi

# 步骤2：提取所有有效样本（排除混样）
echo "【步骤1：提取有效样本】"
VALID_SAMPLES=()
while read -r sample _ _; do
    # 跳过空行和注释行
    if [ -z "$sample" ] || [[ "$sample" =~ ^# ]]; then
        continue
    fi
    # 排除混样，保留有效样本
    if is_mixed_sample "$sample"; then
        VALID_SAMPLES+=("$sample")
    else
        echo "ℹ️  排除混样：$sample"
    fi
done < "$OPTS_FILE"

# 检查有效样本数量
if [ ${#VALID_SAMPLES[@]} -eq 0 ]; then
    echo -e "\n❌ 错误：未从 opts.tsv 提取到有效样本"
    exit 1
fi
echo -e "✅ 共提取到 ${#VALID_SAMPLES[@]} 个有效样本："
printf "  - %s\n" "${VALID_SAMPLES[@]}"
echo -e ""

# 步骤3：批量生成叶绿体BAM（循环处理每个有效样本）
echo "【步骤2：批量生成叶绿体BAM及索引】"
for sample in "${VALID_SAMPLES[@]}"; do
    echo -e "--------------------------------------------------"
    echo "处理样本：$sample"
    
    # 定义样本路径
    bam_dir="$PROJECT_ROOT/$sample/3_bwa"
    raw_bam="$bam_dir/R.sort.bam"          # 前期生成的原始BAM
    raw_bam_bai="$bam_dir/R.sort.bai"             # 原始BAM索引
    temp_bam="$bam_dir/chloroplast_temp.bam"  # 临时叶绿体BAM
    target_bam="$bam_dir/$TARGET_BAM"      # 最终目标BAM
    target_bam_bai="$bam_dir/$TARGET_BAM_BAI"  # 最终目标BAM索引
    log_file="$PROJECT_ROOT/$sample/chloro_bam_generate.log"  # 生成日志
    
    # 检查1：原始BAM是否存在（确保前期比对已完成）
    check_file "$raw_bam" "$sample" "原始BAM（R.sort.bam）"
    check_file "$raw_bam_bai" "$sample" "原始BAM索引（R.sort.bai）"
    
    # 检查2：目标BAM是否已存在（已存在则跳过）
    if is_bam_exist "$sample"; then
        echo "⚠️  目标BAM（$TARGET_BAM）已存在，跳过样本 $sample"
        continue
    fi
    
    # 检查3：任务是否已在运行（避免重复提交）
    if bjobs -w | grep -q "${sample}_chloro_bam"; then
        echo "⚠️  生成BAM的任务（${sample}_chloro_bam）已在运行，跳过样本 $sample"
        continue
    fi
    
    # 清理旧的临时文件（避免残留）
    rm -f "$temp_bam" "$target_bam" "$target_bam_bai"
    > "$log_file"  # 初始化日志
    
    # 提交LSF任务：提取叶绿体区域、排序、生成索引
    bsub -q mpi -n 24 -J "${sample}_chloro_bam" << EOF
        echo "===== 开始生成样本 $sample 的叶绿体BAM =====" | tee -a "$log_file"
        echo "开始时间：\$(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$log_file"
        echo "原始BAM路径：$raw_bam" | tee -a "$log_file"
        echo "目标BAM路径：$target_bam" | tee -a "$log_file"
        
        # 子步骤1：提取叶绿体区域并过滤低质量reads（MAPQ<10，排除次级比对）
        echo -e "\n1. 提取叶绿体区域并过滤低质量reads..." | tee -a "$log_file"
        # -q 10：过滤MAPQ<10的低质量比对
        # -F 2304：排除次级比对（FLAG=2048）和补充比对（FLAG=256）
        samtools view -b -q 10 -F 2304 -o "$temp_bam" "$raw_bam" "$CHLORO_NAME"
        
        # 检查临时BAM是否有效（非空）
        if [ ! -s "$temp_bam" ]; then
            echo "❌ 错误：提取的叶绿体临时BAM为空（可能无叶绿体reads）" | tee -a "$log_file"
            rm -f "$temp_bam"
            exit 1
        fi
        echo "✅ 叶绿体临时BAM生成成功（大小：\$(du -sh "$temp_bam" | cut -f1)）" | tee -a "$log_file"
        
        # 子步骤2：对临时BAM排序（按坐标排序，确保后续分析可用）
        echo -e "\n2. 对叶绿体BAM进行坐标排序..." | tee -a "$log_file"
        samtools sort -@ 8 -o "$target_bam" "$temp_bam"
        # 检查排序后的BAM是否有效
        if [ ! -s "$target_bam" ]; then
            echo "❌ 错误：叶绿体BAM排序失败（输出文件为空）" | tee -a "$log_file"
            rm -f "$temp_bam" "$target_bam"
            exit 1
        fi
        
        # 子步骤3：生成BAM索引（必须，否则后续无法按位点提取reads）
        echo -e "\n3. 生成叶绿体BAM索引..." | tee -a "$log_file"
        samtools index "$target_bam" "$target_bam_bai"
        # 检查索引是否有效
        if [ ! -f "$target_bam_bai" ]; then
            echo "❌ 错误：叶绿体BAM索引生成失败" | tee -a "$log_file"
            rm -f "$temp_bam" "$target_bam" "$target_bam_bai"
            exit 1
        fi
        
        # 子步骤4：清理临时文件并记录结果
        echo -e "\n4. 清理临时文件并确认结果..." | tee -a "$log_file"
        rm -f "$temp_bam"  # 仅保留最终BAM和索引
        # 最终验证
        if [ -f "$target_bam" ] && [ -f "$target_bam_bai" ]; then
            echo "✅ 样本 $sample 的叶绿体BAM及索引生成完成！" | tee -a "$log_file"
            echo "   - BAM大小：\$(du -sh "$target_bam" | cut -f1)" | tee -a "$log_file"
            echo "   - 索引大小：\$(du -sh "$target_bam_bai" | cut -f1)" | tee -a "$log_file"
        else
            echo "❌ 错误：叶绿体BAM或索引生成后缺失" | tee -a "$log_file"
            exit 1
        fi
        
        echo -e "\n===== 样本 $sample 处理结束 =====" | tee -a "$log_file"
        echo "结束时间：\$(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$log_file"
EOF

    echo "✅ 样本 $sample 的BAM生成任务已提交（任务名：${sample}_chloro_bam）"
done

# 步骤4：提示后续操作
echo -e "\n=================================================="
echo "===== 所有样本的BAM生成任务提交完毕！ ====="
echo "=================================================="
echo "后续操作建议："
echo "1. 用 'bjobs' 命令查看任务进度（等待所有任务完成）"
echo "2. 任务完成后，运行 'trace_parentage_final.sh' 进行来源追溯分析"
echo "3. 若需验证BAM是否生成：ls -lh $PROJECT_ROOT/[样本名]/3_bwa/$TARGET_BAM"
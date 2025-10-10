#!/bin/bash
set -euo pipefail

# ====================== 核心配置（已适配你的样本信息，无需修改） ======================
PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross"  # 你的项目根目录
CHLORO_NAME="Pt"                                         # 叶绿体染色体名称
MOTHER="Sample_Col_G"                                    # 母本样本名（固定）
FATHER="Sample_Ler_XL_4"                                 # 父本样本名（固定）
# 排除的样本（母本+父本+4个混样，无需处理）
EXCLUDE_SAMPLES=("Sample_Col_G" "Sample_Ler_XL_4" "Sample_l2c2" "Sample_l2l3" "Sample_l4c1" "Sample_l4l3")
# ======================================================================================

# 函数1：检查文件是否存在（不存在则报错退出，确保依赖文件完整）
check_file() {
    local file_path=$1
    local sample_name=$2
    local file_desc=$3  # 文件描述（如“叶绿体BAM”“变异VCF”）
    if [ ! -f "$file_path" ]; then
        echo -e "\n❌ 错误：样本 $sample_name 的【$file_desc】不存在 - $file_path"
        echo "提示：请先确保该样本已生成叶绿体BAM（chloroplast_Pt_sorted.bam）和筛选后的VCF"
        exit 1
    fi
}

# 函数2：判断样本是否为有效子代（排除母本、父本、混样）
is_valid_offspring() {
    local sample=$1
    for exclude in "${EXCLUDE_SAMPLES[@]}"; do
        if [ "$sample" = "$exclude" ]; then
            return 1  # 是排除样本，返回1
        fi
    done
    return 0  # 是有效子代，返回0
}

# 函数3：解析reads的目标位点碱基（处理纯匹配CIGAR，适配叶绿体无复杂结构特点）
get_read_base() {
    local bam_path=$1
    local chrom=$2
    local pos=$3
    # 提取目标位点的有效reads（排除次级比对、低质量，仅保留纯匹配）
    samtools view -F 2304 -q 10 "$bam_path" "$chrom:$pos" | \
    awk -v target_pos="$pos" '{
        # 仅处理CIGAR为纯匹配（如“150M”）的reads，避免插入缺失干扰
        if ($6 ~ /^[0-9]+M$/) {
            # 计算reads上目标位点的偏移量（start为reads起始位置）
            offset = target_pos - $4 + 1;
            # 确保偏移量在reads长度范围内，提取对应碱基
            if (offset >= 1 && offset <= length($10)) {
                print substr($10, offset, 1);
            }
        }
    }'
}

# ====================== 主分析流程 ======================
echo -e "=================================================="
echo "===== 叶绿体子代reads父母本来源追溯分析 ====="
echo "=================================================="
echo "父本样本：$FATHER"
echo "母本样本：$MOTHER"
echo "项目根目录：$PROJECT_ROOT"
echo -e "==================================================\n"

# 步骤1：检查父本和母本的核心文件（确保父母本已完成前期分析）
echo "【步骤1：检查父本/母本核心文件】"
# 父本核心文件（筛选后的叶绿体变异VCF及索引）
FATHER_VCF="$PROJECT_ROOT/$FATHER/3_bwa/chloroplast_variants_from_gatk.vcf.gz"
FATHER_VCF_TBI="$FATHER_VCF.tbi"
# 母本核心文件
MOTHER_VCF="$PROJECT_ROOT/$MOTHER/3_bwa/chloroplast_variants_from_gatk.vcf.gz"
MOTHER_VCF_TBI="$MOTHER_VCF.tbi"

# 验证父本文件
check_file "$FATHER_VCF" "$FATHER" "叶绿体变异VCF"
check_file "$FATHER_VCF_TBI" "$FATHER" "叶绿体变异VCF索引"
# 验证母本文件
check_file "$MOTHER_VCF" "$MOTHER" "叶绿体变异VCF"
check_file "$MOTHER_VCF_TBI" "$MOTHER" "叶绿体变异VCF索引"
echo -e "✅ 父本/母本核心文件均存在，继续分析\n"

# 步骤2：从opts.tsv提取所有有效子代样本
echo "【步骤2：提取有效子代样本】"
OFFSPRING_SAMPLES=()
# 读取opts.tsv第一列（样本名），过滤排除样本
while read -r sample _ _; do
    # 跳过空行和注释行
    if [ -z "$sample" ] || [[ "$sample" =~ ^# ]]; then
        continue
    fi
    # 判断是否为有效子代
    if is_valid_offspring "$sample"; then
        OFFSPRING_SAMPLES+=("$sample")
    else
        echo "ℹ️  排除非子代样本：$sample（母本/父本/混样）"
    fi
done < "$PROJECT_ROOT/opts.tsv"

# 检查是否有有效子代样本
if [ ${#OFFSPRING_SAMPLES[@]} -eq 0 ]; then
    echo -e "\n❌ 错误：未从opts.tsv提取到有效子代样本"
    exit 1
fi
echo -e "✅ 共提取到 ${#OFFSPRING_SAMPLES[@]} 个有效子代样本："
printf "  - %s\n" "${OFFSPRING_SAMPLES[@]}"
echo -e ""

# 步骤3：循环分析每个子代样本
for OFFSPRING in "${OFFSPRING_SAMPLES[@]}"; do
    echo -e "=================================================="
    echo "【步骤3：分析子代样本：$OFFSPRING】"
    echo "=================================================="
    
    # 定义子代样本的核心文件路径（与前期预处理输出路径一致）
    OFF_BAM="$PROJECT_ROOT/$OFFSPRING/3_bwa/chloroplast_Pt_sorted.bam"       # 保留的叶绿体BAM
    OFF_BAM_TBI="$OFF_BAM.bai"                                               # BAM索引
    OFF_VCF="$PROJECT_ROOT/$OFFSPRING/3_bwa/chloroplast_variants_from_gatk.vcf.gz"  # 筛选变异VCF
    OFF_VCF_TBI="$OFF_VCF.tbi"                                               # VCF索引
    OFF_FREQ="$PROJECT_ROOT/$OFFSPRING/chloroplast_heteroplasmy_freq.txt"    # 异质性频率文件
    # 输出结果文件（单独存放，避免覆盖）
    RESULT_FILE="$PROJECT_ROOT/$OFFSPRING/chloroplast_parentage_tracing.txt"
    LOG_FILE="$PROJECT_ROOT/$OFFSPRING/chloroplast_parentage_log.txt"        # 分析日志
    
    # 步骤3.1：检查子代样本的核心文件
    echo "1. 检查子代样本核心文件..."
    check_file "$OFF_BAM" "$OFFSPRING" "叶绿体BAM（chloroplast_Pt_sorted.bam）"
    check_file "$OFF_BAM_TBI" "$OFFSPRING" "叶绿体BAM索引"
    check_file "$OFF_VCF" "$OFFSPRING" "筛选后叶绿体变异VCF"
    check_file "$OFF_VCF_TBI" "$OFFSPRING" "筛选后叶绿体变异VCF索引"
    check_file "$OFF_FREQ" "$OFFSPRING" "异质性频率文件"
    echo "✅ 子代样本 $OFFSPRING 核心文件齐全"
    
    # 步骤3.2：初始化结果文件和日志
    echo -e "\n2. 初始化结果文件..."
    > "$LOG_FILE"  # 清空旧日志
    echo "分析开始时间：$(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$LOG_FILE"
    echo "子代样本：$OFFSPRING | 父本：$FATHER | 母本：$MOTHER" | tee -a "$LOG_FILE"
    
    # 结果文件表头（制表符分隔，支持Excel直接打开）
    echo -e "染色体\t位点\t参考碱基\t变异碱基\t子代变异频率\t总reads数\t父本基因型\t母本基因型\t父本来源reads数\t母本来源reads数\t父本来源占比\t母本来源占比\t来源判断" > "$RESULT_FILE"
    
    # 步骤3.3：提取子代变异位点并追溯来源
    echo -e "\n3. 提取变异位点并追溯来源（植物叶绿体默认母系遗传）..."
    # 用bcftools提取子代VCF的关键信息：染色体、位点、参考碱基、变异碱基
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" "$OFF_VCF" | while read -r CHROM POS REF ALT; do
        # 冗余检查：确保是叶绿体区域（避免非目标位点）
        if [ "$CHROM" != "$CHLORO_NAME" ]; then
            echo "⚠️  跳过非叶绿体位点：$CHROM:$POS" | tee -a "$LOG_FILE"
            continue
        fi
        
        # 3.3.1 从频率文件提取该位点的变异频率和总reads数
        OFF_INFO=$(awk -v c="$CHROM" -v p="$POS" -v r="$REF" -v a="$ALT" \
            '$1==c && $2==p && $3==r && $4==a {print $5,$6}' "$OFF_FREQ")
        if [ -z "$OFF_INFO" ]; then
            echo "⚠️  子代位点 $CHROM:$POS（$REF→$ALT）未在频率文件中找到，跳过" | tee -a "$LOG_FILE"
            continue
        fi
        OFF_FREQ_VAL=$(echo "$OFF_INFO" | cut -d' ' -f1)  # 变异频率
        TOTAL_READS=$(echo "$OFF_INFO" | cut -d' ' -f2)   # 总reads数
        
        # 3.3.2 提取父本/母本在该位点的基因型（无变异时默认纯合参考0/0）
        echo "   分析位点：$CHROM:$POS（参考：$REF，变异：$ALT）" | tee -a "$LOG_FILE"
        # 父本基因型（若父本无该变异，基因型为“.”，替换为0/0）
        FATHER_GT=$(bcftools query -f "%GT\n" -r "$CHROM:$POS" "$FATHER_VCF" 2>/dev/null || echo "0/0")
        FATHER_GT=${FATHER_GT:-.}  # 处理空值
        if [ "$FATHER_GT" = "." ] || [ -z "$FATHER_GT" ]; then
            FATHER_GT="0/0"
        fi
        # 母本基因型（处理逻辑同父本）
        MOTHER_GT=$(bcftools query -f "%GT\n" -r "$CHROM:$POS" "$MOTHER_VCF" 2>/dev/null || echo "0/0")
        MOTHER_GT=${MOTHER_GT:-.}
        if [ "$MOTHER_GT" = "." ] || [ -z "$MOTHER_GT" ]; then
            MOTHER_GT="0/0"
        fi
        
        # 3.3.3 提取子代该位点的有效reads碱基
        READS_BASES=$(get_read_base "$OFF_BAM" "$CHROM" "$POS")
        # 统计有效reads总数（去重前）
        READS_COUNT=$(echo "$READS_BASES" | wc -w | tr -d ' ')
        echo "   该位点有效reads数：$READS_COUNT" | tee -a "$LOG_FILE"
        
        # 3.3.4 转换父母本基因型为实际碱基（0=参考碱基，1=变异碱基）
        # 父本碱基：0/0→REF，1/1→ALT，0/1→REF+ALT
        case "$FATHER_GT" in
            "0/0") FATHER_BASES="$REF" ;;
            "1/1") FATHER_BASES="$ALT" ;;
            "0/1"|"1/0") FATHER_BASES="$REF $ALT" ;;
            *) FATHER_BASES="$REF" ;;  # 异常基因型默认参考碱基
        esac
        # 母本碱基（逻辑同父本）
        case "$MOTHER_GT" in
            "0/0") MOTHER_BASES="$REF" ;;
            "1/1") MOTHER_BASES="$ALT" ;;
            "0/1"|"1/0") MOTHER_BASES="$REF $ALT" ;;
            *) MOTHER_BASES="$REF" ;;
        esac
        
        # 3.3.5 统计父本/母本来源的reads数（优先匹配母本，符合叶绿体母系遗传）
        FATHER_READS=0
        MOTHER_READS=0
        for BASE in $READS_BASES; do
            # 先匹配母本碱基（避免碱基重叠时误判）
            if echo "$MOTHER_BASES" | grep -q "\b$BASE\b"; then
                ((MOTHER_READS++))
            # 再匹配父本碱基
            elif echo "$FATHER_BASES" | grep -q "\b$BASE\b"; then
                ((FATHER_READS++))
            fi
        done
        
        # 3.3.6 计算来源占比并判断遗传模式
        if [ $((FATHER_READS + MOTHER_READS)) -eq 0 ]; then
            FATHER_RATIO="0.00"
            MOTHER_RATIO="0.00"
            SOURCE="无法判断（无有效reads）"
        else
            # 用bc计算占比（保留2位小数）
            FATHER_RATIO=$(printf "%.2f" $(echo "scale=4; $FATHER_READS / ($FATHER_READS + $MOTHER_READS)" | bc))
            MOTHER_RATIO=$(printf "%.2f" $(echo "scale=4; $MOTHER_READS / ($FATHER_READS + $MOTHER_READS)" | bc))
            # 判断来源（母本占比>0.5为正常母系遗传）
            if [ $(echo "$MOTHER_RATIO > 0.5" | bc) -eq 1 ]; then
                SOURCE="母本（占比$MOTHER_RATIO，符合母系遗传）"
            elif [ $(echo "$FATHER_RATIO > 0.5" | bc) -eq 1 ]; then
                SOURCE="父本（占比$FATHER_RATIO，异常）"
            else
                SOURCE="混合来源（父本$FATHER_RATIO/母本$MOTHER_RATIO，需排查）"
            fi
        fi
        
        # 3.3.7 写入结果文件
        echo -e "$CHROM\t$POS\t$REF\t$ALT\t$OFF_FREQ_VAL\t$TOTAL_READS\t$FATHER_GT\t$MOTHER_GT\t$FATHER_READS\t$MOTHER_READS\t$FATHER_RATIO\t$MOTHER_RATIO\t$SOURCE" >> "$RESULT_FILE"
    done
    
    # 步骤3.4：统计分析结果并写入日志
    echo -e "\n4. 统计分析结果..."
    # 统计总变异位点数、正常/异常位点数
    TOTAL_VARIANTS=$(grep -v "^染色体" "$RESULT_FILE" | wc -l | tr -d ' ')
    NORMAL_VARIANTS=$(grep "符合母系遗传" "$RESULT_FILE" | wc -l | tr -d ' ')
    ABNORMAL_VARIANTS=$(grep "异常" "$RESULT_FILE" | wc -l | tr -d ' ')
    UNKNOWN_VARIANTS=$(grep "无法判断" "$RESULT_FILE" | wc -l | tr -d ' ')
    
    # 写入日志
    echo -e "\n===== 分析结果统计 =====" | tee -a "$LOG_FILE"
    echo "总变异位点数：$TOTAL_VARIANTS" | tee -a "$LOG_FILE"
    echo "符合母系遗传的位点：$NORMAL_VARIANTS" | tee -a "$LOG_FILE"
    echo "父本来源异常位点：$ABNORMAL_VARIANTS" | tee -a "$LOG_FILE"
    echo "无法判断来源的位点：$UNKNOWN_VARIANTS" | tee -a "$LOG_FILE"
    echo "结果文件路径：$RESULT_FILE" | tee -a "$LOG_FILE"
    echo "分析结束时间：$(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$LOG_FILE"
    
    # 预览前5行结果（验证格式）
    echo -e "\n5. 结果预览（前5行，含表头）："
    head -n 6 "$RESULT_FILE" | column -t -s $'\t'
    echo -e "✅ 子代样本 $OFFSPRING 分析完成！\n"
done

# 步骤
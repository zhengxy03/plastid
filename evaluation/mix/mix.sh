#!/bin/bash
set -euo pipefail  # 严格模式，报错立即退出

# ==============================
# 全局参数配置（核心：复用Col100倍BAM）
# ==============================
COL_FOLD=100                          # Col固定100倍
LER_FOLD_LIST=(2 4 8 16)          # 要运行的Ler倍数列表
SEED=123                               # 随机种子（可重复）
WORK=~/zxy/plastid/mix2/mix_results    # 主结果目录
GENOME=~/zxy/plastid/mix2/genome/genome.fa  # 基因组文件
# 关键：复用已生成的Col100倍BAM（之前Col100x_Ler1x目录下的结果）
COL_100x_BAM=~/zxy/plastid/mix2/mix_results/Col100x_Ler1x/col_100x.sm.bam
COL_100x_BAM_INDEX=${COL_100x_BAM}.bai  # Col100倍BAM索引
LER_BAM=~/zxy/plastid/mix2/SRR616965_2/3_bwa/R.sort.bam  # Ler原始BAM
# 复用已生成的位点文件（无需重新生成）
COL_SITES=~/zxy/plastid/mix2/mix_results/Col100x_Ler1x/Col_sites.txt
LER_UNIQUE_SITES=~/zxy/plastid/mix2/mix_results/Col100x_Ler1x/Ler_unique_sites.txt

# ==============================
# 预处理：检查复用文件是否存在
# ==============================
echo "=== 开始批量分析：Col${COL_FOLD}x（复用） + Ler[${LER_FOLD_LIST[*]}]x ==="
echo "依赖文件检查..."
# 检查Col100倍BAM及索引（核心复用文件）
[ ! -f "$COL_100x_BAM" ] && echo "错误：Col100倍BAM文件不存在！路径：$COL_100x_BAM" && exit 1
[ ! -f "$COL_100x_BAM_INDEX" ] && echo "错误：Col100倍BAM索引不存在！" && exit 1
# 检查其他依赖文件
[ ! -f "$LER_BAM" ] && echo "错误：Ler原始BAM文件不存在！" && exit 1
[ ! -f "$GENOME" ] && echo "错误：基因组文件不存在！" && exit 1
[ ! -f "$COL_SITES" ] && echo "错误：Col位点文件不存在！" && exit 1
[ ! -f "$LER_UNIQUE_SITES" ] && echo "错误：Ler特有位点文件不存在！" && exit 1

# 生成基因组字典（若未生成）
[ ! -f "${GENOME%.fa}.dict" ] && {
  echo "生成基因组序列字典..."
  cd "$(dirname "$GENOME")"
  gatk CreateSequenceDictionary -R "$GENOME" -O "${GENOME%.fa}.dict"
}

echo "复用Col100倍BAM：$COL_100x_BAM"
echo "Col100倍BAM有效reads数：$(samtools view -c "$COL_100x_BAM")"

# ==============================
# 循环遍历每个Ler倍数，批量运行（仅处理Ler和后续步骤）
# ==============================
for LER_FOLD in "${LER_FOLD_LIST[@]}"; do
  # 生成当前倍数的结果目录（例如：Col100x_Ler0.5x）
  RESULT_DIR="${WORK}/Col${COL_FOLD}x_Ler${LER_FOLD}x"
  mkdir -p "$RESULT_DIR"
  echo -e "\n=================================================="
  echo "正在处理：Col${COL_FOLD}x（复用） + Ler${LER_FOLD}x → 结果目录：$RESULT_DIR"
  echo "=================================================="
  cd "$RESULT_DIR" || exit 1

  # ==============================
  # 步骤1：仅处理Ler 目标倍数（0.5/2/4/8/16倍）
  # ==============================
  echo -e "\n[1/3] 生成Ler${LER_FOLD}x Pt区域BAM..."
  # 提取Ler Pt区域有效reads（过滤低质量）
  samtools view -b -F 0x904 -q 1 "$LER_BAM" Pt > ler_pt_raw.bam
  samtools index ler_pt_raw.bam
  ler_pt_raw_reads=$(samtools view -c ler_pt_raw.bam)
  echo "Ler原始Pt区域有效reads：$ler_pt_raw_reads"

  # 根据倍数处理：
  if [[ "$LER_FOLD" == "0.5" ]]; then
    # 0.5倍 = 随机抽取50%的reads（可重复抽样）
    samtools view -s ${SEED}.5 -b ler_pt_raw.bam -o ler_${LER_FOLD}x_temp.bam
    echo "Ler 0.5倍：随机抽取50% reads"
  else
    # ≥2倍 = 复制对应次数，添加唯一RG ID
    for ((i=1; i<=LER_FOLD; i++)); do
      samtools addreplacerg -r "@RG\tID:Ler_$i\tSM:Ler_${LER_FOLD}x\tPL:ILLUMINA\tLB:Ler\tPU:Ler_$i" \
        ler_pt_raw.bam -o ler_pt_copy_$i.bam
    done
    # 合并复制后的BAM
    samtools cat ler_pt_copy_*.bam -o ler_${LER_FOLD}x_merged.bam
    samtools sort ler_${LER_FOLD}x_merged.bam -o ler_${LER_FOLD}x_temp.bam
  fi

  # 统一添加最终RG ID、排序、建索引
  samtools addreplacerg -r "@RG\tID:Ler_Final\tSM:Ler_${LER_FOLD}x\tPL:ILLUMINA\tLB:Ler\tPU:Ler_Final" \
    ler_${LER_FOLD}x_temp.bam -o ler_${LER_FOLD}x.sm.bam
  samtools index ler_${LER_FOLD}x.sm.bam

  # 验证Ler倍数效果（修复：用bc处理浮点数计算）
  ler_target_reads=$(samtools view -c ler_${LER_FOLD}x.sm.bam)
  echo "Ler${LER_FOLD}x Pt区域有效reads：$ler_target_reads"
  # 用bc计算理论目标reads（兼容浮点数和整数）
  ler_theory_reads=$(echo "$ler_pt_raw_reads * $LER_FOLD" | bc | xargs printf "%.0f")  # 四舍五入为整数
  echo "Ler理论目标reads：$ler_theory_reads"

  # ==============================
  # 步骤2：混合样本（复用Col100倍BAM） + GATK变异检测
  # ==============================
  echo -e "\n[2/3] 混合样本（复用Col100x） + GATK变异检测..."
  # 直接复用Col100倍BAM，仅合并Ler目标倍数BAM
  mix_bam="mix_Col${COL_FOLD}x_Ler${LER_FOLD}x.bam"
  samtools merge -c "${mix_bam%.bam}_temp.bam" \
    "$COL_100x_BAM" \
    ler_${LER_FOLD}x.sm.bam
  samtools sort "${mix_bam%.bam}_temp.bam" -o "$mix_bam"
  samtools index "$mix_bam"
  mix_reads=$(samtools view -c "$mix_bam")
  echo "混合样本Pt区域有效reads：$mix_reads"

  # Mutect2变异检测
  gatk --java-options "-Xmx40g" \
    Mutect2 \
    --native-pair-hmm-threads 16 \
    -R "$GENOME" \
    -I "$mix_bam" \
    -L Pt \
    --mitochondria-mode \
    --max-reads-per-alignment-start 100 \
    --max-mnp-distance 0 \
    -O mix_Pt.raw.vcf

  # FilterMutectCalls过滤
  gatk --java-options "-Xmx80g" \
    FilterMutectCalls \
    -R "$GENOME" \
    -V mix_Pt.raw.vcf \
    --max-alt-allele-count 4 \
    --mitochondria-mode \
    -O mix.filtered.vcf

  # 提取PASS、单ALT的Pt区域变异
  bcftools view --apply-filters PASS --max-alleles 2 --targets Pt -O v mix.filtered.vcf > mix_Pt.vcf
  echo "变异检测完成：生成 mix_Pt.vcf"

  # ==============================
  # 步骤3：检测率 + 异质性频率计算（与之前逻辑一致）
  # ==============================
  echo -e "\n[3/3] 计算检测率 + 异质性频率..."
  mix_vcf="mix_Pt.vcf"
  total_ler=$(wc -l < "$LER_UNIQUE_SITES")

  # 计算Ler特有位点检测率
  detected_all=$(grep -Ff "$LER_UNIQUE_SITES" "$mix_vcf" | wc -l)
  grep -Ff "$LER_UNIQUE_SITES" "$mix_vcf" | awk '!/#/ {print $1 "\t" $2}' | sort > detected_ler_sites.tmp
  detected_true=$(comm -23 detected_ler_sites.tmp <(sort "$COL_SITES") | wc -l)
  detection_rate_all=$(echo "scale=2; $detected_all / $total_ler * 100" | bc)
  detection_rate_true=$(echo "scale=2; $detected_true / $total_ler * 100" | bc)

  # 输出检测率结果
  echo -e "\nLer特有位点总数：$total_ler"
  echo -e "混合样本检测到的Ler位点（含假阳性）：$detected_all → 检测率：$detection_rate_all%"
  echo -e "排除Col假阳性后，真实Ler检测位点：$detected_true → 真实检测率：$detection_rate_true%"

  # 提取真实Ler位点列表
  comm -23 detected_ler_sites.tmp <(sort "$COL_SITES") > true_ler_sites_final.txt

  # 从VCF中提取真实Ler位点的完整数据
  awk '
  BEGIN {
    while ((getline < "true_ler_sites_final.txt") > 0) {
      key = $1 "\t" $2
      ler_sites[key] = 1
    }
  }
  !/#/ {
    key = $1 "\t" $2
    if (ler_sites[key] == 1) {
      print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10 "\t" $11
    }
  }
  ' "$mix_vcf" > true_ler_vcf_data.txt
  echo "从VCF中提取到的真实Ler位点数据行数：$(wc -l < true_ler_vcf_data.txt)"

  # 生成异质性频率TSV（仅核心字段）
  echo -e "CHROM\tPOS\tREF\tLer_ALT\tCol_REF_reads\tLer_ALT_reads\t混合总reads\t异质性频率(%)" > ler_heteroplasmy_final.tsv

  # 循环计算异质性频率
  while IFS=$'\t' read -r chrom pos ref alt col_sample ler_sample; do
    [[ -z "$chrom" || -z "$pos" ]] && continue
    # 提取reads数并转换为数字
    col_ref_reads=$(echo "$col_sample" | cut -d':' -f2 | cut -d',' -f1 | awk '{print $1+0}')
    ler_alt_reads=$(echo "$ler_sample" | cut -d':' -f2 | cut -d',' -f2 | awk '{print $1+0}')
    total_reads=$((col_ref_reads + ler_alt_reads))
    # 计算异质性频率（避免除零）
    if [[ $total_reads -eq 0 || $ler_alt_reads -eq 0 ]]; then
      hetero_freq="0.0000"
    else
      hetero_freq=$(echo "scale=4; $ler_alt_reads / $total_reads * 100" | bc 2>/dev/null)
    fi
    # 写入结果
    echo -e "$chrom\t$pos\t$ref\t$alt\t$col_ref_reads\t$ler_alt_reads\t$total_reads\t$hetero_freq" >> ler_heteroplasmy_final.tsv
  done < true_ler_vcf_data.txt

  # 统计核心指标
  valid_sites=$(awk -F'\t' 'NR>1 && $8!="0.0000" {count++} END {print count+0}' ler_heteroplasmy_final.tsv)
  avg_hetero=$(awk -F'\t' 'NR>1 {sum+=$8} END {if(NR>1) print sum/(NR-1); else print 0}' ler_heteroplasmy_final.tsv | bc -l | xargs printf "%.4f")

  # 输出当前倍数的最终报告
  echo -e "\n========================================"
  echo "Col${COL_FOLD}x（复用） + Ler${LER_FOLD}x 最终统计"
  echo "========================================"
  echo "1. 检测率：真实检测率 = $detection_rate_true%"
  echo "2. 异质性统计：有效位点 = $valid_sites | 平均异质性频率 = $avg_hetero%"
  echo "3. 关键文件："
  echo "   - 混合BAM：$mix_bam"
  echo "   - 变异VCF：mix_Pt.vcf"
  echo "   - 异质性数据：ler_heteroplasmy_final.tsv"
  echo "========================================\n"

  # 清理临时文件（保留核心结果）
  rm -f ler_pt_raw.bam ler_pt_copy_*.bam ler_${LER_FOLD}x_merged.bam \
        ler_${LER_FOLD}x_temp.bam "${mix_bam%.bam}_temp.bam" \
        detected_ler_sites.tmp true_ler_vcf_data.txt
done

echo -e "\n=== 所有倍数分析完成！ ==="
echo "结果目录汇总："
for LER_FOLD in "${LER_FOLD_LIST[@]}"; do
  echo " - Col${COL_FOLD}x_Ler${LER_FOLD}x：${WORK}/Col${COL_FOLD}x_Ler${LER_FOLD}x"
done
echo "核心复用说明：Col100倍BAM未重复生成，直接引用已有文件，节省大量时间！"
echo "每个目录下的 ler_heteroplasmy_final.tsv 可用于后续倍数-异质性关联分析"
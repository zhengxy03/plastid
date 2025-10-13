#!/bin/bash
set -euo pipefail

echo "===== 开始计算异质性父母本来源 (直接用 heteroplasmy_freq.txt) ====="

WORKDIR="/share/home/wangq/zxy/plastid/Atha_cross"
MOTHER_FREQ="${WORKDIR}/Sample_Col_G/chloroplast_heteroplasmy_freq.txt"
FATHER_FREQ="${WORKDIR}/Sample_Ler_XL_4/chloroplast_heteroplasmy_freq.txt"

check_file() {
    local f=$1
    if [[ ! -f "$f" ]]; then
        echo "❌ 缺少文件: $f"
        exit 1
    fi
}

check_file "$MOTHER_FREQ"
check_file "$FATHER_FREQ"

# 读取父母本频率文件到关联数组
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

    if [[ ! -s "$FREQ_FILE" ]]; then
        echo "⚠️ 缺少频率文件，跳过 ${SAMPLE}"
        continue
    fi

    echo -e "染色体\t位点\t参考碱基\t变异碱基\t异质性频率\t总reads数\t母本频率\t父本频率\t来源判断" > "$OUT_FILE"

    tail -n +1 "$FREQ_FILE" | while read -r chr pos ref alt hetfreq depth; do
        key="${chr}_${pos}_${alt}"
        mom_f=${mom_map["$key"]:-0}
        dad_f=${dad_map["$key"]:-0}

        # 判断来源
        if (( $(echo "$mom_f>0" | bc -l) )) && (( $(echo "$dad_f==0" | bc -l) )); then
            origin="母本偏向"
        elif (( $(echo "$dad_f>0" | bc -l) )) && (( $(echo "$mom_f==0" | bc -l) )); then
            origin="父本偏向"
        elif (( $(echo "$mom_f>0" | bc -l) )) && (( $(echo "$dad_f>0" | bc -l) )); then
            origin="混合偏向（父母均含）"
        else
            origin="自发变异（父母均无）"
        fi

        echo -e "${chr}\t${pos}\t${ref}\t${alt}\t${hetfreq}\t${depth}\t${mom_f}\t${dad_f}\t${origin}" >> "$OUT_FILE"
    done

    echo "✅ 已完成 ${SAMPLE} -> $(basename "$OUT_FILE")"
done

echo "===== 全部样本处理完成 ====="

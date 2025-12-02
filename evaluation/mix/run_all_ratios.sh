#!/bin/bash
WORK=~/zxy/plastid/mix3/mix_results
mkdir -p "$WORK"
cd "$WORK"

# 设置参数
SEED=123
GENOME=~/zxy/plastid/mix2/genome/genome.fa
COL_BAM=~/zxy/plastid/mix3/SRR611086_2/3_bwa/R.sort.bam
LER_BAM=~/zxy/plastid/mix3/SRR611087_2/3_bwa/R.sort.bam

# Ler比例列表
LER_RATIOS=("0.5" "2" "4" "8" "16" "32" "64")

for ratio in "${LER_RATIOS[@]}"; do
    echo "=========================================="
    echo "处理 Col100x_Ler${ratio}x"
    echo "=========================================="
    
    # 创建结果目录
    RESULT_DIR=$WORK/Col100x_Ler${ratio}x
    mkdir -p $RESULT_DIR
    cd $RESULT_DIR
    
    # 运行mix流程（只修改比例）
    echo "运行mix流程..."
    
    # Col - 直接使用原始BAM
    samtools view -b -F 0x904 -q 1 $COL_BAM Pt > col_pt.bam
    samtools index col_pt.bam
    col_reads=$(samtools view -c col_pt.bam)
    echo "Col 原始 Pt 区域有效 reads：$col_reads"
    
    # 计算目标Ler reads数（根据比例调整）
    if [[ $ratio == "0.5" ]]; then
        target_ler_reads=$((col_reads / 200))  # 0.5% = 1/200
    else
        target_ler_reads=$(echo "scale=0; $col_reads * $ratio / 100" | bc)
    fi
    echo "目标Ler reads数（Col的${ratio}%）：$target_ler_reads"
    
    # Ler - 采样到目标深度
    echo "处理Ler样本..."
    samtools view -b -F 0x904 -q 1 $LER_BAM Pt > ler_pt_full.bam
    ler_full_reads=$(samtools view -c ler_pt_full.bam)
    echo "Ler 原始 Pt 区域有效 reads：$ler_full_reads"
    
    if [[ $ler_full_reads -gt 0 ]]; then
        # 计算目标read pairs数（因为paired-end，所以要除以2）
        target_pairs=$((target_ler_reads / 2))
        echo "目标read pairs数：$target_pairs"
        
        # ========== 直接导出所有reads为单端，再按总量采样 ==========
        echo "导出Ler BAM所有reads为单端FASTQ..."
        # 导出所有reads（不区分R1/R2，避免配对过滤）
        samtools fastq ler_pt_full.bam > ler_all.fq
        # 计算目标总reads数（单端，无需除以2）
        target_total_reads=$target_ler_reads
        echo "目标总reads数：$target_total_reads"
        
        # 检查原始FASTQ的reads数
        original_total_reads=$(wc -l < ler_all.fq | awk '{print int($1/4)}')
        echo "原始总reads数：$original_total_reads"
        
        # 对单端FASTQ采样（直接取目标数量）
        echo "使用seqtk对单端reads采样..."
        seqtk sample -s $SEED ler_all.fq $target_total_reads > ler_sampled_single.fq
        
        # 检查采样结果
        sampled_total_reads=$(wc -l < ler_sampled_single.fq | awk '{print int($1/4)}')
        echo "采样后总reads数：$sampled_total_reads"
        
        # ========== 用单端模式比对（避免配对错误） ==========
        echo "将采样后的FASTQ转回BAM（单端模式）..."
        bwa mem -t 8 -p $GENOME ler_sampled_single.fq > ler_sampled.sam
        
        # 检查比对结果
        if [[ ! -s ler_sampled.sam ]]; then
            echo "错误：比对后SAM文件为空"
        fi
        
        # 转换SAM为BAM并过滤（保持原过滤条件）
        samtools view -b -F 0x904 -q 1 ler_sampled.sam > ler_sampled.bam
        samtools sort ler_sampled.bam -o ler_sampled_sorted.bam
        samtools index ler_sampled_sorted.bam
        
        # 清理临时文件
        rm -f ler_sampled.sam ler_sampled.bam ler_all.fq ler_sampled_single.fq
        
    else
        echo "错误：Ler样本没有reads"
        continue
    fi
    
    ler_sampled_reads=$(samtools view -c ler_sampled_sorted.bam)
    echo "采样后Ler Pt 区域有效 reads：$ler_sampled_reads"
    
    # 添加read group信息
    samtools addreplacerg -r "@RG\tID:Ler\tSM:Ler_${ratio}percent\tPL:ILLUMINA\tLB:Ler\tPU:Ler" \
        ler_sampled_sorted.bam -o ler_final.bam
    samtools index ler_final.bam
    
    samtools addreplacerg -r "@RG\tID:Col\tSM:Col_original\tPL:ILLUMINA\tLB:Col\tPU:Col" \
        col_pt.bam -o col_final.bam
    samtools index col_final.bam
    
    # 合并样本
    echo "合并Col和Ler样本..."
    samtools merge -c mix_Col100x_Ler${ratio}x_temp.bam \
        col_final.bam \
        ler_final.bam
        
    samtools sort mix_Col100x_Ler${ratio}x_temp.bam -o mix_Col100x_Ler${ratio}x.bam
    samtools index mix_Col100x_Ler${ratio}x.bam
    
    # 最终统计
    mix_reads=$(samtools view -c mix_Col100x_Ler${ratio}x.bam)
    col_final_reads=$(samtools view -c col_final.bam)
    ler_final_reads=$(samtools view -c ler_final.bam)
    
    echo "=== 最终统计 ==="
    echo "Col 原始深度 reads：$col_final_reads"
    echo "Ler ${ratio}%深度 reads：$ler_final_reads" 
    echo "混合样本总 reads：$mix_reads"
    echo "实际深度比例 Col:Ler = $col_final_reads : $ler_final_reads"
    
    # 计算实际比例
    if [[ $ler_final_reads -gt 0 ]]; then
        actual_ratio=$(echo "scale=2; $col_final_reads / $ler_final_reads" | bc)
        echo "实际比例 Col:Ler ≈ $actual_ratio:1"
        echo "目标比例 Col:Ler = 100:${ratio}"
        
        # 计算偏差
        ratio_diff=$(echo "scale=2; ($actual_ratio - 100) / 100 * 100" | bc)
        echo "与目标比例的偏差：${ratio_diff}%"
    fi
    
    # GATK变异检测
    echo "运行GATK变异检测..."
    gatk --java-options "-Xmx40g" \
        Mutect2 \
        --native-pair-hmm-threads 16 \
        -R $GENOME \
        -I mix_Col100x_Ler${ratio}x.bam \
        -L Pt \
        --mitochondria-mode \
        --max-reads-per-alignment-start 100 \
        --max-mnp-distance 0 \
        -O mix_Pt.raw.vcf
    
    gatk --java-options "-Xmx80g" \
        FilterMutectCalls \
        -R $GENOME \
        -V mix_Pt.raw.vcf \
        --max-alt-allele-count 4 \
        --mitochondria-mode \
        -O mix.filtered.vcf
    
    bcftools view --apply-filters PASS --types snps --max-alleles 2 --targets Pt --include "AF>0.01" -Oz mix.filtered.vcf > mix_Pt.vcf.gz
    gunzip mix_Pt.vcf.gz
    
    # 检测率分析
    echo "运行检测率分析..."
    mix_vcf="mix_Pt.vcf"
    # 使用1x的位点列表作为参考
    ler_unique="../Col100x_Ler1x/Ler_unique_sites.txt"
    col_sites="../Col100x_Ler1x/Col_sites.txt"
    
    ## 混合样本中检测到的Ler特有位点
    detected_all=$(grep -Ff $ler_unique $mix_vcf | wc -l)
    grep -Ff $ler_unique $mix_vcf | awk '!/#/ {print $1 "\t" $2}' | sort > detected_ler_sites.tmp
    detected_true=$(comm -23 detected_ler_sites.tmp <(sort $col_sites) | wc -l)
    total_ler=$(wc -l < $ler_unique)
    
    detection_rate_all=$(echo "scale=2; $detected_all / $total_ler * 100" | bc)
    detection_rate_true=$(echo "scale=2; $detected_true / $total_ler * 100" | bc)
    
    echo -e "Ler特有位点总数：$total_ler"
    echo -e "混合样本检测到的Ler位点（含Col假阳性）：$detected_all → 检测率：$detection_rate_all%"
    echo -e "排除Col假阳性后，真实Ler检测位点：$detected_true → 真实检测率：$detection_rate_true%"
    
    comm -23 detected_ler_sites.tmp <(sort $col_sites) > true_ler_sites_final.txt  # 真实Ler位点列表
    
    awk '
    BEGIN {
      # 读取真实Ler位点到字典（key: chrom\tpos，value: 1）
      while ((getline < "true_ler_sites_final.txt") > 0) {
        key = $1 "\t" $2
        ler_sites[key] = 1
      }
    }
    !/#/ {
      # 遍历VCF，匹配真实Ler位点
      key = $1 "\t" $2
      if (ler_sites[key] == 1) {
        print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10 "\t" $11  # 输出：chrom+pos+REF+ALT+Col样本+Ler样本
      }
    }
    ' $mix_vcf > true_ler_vcf_data.txt
    
    echo "从VCF中提取到的真实Ler位点数据行数：$(wc -l < true_ler_vcf_data.txt)"
    head -n 3 true_ler_vcf_data.txt
    
    
    echo -e "CHROM\tPOS\tREF\tLer_ALT\tCol_REF_reads\tLer_ALT_reads\t混合总reads\t异质性频率(%)" > ler_heteroplasmy_final.tsv
    while IFS=$'\t' read chrom pos ref alt col_sample ler_sample; do
      [[ -z "$chrom" || -z "$pos" ]] && continue
    
      col_ref_reads=$(echo "$col_sample" | cut -d':' -f2 | cut -d',' -f1 | awk '{print $1+0}')
      ler_alt_reads=$(echo "$ler_sample" | cut -d':' -f2 | cut -d',' -f2 | awk '{print $1+0}')
      total_reads=$((col_ref_reads + ler_alt_reads))
    
      # 计算异质性频率
      if [[ $total_reads -eq 0 || $ler_alt_reads -eq 0 ]]; then
        hetero_freq="0.0000"
      else
        hetero_freq=$(echo "scale=4; $ler_alt_reads / $total_reads * 100" | bc 2>/dev/null)
      fi
      
      echo -e "$chrom\t$pos\t$ref\t$alt\t$col_ref_reads\t$ler_alt_reads\t$total_reads\t$hetero_freq" >> ler_heteroplasmy_final.tsv
    
    done < true_ler_vcf_data.txt
    
    valid_sites=$(awk -F'\t' 'NR>1 && $8!="0.0000" {count++} END {print count+0}' ler_heteroplasmy_final.tsv)
    avg_hetero=$(awk -F'\t' 'NR>1 {sum+=$8} END {if(NR>1) print sum/(NR-1); else print 0}' ler_heteroplasmy_final.tsv | bc -l | xargs printf "%.4f")
    
    
    # 提取未检测到的Ler异质性位点
    echo "提取未检测到的Ler异质性位点..."
    
    # 生成未检测到的位点列表
    comm -23 <(sort $ler_unique) <(sort true_ler_sites_final.txt) > missed_ler_sites.txt
    
    echo -e "未检测到的Ler位点数量：$(wc -l < missed_ler_sites.txt)"
    
    # 从mix.filtered.vcf中提取这些未检测到位点的完整VCF信息
    awk '
    BEGIN {
      # 读取未检测到的Ler位点到字典
      while ((getline < "missed_ler_sites.txt") > 0) {
        key = $1 "\t" $2
        missed_sites[key] = 1
      }
    }
    /^#/ {
      # 输出VCF头信息
      print
      next
    }
    {
      # 检查当前位点是否在未检测到位点列表中
      key = $1 "\t" $2
      if (missed_sites[key] == 1) {
        print  # 输出完整的VCF行
      }
    }
    ' mix.filtered.vcf > missed_ler_sites_detailed.vcf
    
    echo "未检测到位点的详细VCF信息已保存到：missed_ler_sites_detailed.vcf"
    echo -e "提取到的未检测到位点VCF行数：$(grep -v '^#' missed_ler_sites_detailed.vcf | wc -l)"
    
    # 统计混合后新检测到的但原本Ler中没有的位点
    echo "统计混合后新检测到的位点..."
    
    # 从混合VCF中提取所有非头部的位点
    grep -v '^#' $mix_vcf | awk '{print $1 "\t" $2}' | sort > all_detected_sites.tmp
    
    # 生成新检测到的位点列表（在混合VCF中但不在原始Ler特有位点列表中）
    comm -23 <(sort all_detected_sites.tmp) <(sort $ler_unique) > newly_detected_sites.txt
    
    echo -e "混合后新检测到的位点数量：$(wc -l < newly_detected_sites.txt)"
    
    # 从mix.filtered.vcf中提取这些新检测到位点的完整VCF信息
    awk '
    BEGIN {
      # 读取新检测到的位点到字典
      while ((getline < "newly_detected_sites.txt") > 0) {
        key = $1 "\t" $2
        new_sites[key] = 1
      }
    }
    /^#/ {
      # 输出VCF头信息
      print
      next
    }
    {
      # 检查当前位点是否在新检测到位点列表中
      key = $1 "\t" $2
      if (new_sites[key] == 1) {
        print  # 输出完整的VCF行
      }
    }
    ' mix.filtered.vcf > newly_detected_sites_detailed.vcf
    
    echo "新检测到位点的详细VCF信息已保存到：newly_detected_sites_detailed.vcf"
    echo -e "提取到的新检测到位点VCF行数：$(grep -v '^#' newly_detected_sites_detailed.vcf | wc -l)"
    
    # 清理临时文件
    rm -f detected_ler_sites.tmp all_detected_sites.tmp
    
    echo "完成 Col100x_Ler${ratio}x 分析"
    cd ..
done

echo "所有比例分析完成！"
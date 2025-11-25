# preperation
```
mkdir -p mix2/genome
cp genome/col0/* mix2/genome
mkdir -p mix2/ena
cp ena/SRR61696* mix2/ena
mkdir mix2/SRR616966
cd mix2/SRR616966
mkdir 1_genome
cd 1_genome
cp ../../genome/genome.fa genome.fa
cp ../../genome/chr.sizes chr.sizes
cd ..
mkdir 2_illumina
cd 2_illumina
cp ../../ena/SRR616966_1.fastq.gz R1.fq.gz
cp ../../ena/SRR616966_2.fastq.gz R2.fq.gz
cd ../../ena
(head -1 ../../ena/ena_info.tsv && grep -E "(SRR616966|SRR616965)" ../../ena/ena_info.tsv) > ena_info.tsv
cd ..
export FOLD=2
export GENOME_SIZE=$(
    cat genome/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.tsv |
    tsv-select -H -f name,srr,bases |
    grep -E "(SRR616966|SRR616965)" |
    perl -nla -F'\t' -e '
        BEGIN { our %seen }
        /^name/ and next;
        $seen{$F[0]} and next;
        my $bases = $F[2];
        $bases =~ s/G$//;
        my $cutoff = $bases * 1000 * 1000 * 1000 / $ENV{GENOME_SIZE} * $ENV{FOLD};
        $cutoff = int $cutoff;
        print join qq(\t), ($F[0], $F[1], $cutoff, $ENV{FOLD});
        $seen{$F[0]}++;
    ' \
    > opts.tsv



#syn
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

        mkdir -p {2}_{4}/1_genome
        pushd {2}_{4}/1_genome

        cp ../../genome/genome.fa genome.fa
        cp ../../genome/chr.sizes chr.sizes
        popd > /dev/null

        mkdir -p {2}_{4}/2_illumina
        pushd {2}_{4}/2_illumina

        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd > /dev/null
    '
```

## run anchr 
```
cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        if [ ! -d {2}_{4} ]; then
            exit;
        fi

        if [ ! -e {2}_{4}/2_illumina/R1.fq.gz ]; then
            exit;
        fi
        if [ ! -e {2}_{4}/2_illumina/R2.fq.gz ]; then
            exit;
        fi

        if bjobs -w | tr -s " " | cut -d " " -f 7 | grep -w "^{2}$"; then
            echo Job {2}_{4} exists
            exit;
        fi

        cd {2}_{4}

        echo {2}_{4}

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

        bsub -q mpi -n 24 -J "{2}_{4}" "
            bash 0_script/0_master.sh
        "
    '
```
# vcf
```
cd ~/zxy/plastid/mix2/SRR616966_2/3_gatk
bcftools view --apply-filters PASS --max-alleles 2 --types snps --targets Pt --include "AF>0.01" -Oz R.filtered.vcf > Col_Pt.vcf.gz
gunzip Col_Pt.vcf.gz

cd ~/zxy/plastid/mix2/SRR616965_2/3_gatk
bcftools view --apply-filters PASS --max-alleles 2 --types snps --targets Pt --include "AF>0.01" -Oz R.filtered.vcf > Ler_Pt.vcf.gz
gunzip Ler_Pt.vcf.gz
# 提取Col 100×和Ler 1×的变异位点（仅保留位置）
cd ~/zxy/plastid/mix2/mix_results/Col100x_Ler1x
awk '!/#/ {print $1 "\t" $2}' ~/zxy/plastid/mix2/SRR616966_2/3_gatk/Col_Pt.vcf > Col_sites.txt
awk '!/#/ {print $1 "\t" $2}' ~/zxy/plastid/mix2/SRR616965_2/3_gatk/Ler_Pt.vcf > Ler_sites.txt

# 筛选Ler特有位点（Ler有而Col没有的位点）
comm -23 <(sort Ler_sites.txt) <(sort Col_sites.txt) > Ler_unique_sites.txt

echo "Col 100×的变异位点数量：$(wc -l < Col_sites.txt)"  # 0个
echo "Ler 1×的变异位点数量：$(wc -l < Ler_sites.txt)"    # 41个
echo "Ler特有的位点数量（Ler有、Col没有）：$(wc -l < Ler_unique_sites.txt)"

```
# mix
```
WORK=~/zxy/plastid/mix2/mix_results
mkdir -p "$WORK"
cd "$WORK"


COL_FOLD=100
LER_FOLD=1
SEED=123
WORK=~/zxy/plastid/mix2/mix_results
GENOME=~/zxy/plastid/mix2/genome/genome.fa
COL_BAM=~/zxy/plastid/mix2/SRR616966_2/3_bwa/R.sort.bam
LER_BAM=~/zxy/plastid/mix2/SRR616965_2/3_bwa/R.sort.bam
COL_VCF=~/zxy/plastid/mix2/SRR616966_2/3_gatk/R.filtered.vcf
LER_VCF=~/zxy/plastid/mix2/SRR616965_2/3_gatk/R.filtered.vcf


RESULT_DIR=$WORK/Col${COL_FOLD}x_Ler${LER_FOLD}x
mkdir -p $RESULT_DIR
cd $RESULT_DIR

cd ~/zxy/plastid/mix2/mix_results/Col100x_Ler1x
#col100x
samtools view -b -F 0x904 -q 1 $COL_BAM Pt > col_pt_raw.bam
samtools index col_pt_raw.bam
col_pt_raw_reads=$(samtools view -c col_pt_raw.bam)
echo "Col 原始 Pt 区域有效 reads：$col_pt_raw_reads"

#复制 100 次
for ((i=1; i<=COL_FOLD; i++)); do
  samtools addreplacerg -r "@RG\tID:Col_$i\tSM:Col_100x\tPL:ILLUMINA\tLB:Col\tPU:Col_$i" \
    col_pt_raw.bam -o col_pt_copy_$i.bam
done

samtools cat col_pt_copy_*.bam -o col_100x_merged.bam

samtools sort col_100x_merged.bam -o col_100x.sm.bam
samtools index col_100x.sm.bam

col_pt_raw_reads=$(samtools view -c col_pt_raw.bam)
col_100x_reads=$(samtools view -c col_100x.sm.bam)
echo "Col 原始 Pt 区域有效 reads：$col_pt_raw_reads"
echo "Col 100× Pt 区域有效 reads：$col_100x_reads"
echo "实际放大倍数：$((col_100x_reads / col_pt_raw_reads))"

#ler1x
samtools view -b -F 0x904 -q 1 $LER_BAM Pt -o ler_pt_temp.bam
samtools addreplacerg -r "@RG\tID:Ler\tSM:Ler_1x\tPL:ILLUMINA\tLB:Ler\tPU:Ler" \
  -o ler_1x_temp.bam \
  ler_pt_temp.bam

samtools sort ler_1x_temp.bam -o ler_1x.sm.bam
samtools index ler_1x.sm.bam
ler_1x_reads=$(samtools view -c ler_1x.sm.bam)
echo "Ler 1× Pt 区域有效 reads：$ler_1x_reads"

samtools merge -c mix_Col100x_Ler1x_temp.bam \
  col_100x.sm.bam \
  ler_1x.sm.bam
samtools sort mix_Col100x_Ler1x_temp.bam -o mix_Col100x_Ler1x.bam
samtools index mix_Col100x_Ler1x.bam
mix_reads=$(samtools view -c mix_Col100x_Ler1x.bam)
echo "混合样本 Pt 区域有效 reads：$mix_reads"
```
# gatk
```
cd ../genome
gatk CreateSequenceDictionary \
  -R genome.fa \
  -O genome.dict 
cd ~/zxy/plastid/mix2/mix_results/Col100x_Ler1x
gatk --java-options "-Xmx40g" \
    Mutect2 \
    --native-pair-hmm-threads 16 \
    -R ../../genome/genome.fa \
    -I mix_Col100x_Ler1x.bam \
    -L Pt \
    --mitochondria-mode \
    --max-reads-per-alignment-start 100 \
    --max-mnp-distance 0 \
    -O mix_Pt.raw.vcf

gatk --java-options "-Xmx80g" \
    FilterMutectCalls \
    -R ../../genome/genome.fa \
    -V mix_Pt.raw.vcf \
    --max-alt-allele-count 4 \
    --mitochondria-mode \
    -O mix.filtered.vcf

bcftools view --apply-filters PASS --types snps --max-alleles 2 --targets Pt --include "AF>0.01" -Oz mix.filtered.vcf > mix_Pt.vcf.gz
gunzip mix_Pt.vcf.gz

```
# detective rate
```
cd ~/zxy/plastid/mix2/mix_results/Col100x_Ler1x
mix_vcf="mix_Pt.vcf"
ler_unique="/share/home/wangq/zxy/plastid/mix2/mix_results/Col100x_Ler1x/Ler_unique_sites.txt"
col_sites="/share/home/wangq/zxy/plastid/mix2/mix_results/Col100x_Ler1x/Col_sites.txt"

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
```
# 遍历ler0.5x/2x/4x/8x/16x
```
bash mix.sh
```
# summary
```
WORK=~/zxy/plastid/mix2/mix_results  # 结果目录根路径
SUMMARY_FILE="heteroplasmy_summary.tsv"  # 汇总结果文件
LER_UNIQUE_SITES=~/zxy/plastid/mix2/mix_results/Col100x_Ler1x/Ler_unique_sites.txt  # Ler特有位点总数文件

echo -e "样本名称\tLer倍数\tLer特有位点总数\t真实检测位点\t真实检测率(%)\t有效异质位点\t平均异质性频率(%)" > "$SUMMARY_FILE"

total_ler_sites=$(wc -l < "$LER_UNIQUE_SITES")

for dir in "$WORK"/Col100x_Ler*; do
  [ ! -d "$dir" ] && continue
  dir_name=$(basename "$dir")
  echo -e "\n处理目录：$dir_name"

  # 提取Ler倍数
  ler_fold=$(echo "$dir_name" | sed 's/Col100x_Ler//; s/x//')

  heteroplasmy_tsv="$dir/ler_heteroplasmy_final.tsv"
  true_ler_sites="$dir/true_ler_sites_final.txt"

  true_detected_sites=$(wc -l < "$true_ler_sites")
  true_detection_rate=$(echo "scale=4; $true_detected_sites / $total_ler_sites * 100" | bc)
  valid_hetero_sites=$(awk -F"\t" 'NR>1 && $8!="0.0000" {count++} END {print count+0}' "$heteroplasmy_tsv")

  avg_hetero_freq=$(awk -F"\t" '
    NR>1 && $8!="0.0000" {
      sum+=$8; 
      count++
    } 
    END {
      if(count>0) print sum/count; 
      else print 0
    }
  ' "$heteroplasmy_tsv" | bc -l | xargs printf "%.4f")

  echo -e "$dir_name\t$ler_fold\t$total_ler_sites\t$true_detected_sites\t$true_detection_rate\t$valid_hetero_sites\t$avg_hetero_freq" >> "$SUMMARY_FILE"

  # 打印当前目录的统计结果（实时查看）
  echo "  Ler倍数：$ler_fold"
  echo "  真实检测位点：$true_detected_sites"
  echo "  真实检测率：$true_detection_rate%"
  echo "  有效异质位点：$valid_hetero_sites"
  echo "  平均异质性频率：$avg_hetero_freq%"
done
```
# plot
```
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(6)

work_dir <- "~/zxy/plastid/mix2/mix_results"
summary_file <- file.path(work_dir, "heteroplasmy_summary.tsv")
summary_df <- read.delim(summary_file, sep = "\t", stringsAsFactors = FALSE)

raw_data_list <- list()  # 存储所有样本的原始数据
sample_dirs <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)  # 遍历所有样本目录

for (dir in sample_dirs) {
  # 筛选Col100x_LerXx格式的目录
  if (grepl("Col100x_Ler", basename(dir))) {
    # 异质性原始数据文件路径
    hetero_file <- file.path(dir, "ler_heteroplasmy_final.tsv")
    if (file.exists(hetero_file)) {
      # 读取当前样本的原始数据
      df <- read.delim(hetero_file, sep = "\t", stringsAsFactors = FALSE)
      # 提取样本名称和Ler倍数（从目录名解析）
      sample_name <- basename(dir)
      ler_fold <- as.numeric(gsub("Col100x_Ler|x", "", sample_name))
      # 添加样本信息列
      df$样本名称 <- sample_name
      df$Ler倍数 <- ler_fold
      # 过滤无效位点（异质性频率=0的位点，可选保留）
      df <- df %>% filter(异质性频率... != 0)  # 只保留有效异质位点
      # 存入列表
      raw_data_list[[sample_name]] <- df
    }
  }
}

raw_combined_df <- bind_rows(raw_data_list)
ler_order <- sort(unique(raw_combined_df$Ler倍数))  # 提取有序的Ler倍数
raw_combined_df$Ler倍数 <- factor(raw_combined_df$Ler倍数, levels = ler_order, ordered = TRUE)  # 有序因子
summary_df$Ler倍数 <- factor(summary_df$Ler倍数, levels = ler_order, ordered = TRUE)
#boxplot
plot_box_dot <- ggplot(
  data = raw_combined_df,
  aes(x = Ler倍数, y = 异质性频率..., fill = Ler倍数, color = Ler倍数)
) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.6) +
  labs(x = "Ler Fold", y = "Heteroplasmy Frequency (%)") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid = element_blank(), 
    panel.border = element_blank(),
    axis.line = element_line(size = 0.8)
  ) +
  scale_fill_manual(values = npg_extended) +
  scale_color_manual(values = npg_extended)

ggsave("heteroplasmy_box_jitter.png", plot_box_dot, width = 6, height = 4, dpi = 300)

#lineplot
plot_line <- ggplot(
  data = summary_df,
  aes(x = Ler倍数, y = 有效异质位点, group = 1)
) +
  geom_line(linewidth = 2, color = npg_pal[2]) +
  geom_point(size = 4, shape = 21, fill = "white", color = npg_pal[2]) +
  geom_text(aes(label = 有效异质位点), vjust = -1, size = 3.5) +
  labs(x = "Ler Fold", y = "Number of Heteroplasmic Sites") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid = element_blank(),      # ❗ 去网格
    panel.border = element_blank(),    # ❗ 去边框
    axis.line = element_line(size = 0.8)
  )


ggsave("heteroplasmy_count_line.png", plot_line, width = 6, height = 4, dpi = 300)
```
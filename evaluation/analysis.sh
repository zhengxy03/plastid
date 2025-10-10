#!/bin/bash
# ================================
# 基因组残留reads过滤与注释分析（不含嵌合体部分）
# ================================

# ------------------------------
# 1️⃣ pretreatment
# ------------------------------
SAMPLE="SRR616966"
FOLDS=(1 2 4 8 16 32 64)
BASE_DIR="$HOME/zxy/plastid/evaluation"
RESULT_DIR="$BASE_DIR/summary_results"
mkdir -p $RESULT_DIR

# 初始化汇总表格（删除嵌合体相关列）
echo -e "fold\t总残留reads数\t转座元件reads数\t转座元件占比(%)\t多比对序列总数\t转座子多比对数\t转座子多比对占比(%)\t核平均read长度\t核chr1-read数\t核chr2-read数\t核chr3-read数\t核chr4-read数\t核chr5-read数\t线粒体平均覆盖度\t串联重复reads数\t串联重复占比(%)\t其他核残留reads数\t其他核残留占比(%)" > $RESULT_DIR/summary_stats.tsv

# 检查目录与bam文件
for fold in "${FOLDS[@]}"; do
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    BWA_DIR="$FOLD_DIR/3_bwa"
    if [ ! -d "$FOLD_DIR" ]; then
        echo "错误：未找到目录 $FOLD_DIR"
        exit 1
    fi
    if [ ! -f "$BWA_DIR/R.sort.bam" ]; then
        echo "错误：$BWA_DIR 中未找到 R.sort.bam 文件"
        exit 1
    fi
done

# 公共索引与注释准备
SHARED_GENOME_DIR="$BASE_DIR/${SAMPLE}_1/1_genome"
if [ ! -f "$SHARED_GENOME_DIR/genome.fa.fai" ]; then
    samtools faidx "$SHARED_GENOME_DIR/genome.fa"
fi

if [ ! -f "$SHARED_GENOME_DIR/nuclear_regions.bed" ]; then
    echo "生成区域BED文件..."
    cut -f1,2 "$SHARED_GENOME_DIR/genome.fa.fai" > "$SHARED_GENOME_DIR/genome.chrom.sizes"
    perl -ne 'chomp; @f=split(/\t/); print "$f[0]\t0\t$f[1]\n" if $f[0]=~/^[1-5]$/' "$SHARED_GENOME_DIR/genome.chrom.sizes" > "$SHARED_GENOME_DIR/nuclear_regions.bed"
    perl -ne 'chomp; @f=split(/\t/); print "$f[0]\t0\t$f[1]\n" if $f[0] eq "Mt"' "$SHARED_GENOME_DIR/genome.chrom.sizes" > "$SHARED_GENOME_DIR/mito_regions.bed"
    perl -ne 'chomp; @f=split(/\t/); print "$f[0]\t0\t$f[1]\n" if $f[0] eq "Pt"' "$SHARED_GENOME_DIR/genome.chrom.sizes" > "$SHARED_GENOME_DIR/chloroplast_regions.bed"
fi

if [ ! -f "$SHARED_GENOME_DIR/nuclear_te.bed" ]; then
    echo "生成转座元件注释..."
    perl -ne 'chomp;
              next if /^Transposon_Name/;
              @f = split(/\t/);
              $te_id = $f[0];
              if ($te_id =~ /^AT(\d)TE/) { $chr = $1; }
              else { next; }
              $start = $f[2] - 1;
              $end = $f[3];
              $family = $f[5];
              $super_family = $f[6];
              $type = "$super_family($family)";
              if ($chr =~ /^[1-5]$/) {
                  print "$chr\t$start\t$end\t$type\t$te_id\n";
              }' "$SHARED_GENOME_DIR/TAIR10_Transposable_Elements.txt" > "$SHARED_GENOME_DIR/nuclear_te.bed"
fi

if [ ! -f "$SHARED_GENOME_DIR/nuclear_tandem.bed" ]; then
    echo "警告：未找到串联重复注释文件，将使用空文件替代"
    touch "$SHARED_GENOME_DIR/nuclear_tandem.bed"
fi

# ------------------------------
# 2️⃣ 转座子分析
# ------------------------------
for fold in "${FOLDS[@]}"; do
    echo "----- 处理fold $fold（目录：${SAMPLE}_${fold}） -----"
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    WORK_DIR="$FOLD_DIR/3_bwa"
    mkdir -p $WORK_DIR
    cd $WORK_DIR

    samtools view -b -L "$SHARED_GENOME_DIR/nuclear_regions.bed" "$WORK_DIR/R.sort.bam" > nuclear_reads.bam
    samtools view -b -L "$SHARED_GENOME_DIR/mito_regions.bed" "$WORK_DIR/R.sort.bam" > mito_reads.bam
    samtools index nuclear_reads.bam
    samtools index mito_reads.bam

    bedtools bamtobed -i nuclear_reads.bam > nuclear_reads.bed
    bedtools bamtobed -i mito_reads.bam > mito_reads.bed

    # 核基因组分布
    perl -ne 'chomp; @f=split(/\t/); $count{$f[0]}++; END{
        foreach my $chr (1..5) { print "$chr\t".($count{$chr}//0)."\n"; }
    }' nuclear_reads.bed > $WORK_DIR/nuclear_chrom_dist.txt
    
    chr1=$(grep "^1" nuclear_chrom_dist.txt | cut -f2); chr1=${chr1:-0}
    chr2=$(grep "^2" nuclear_chrom_dist.txt | cut -f2); chr2=${chr2:-0}
    chr3=$(grep "^3" nuclear_chrom_dist.txt | cut -f2); chr3=${chr3:-0}
    chr4=$(grep "^4" nuclear_chrom_dist.txt | cut -f2); chr4=${chr4:-0}
    chr5=$(grep "^5" nuclear_chrom_dist.txt | cut -f2); chr5=${chr5:-0}

    # 平均长度
    avg_length=$(perl -ne 'chomp; @f=split(/\t/); $len=$f[2]-$f[1]; $sum+=$len; $count++; END{printf "%.2f", $sum/$count}' nuclear_reads.bed)

    # 线粒体覆盖度
    bedtools genomecov -i mito_reads.bed -g "$SHARED_GENOME_DIR/genome.chrom.sizes" -d | \
        perl -ne 'chomp; @f=split(/\t/); if($f[0] eq "Mt"){ $sum+=$f[2]; $count++ } END{printf "%.2f",$count>0?$sum/$count:0}' > mito_avg_coverage.txt
    mito_avg_cov=$(cat mito_avg_coverage.txt)

    # 转座子与多比对统计
    samtools view -f 256 nuclear_reads.bam | awk '{print $1}' | sort | uniq > nuclear_multimap_reads.txt
    bedtools intersect -abam nuclear_reads.bam -b "$SHARED_GENOME_DIR/nuclear_te.bed" -u > nuclear_te_reads.bam
    bedtools intersect -abam nuclear_reads.bam -b "$SHARED_GENOME_DIR/nuclear_te.bed" -v > nuclear_non_te_reads.bam
    bedtools bamtobed -i nuclear_reads.bam | grep -Ff nuclear_multimap_reads.txt > nuclear_multimap.bed
    bedtools intersect -a nuclear_multimap.bed -b "$SHARED_GENOME_DIR/nuclear_te.bed" -wa -wb | cut -f4 | sort | uniq > te_multimap_ids.txt
    
    total_reads=$(samtools view -c nuclear_reads.bam)
    te_reads=$(samtools view -c nuclear_te_reads.bam)
    te_ratio=$(echo "scale=2; $te_reads/$total_reads*100" | bc)
    multi_total=$(wc -l < nuclear_multimap_reads.txt)
    multi_te=$(wc -l < te_multimap_ids.txt)
    multi_te_ratio=$(echo "scale=2; $multi_te/$multi_total*100" | bc)

    echo -e "$fold\t$total_reads\t$te_reads\t$te_ratio\t$multi_total\t$multi_te\t$multi_te_ratio\t$avg_length\t$chr1\t$chr2\t$chr3\t$chr4\t$chr5\t$mito_avg_cov" > repeat_stats.tmp
done

# ------------------------------
# 3️⃣ 串联重复分析
# ------------------------------
for fold in "${FOLDS[@]}"; do
    echo "----- 处理fold $fold（目录：${SAMPLE}_${fold}） -----"
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    WORK_DIR="$FOLD_DIR/3_bwa"
    cd $WORK_DIR

    if [ ! -f "nuclear_reads.bam" ] || [ ! -f "nuclear_non_te_reads.bam" ]; then
        echo "错误：未找到重复序列分析的中间文件，请先运行转座子分析"
        exit 1
    fi

    bedtools intersect -abam nuclear_non_te_reads.bam -b "$SHARED_GENOME_DIR/nuclear_tandem.bed" -u > nuclear_tandem_reads.bam
    tandem_reads=$(samtools view -c nuclear_tandem_reads.bam)
    total_reads=$(samtools view -c nuclear_reads.bam)
    te_reads=$(samtools view -c nuclear_te_reads.bam)
    tandem_ratio=$(echo "scale=2; $tandem_reads/$total_reads*100" | bc)
    other_reads=$((total_reads - te_reads - tandem_reads))
    te_ratio=$(cut -f4 repeat_stats.tmp)
    other_ratio=$(echo "scale=2; 100 - $te_ratio - $tandem_ratio" | bc)
    echo -e "$tandem_reads\t$tandem_ratio\t$other_reads\t$other_ratio" > other_nuclear_stats.tmp
done

# ------------------------------
# 4️⃣ 汇总结果
# ------------------------------
for fold in "${FOLDS[@]}"; do
    WORK_DIR="$BASE_DIR/${SAMPLE}_${fold}/3_bwa"
    repeat_stats=$(cat $WORK_DIR/repeat_stats.tmp)
    other_nuclear_stats=$(cat $WORK_DIR/other_nuclear_stats.tmp)
    echo -e "$repeat_stats\t$other_nuclear_stats" >> $RESULT_DIR/summary_stats.tsv
    rm -f $WORK_DIR/repeat_stats.tmp $WORK_DIR/other_nuclear_stats.tmp
done

echo "✅ 所有fold分析完成，结果已保存至：$RESULT_DIR/summary_stats.tsv"

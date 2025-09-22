#!/bin/bash
# 定义样本和fold参数（与目录名完全对应）
SAMPLE="SRR616966"
FOLDS=(1 2 4 8 16 32 64)  # 与目录名中的fold值值完全匹配
BASE_DIR="$HOME/zxy/plastid/evaluation"
RESULT_DIR="$BASE_DIR/summary_results"
mkdir -p $RESULT_DIR

# 初始化汇总表格（调整统计项，移除高拷贝基因相关）
echo -e "fold\t总残留reads数\t转座元件reads数\t转座元件占比(%)\t多比对序列总数\t转座子多比对数\t转座子多比对占比(%)\t核平均read长度\t核chr1-read数\t核chr2-read数\t核chr3-read数\t核chr4-read数\t核chr5-read数\t线粒体平均覆盖度\t串联重复reads数\t串联重复占比(%)\t其他核残留reads数\t其他核残留占比(%)\t核-叶绿体嵌合体数\t核-线粒体嵌合体数\t叶绿体-线粒体嵌合体数\t总嵌合体数" > $RESULT_DIR/summary_stats.tsv


# ==============================================
# 预处理：检查所有fold目录并准备共享基因组数据
# ==============================================
echo "===== 验证目录结构并准备共享数据 ====="

# 检查所有fold目录及3_bwa子目录中的BAM文件
for fold in "${FOLDS[@]}"; do
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    BWA_DIR="$FOLD_DIR/3_bwa"  # BAM文件所在的子目录
    
    if [ ! -d "$FOLD_DIR" ]; then
        echo "错误：未找到目录 $FOLD_DIR，请检查目录命名是否正确"
        exit 1
    fi
    
    # 检查3_bwa子目录中是否有R.sort.bam
    if [ ! -f "$BWA_DIR/R.sort.bam" ]; then
        echo "错误：$BWA_DIR 中未找到 R.sort.bam 文件"
        exit 1
    fi
done

# 选择fold=1的基因组数据作为共享资源
SHARED_GENOME_DIR="$BASE_DIR/${SAMPLE}_1/1_genome"

# 确保基因组索引存在
if [ ! -f "$SHARED_GENOME_DIR/genome.fa.fai" ]; then
    echo "生成基因组索引..."
    samtools faidx "$SHARED_GENOME_DIR/genome.fa"
fi

# 确保区域BED文件存在
if [ ! -f "$SHARED_GENOME_DIR/nuclear_regions.bed" ]; then
    echo "生成区域BED文件..."
    cut -f1,2 "$SHARED_GENOME_DIR/genome.fa.fai" > "$SHARED_GENOME_DIR/genome.chrom.sizes"
    perl -ne 'chomp; @f=split(/\t/); print "$f[0]\t0\t$f[1]\n" if $f[0]=~/^[1-5]$/' "$SHARED_GENOME_DIR/genome.chrom.sizes" > "$SHARED_GENOME_DIR/nuclear_regions.bed"
    perl -ne 'chomp; @f=split(/\t/); print "$f[0]\t0\t$f[1]\n" if $f[0] eq "Mt"' "$SHARED_GENOME_DIR/genome.chrom.sizes" > "$SHARED_GENOME_DIR/mito_regions.bed"
    perl -ne 'chomp; @f=split(/\t/); print "$f[0]\t0\t$f[1]\n" if $f[0] eq "Pt"' "$SHARED_GENOME_DIR/genome.chrom.sizes" > "$SHARED_GENOME_DIR/chloroplast_regions.bed"
fi

# 确保转座元件和串联重复注释存在
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

# 添加串联重复区域注释
if [ ! -f "$SHARED_GENOME_DIR/nuclear_tandem.bed" ]; then
    echo "警告：未找到串联重复注释文件，将使用空文件替代"
    touch "$SHARED_GENOME_DIR/nuclear_tandem.bed"
    # 实际使用时建议替换为：
    # trf "$SHARED_GENOME_DIR/genome.fa" 2 5 7 80 10 50 500 -h -ngs > "$SHARED_GENOME_DIR/tandem_repeats.txt"
    # perl trf.pl "$SHARED_GENOME_DIR/tandem_repeats.txt" > "$SHARED_GENOME_DIR/nuclear_tandem.bed"
fi

if [ ! -f "$SHARED_GENOME_DIR/combined_ref_db.nhr" ]; then
    echo "构建BLAST数据库..."
    makeblastdb -in "$SHARED_GENOME_DIR/genome.fa" -dbtype nucl -out "$SHARED_GENOME_DIR/combined_ref_db" -parse_seqids
fi


# ==============================================
# 第一个循环：处理所有fold的重复序列分析
# ==============================================
echo "===== 开始所有fold的重复序列分析 ====="
for fold in "${FOLDS[@]}"; do
    echo "----- 处理fold $fold（目录：${SAMPLE}_${fold}） -----"
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    WORK_DIR="$FOLD_DIR/3_bwa"  # 工作目录和BAM文件目录一致
    mkdir -p $WORK_DIR
    cd $WORK_DIR

    # 从3_bwa子目录提取残留reads（BAM文件就在当前目录）
    samtools view -b -L "$SHARED_GENOME_DIR/nuclear_regions.bed" "$WORK_DIR/R.sort.bam" > nuclear_reads.bam
    samtools view -b -L "$SHARED_GENOME_DIR/mito_regions.bed" "$WORK_DIR/R.sort.bam" > mito_reads.bam
    samtools index nuclear_reads.bam
    samtools index mito_reads.bam

    # 转成bed文件用于覆盖分析
    bedtools bamtobed -i nuclear_reads.bam > nuclear_reads.bed
    bedtools bamtobed -i mito_reads.bam > mito_reads.bed

    # 统计残留reads的覆盖区域特征
    # 1. 核基因组残留reads的染色体分布
    perl -ne 'chomp; @f=split(/\t/); $count{$f[0]}++; END{
        foreach my $chr (1..5) {
            print "$chr\t".($count{$chr}//0)."\n";
        }
    }' nuclear_reads.bed > $WORK_DIR/nuclear_chrom_dist.txt
    
    chr1=$(grep "^1" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr2=$(grep "^2" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr3=$(grep "^3" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr4=$(grep "^4" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr5=$(grep "^5" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr1=${chr1:-0}; chr2=${chr2:-0}; chr3=${chr3:-0}; chr4=${chr4:-0}; chr5=${chr5:-0}
    
    # 2. 核基因组残留reads的平均长度
    avg_length=$(perl -ne 'chomp; @f=split(/\t/); $len=$f[2]-$f[1]; $sum+=$len; $count++; END{printf "%.2f", $sum/$count}' nuclear_reads.bed)
    
    # 3. 线粒体残留reads的平均覆盖度
    bedtools genomecov -i mito_reads.bed -g "$SHARED_GENOME_DIR/genome.chrom.sizes" -d | \
        perl -ne 'chomp; @f=split(/\t/); 
                 if($f[0] eq "Mt") { $sum += $f[2]; $count++ }
                 END{ printf "%.2f", $count>0 ? $sum/$count : 0 }' > $WORK_DIR/mito_avg_coverage.txt
    mito_avg_cov=$(cat $WORK_DIR/mito_avg_coverage.txt)

    # 重复序列核心分析
    samtools view -f 256 nuclear_reads.bam | perl -ne 'chomp; @f=split(/\t/); print "$f[0]\n"' | sort | uniq > nuclear_multimap_reads.txt
    
    bedtools intersect -abam nuclear_reads.bam -b "$SHARED_GENOME_DIR/nuclear_te.bed" -u > nuclear_te_reads.bam
    bedtools intersect -abam nuclear_reads.bam -b "$SHARED_GENOME_DIR/nuclear_te.bed" -v > nuclear_non_te_reads.bam
    
    bedtools bamtobed -i nuclear_reads.bam | grep -Ff nuclear_multimap_reads.txt > nuclear_multimap.bed
    bedtools intersect -a nuclear_multimap.bed -b "$SHARED_GENOME_DIR/nuclear_te.bed" -wa -wb | cut -f4 | sort | uniq > te_multimap_ids.txt
    
    # 计算统计值（确保保留两位小数）
    total_reads=$(samtools view -c nuclear_reads.bam)
    te_reads=$(samtools view -c nuclear_te_reads.bam)
    te_ratio=$(echo "scale=4; $te_reads/$total_reads*100" | bc | xargs printf "%.2f")
    multi_total=$(wc -l < nuclear_multimap_reads.txt)
    multi_te=$(wc -l < te_multimap_ids.txt)
    multi_te_ratio=$(echo "scale=4; $multi_te/$multi_total*100" | bc | xargs printf "%.2f")
    
    # 保存到临时文件
    echo -e "$fold\t$total_reads\t$te_reads\t$te_ratio\t$multi_total\t$multi_te\t$multi_te_ratio\t$avg_length\t$chr1\t$chr2\t$chr3\t$chr4\t$chr5\t$mito_avg_cov" > $WORK_DIR/repeat_stats.tmp
done


# ==============================================
# 第二个循环：处理所有fold的其他核残留reads分析
# ==============================================
echo "===== 开始所有fold的其他核残留reads分析 ====="
for fold in "${FOLDS[@]}"; do
    echo "----- 处理fold $fold（目录：${SAMPLE}_${fold}） -----"
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    WORK_DIR="$FOLD_DIR/3_bwa"
    mkdir -p $WORK_DIR
    cd $WORK_DIR

    # 检查必要的中间文件
    if [ ! -f "nuclear_reads.bam" ] || [ ! -f "nuclear_non_te_reads.bam" ]; then
        echo "错误：未找到重复序列分析的中间文件，请先运行第一个循环"
        exit 1
    fi

    # 1. 串联重复区域的reads数（排除已归类为转座子的reads）
    bedtools intersect -abam nuclear_non_te_reads.bam -b "$SHARED_GENOME_DIR/nuclear_tandem.bed" -u > nuclear_tandem_reads.bam
    tandem_reads=$(samtools view -c nuclear_tandem_reads.bam)
    total_reads=$(samtools view -c nuclear_reads.bam)
    te_reads=$(samtools view -c nuclear_te_reads.bam)
    
    # 计算占比（确保保留两位小数）
    tandem_ratio=$(echo "scale=4; $tandem_reads/$total_reads*100" | bc | xargs printf "%.2f")
    
    # 2. 其他核残留reads数（未被上述类别包含）
    other_reads=$((total_reads - te_reads - tandem_reads))
    
    # 其他核残留占比 = 100% - 转座子占比 - 串联重复占比（确保保留两位小数）
    te_ratio=$(grep -w "$fold" $WORK_DIR/repeat_stats.tmp | cut -f4)
    other_ratio=$(echo "scale=4; 100 - $te_ratio - $tandem_ratio" | bc | xargs printf "%.2f")
    
    # 保存到临时文件
    echo -e "$tandem_reads\t$tandem_ratio\t$other_reads\t$other_ratio" > $WORK_DIR/other_nuclear_stats.tmp
done


# ==============================================
# 第三个循环：处理所有fold的嵌合体序列分析
# ==============================================
echo "===== 开始所有fold的嵌合体序列分析 ====="
for fold in "${FOLDS[@]}"; do
    echo "----- 处理fold $fold（目录：${SAMPLE}_${fold}） -----"
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    WORK_DIR="$FOLD_DIR/3_bwa"
    mkdir -p $WORK_DIR
    cd $WORK_DIR

    # 检查必要的中间文件
    if [ ! -f "nuclear_multimap_reads.txt" ] || [ ! -f "te_multimap_ids.txt" ]; then
        echo "错误：未找到重复序列分析的中间文件，请先运行第一个循环"
        exit 1
    fi

    # 筛选非转座元件的多比对ID
    comm -23 <(sort nuclear_multimap_reads.txt | uniq | sed 's/[[:space:]]//g') \
             <(sort te_multimap_ids.txt | uniq | sed 's/[[:space:]]//g') \
             > non_te_multimap_ids_clean.txt
    
    # 提取潜在嵌合体序列
    samtools view -h nuclear_reads.bam > full_reads.sam
    perl -ne '
        BEGIN {
            open(my $fh, "<", "non_te_multimap_ids_clean.txt") or die $!;
            while (my $line = <$fh>) {
                chomp $line;
                $ids{$line} = 1 if $line ne "";
            }
            close $fh;
        }
        if ($. == 1 || /^@/) { print; next; }
        my @fields = split(/\t/);
        print if exists $ids{$fields[0]};
    ' full_reads.sam > target_reads.sam
    samtools fasta target_reads.sam > non_te_multimap.fasta
    rm -f full_reads.sam target_reads.sam
    
    # BLAST筛选真实嵌合体
    blastn -query non_te_multimap.fasta \
           -db "$SHARED_GENOME_DIR/combined_ref_db" \
           -evalue 1e-20 \
           -outfmt "6 qseqid sseqid pident length evalue" \
           -max_target_seqs 3 \
           -num_threads 4 > blast_temp.txt
    
    # 筛选嵌合体ID
    perl -ne 'chomp; @f=split(/\t/); next if @f<5;
             $id=$f[0]; $chr=$f[1]; $p=$f[2]; $e=$f[4];
             $g = ($chr=~/^[1-5]$/) ? "nuc" :
                  ($chr eq "Pt") ? "cp" :
                  ($chr eq "Mt") ? "mt" : "other";
             next if $g eq "other";
             if($p>=70 && $e<=1e-10) { push @{$h{$id}}, $g; }
             END{
                 for my $id (keys %h) {
                     my %unique = map { $_ => 1 } @{$h{$id}};
                     my @types = keys %unique;
                     if(@types >= 2) {
                         if( ($types[0] eq "nuc" && $types[1] eq "cp") || ($types[0] eq "cp" && $types[1] eq "nuc") ) {
                             print "$id\n";
                         } elsif( ($types[0] eq "nuc" && $types[1] eq "mt") || ($types[0] eq "mt" && $types[1] eq "nuc") ) {
                             print "$id\n";
                         } elsif( ($types[0] eq "cp" && $types[1] eq "mt") || ($types[0] eq "mt" && $types[1] eq "cp") ) {
                             print "$id\n";
                         }
                     }
                 }
             }' blast_temp.txt > candidate_chimeric_ids.txt
    
    # 提取并验证嵌合体序列
    sed 's/\/[12]$//' candidate_chimeric_ids.txt > candidate_ids_fixed.txt
    samtools view -h nuclear_reads.bam > full_reads2.sam
    perl -ne '
        BEGIN {
            open(my $fh, "<", "candidate_ids_fixed.txt") or die $!;
            while (my $line = <$fh>) {
                chomp $line;
                $ids{$line} = 1 if $line ne "";
            }
            close $fh;
        }
        if ($. == 1 || /^@/) { print; next; }
        my @fields = split(/\t/);
        print if exists $ids{$fields[0]};
    ' full_reads2.sam > target_reads2.sam
    samtools fasta target_reads2.sam > true_chimeric_seqs.fasta
    rm -f full_reads2.sam target_reads2.sam candidate_ids_fixed.txt
    
    # 统计嵌合体类型
    blastn -query true_chimeric_seqs.fasta \
           -db "$SHARED_GENOME_DIR/combined_ref_db" \
           -evalue 1e-10 \
           -outfmt "6 qseqid sseqid length pident evalue" \
           -max_target_seqs 3 \
           -num_threads 4 > chimeric_blast_verify.txt
    
    perl -ne 'chomp; @f=split(/\t/); next if @f<5;
             $id=$f[0]; $chr=$f[1]; $len=$f[2]; $p=$f[3];
             $g = ($chr=~/^[1-5]$/) ? "nuc" :
                  ($chr eq "Pt") ? "cp" :
                  ($chr eq "Mt") ? "mt" : "other";
             next if $g eq "other";
             $len_sum{$id}{$g} += $len;
             $p_sum{$id} += $p;
             $p_cnt{$id}++;
             END{
                 my ($nuc_cp, $nuc_mt, $cp_mt) = (0,0,0);
                 for my $id (keys %len_sum) {
                     my $nuc = $len_sum{$id}{nuc} // 0;
                     my $cp = $len_sum{$id}{cp} // 0;
                     my $mt = $len_sum{$id}{mt} // 0;
                     next if $nuc+$cp+$mt == 0;
                     
                     if ($nuc>0 && $cp>0) {
                         $nuc_cp++;
                     } elsif ($nuc>0 && $mt>0) {
                         $nuc_mt++;
                     } elsif ($cp>0 && $mt>0) {
                         $cp_mt++;
                     }
                 }
                 my $total = $nuc_cp + $nuc_mt + $cp_mt;
                 print "$nuc_cp\t$nuc_mt\t$cp_mt\t$total";
             }' chimeric_blast_verify.txt > $WORK_DIR/chimera_stats.tmp
done


# ==============================================
# 汇总结果
# ==============================================
echo "===== 汇总所有分析结果 ====="
for fold in "${FOLDS[@]}"; do
    WORK_DIR="$BASE_DIR/${SAMPLE}_${fold}/3_bwa"
    repeat_stats=$(cat $WORK_DIR/repeat_stats.tmp)
    other_nuclear_stats=$(cat $WORK_DIR/other_nuclear_stats.tmp)
    chimera_stats=$(cat $WORK_DIR/chimera_stats.tmp)
    echo -e "$repeat_stats\t$other_nuclear_stats\t$chimera_stats" >> $RESULT_DIR/summary_stats.tsv
    rm -f $WORK_DIR/repeat_stats.tmp $WORK_DIR/other_nuclear_stats.tmp $WORK_DIR/chimera_stats.tmp
done

echo "所有分析完成！汇总结果位于：$RESULT_DIR/summary_stats.tsv"


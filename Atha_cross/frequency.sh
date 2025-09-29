#!/bin/bash
set -euo pipefail

# 项目根目录
PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross"

# 叶绿体染色体名称
CHLORO_NAME="Pt"

# 验证项目根目录
if [ ! -d "$PROJECT_ROOT" ]; then
    echo "错误: 项目根目录不存在 - $PROJECT_ROOT"
    exit 1
fi

# 验证opts.tsv
opts_file="$PROJECT_ROOT/opts.tsv"
if [ ! -f "$opts_file" ]; then
    echo "错误: 未找到opts.tsv - $opts_file"
    exit 1
fi

# 提取样本名（去除前后空格和空行）
samples=$(awk -F'\t' 'NR>0 && $1!="" {gsub(/^[ \t\r]+|[ \t\r]+$/, ""); print $1}' "$opts_file")
if [ -z "$samples" ]; then
    echo "错误: 未提取到样本名"
    exit 1
fi

# 循环处理样本
for sample in $samples; do
    echo -e "\n===== 处理样本: $sample ====="
    
    # 定义文件路径（统一使用变量，避免硬编码错误）
    bam_file="$PROJECT_ROOT/$sample/3_bwa/R.sort.bam"
    bam_index="$PROJECT_ROOT/$sample/3_bwa/R.sort.bai"
    gatk_vcf="$PROJECT_ROOT/$sample/3_gatk/R.filtered.vcf"
    gatk_vcf_gz="$PROJECT_ROOT/$sample/3_gatk/R.filtered.vcf.gz"
    result_file="$PROJECT_ROOT/$sample/chloroplast_heteroplasmy_freq.txt"
    chloro_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_from_gatk.vcf.gz"
    chloro_vcf_tbi="$chloro_vcf.tbi"
    temp_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_temp.vcf"
    filtered_temp_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_filtered_temp.vcf"  # 新增：筛选后的临时VCF
    log_file="$PROJECT_ROOT/$sample/chloro_extraction.log"
    
    # 检查必要文件（BAM、索引、GATK原始VCF）
    if [ ! -f "$bam_file" ] || [ ! -f "$bam_index" ] || [ ! -f "$gatk_vcf" ]; then
        echo "警告: 缺少必要文件（BAM/索引/GATK VCF），跳过样本 $sample"
        continue
    fi

    # 检查任务是否已在运行（避免重复提交）
    if bjobs -w | grep -q "${sample}_chloro"; then
        echo "任务已存在（${sample}_chloro），跳过: $sample"
        continue
    fi
    
    # 创建临时目录和空日志
    tmp_dir="$PROJECT_ROOT/$sample/tmp_scripts"
    mkdir -p "$tmp_dir"
    > "$log_file"  # 清空旧日志，避免累积
    
    # 1. 覆盖度计算脚本（不变）
    cat > "$tmp_dir/calc_coverage.pl" << 'EOF'
#!/usr/bin/perl
use strict;
use warnings;
my ($sum, $count) = (0, 0);
while (<STDIN>) {
    my @f = split;
    next unless @f >= 3 && $f[2] =~ /^\d+$/;
    $sum += $f[2];
    $count++;
}
printf "%.2f\n", $count > 0 ? $sum/$count : 0;
EOF
    chmod +x "$tmp_dir/calc_coverage.pl"
    
    # 2. 频率计算脚本（不变）
    cat > "$tmp_dir/calc_frequency.pl" << 'EOF'
#!/usr/bin/perl
use strict;
use warnings;
while (<STDIN>) {
    chomp;
    my @f = split(/\t/);
    next unless @f >= 6;
    
    my ($chrom, $pos, $ref, $alt, $ad_str, $qual) = @f[0..5];
    
    $qual = 0 unless (defined $qual && $qual =~ /^[-+]?\d+(\.\d+)?$/);
    
    my @ad = split(/,/, $ad_str // "0,0");
    my $ref_cnt = (defined $ad[0] && $ad[0] =~ /^\d+$/) ? $ad[0] : 0;
    my $alt_cnt = (defined $ad[1] && $ad[1] =~ /^\d+$/) ? $ad[1] : 0;
    my $total = $ref_cnt + $alt_cnt;
    
    if ($alt_cnt >= 1) {
        my $freq = $total > 0 ? $alt_cnt / $total : 0;
        printf "%s\t%s\t%s\t%s\t%.6f\t%d\n", 
            $chrom, $pos, $ref, $alt, $freq, $total;
    }
}
EOF
    chmod +x "$tmp_dir/calc_frequency.pl"
    
    # 3. 提交LSF任务（核心修改：加入PASS和AF>0.01筛选）
    bsub -q mpi -n 24 -J "${sample}_chloro" <<EOT
        echo '开始处理样本 $sample' | tee -a '$log_file'
        
        # 步骤1：从BAM提取叶绿体区域reads并过滤低质量
        echo '处理BAM文件：提取叶绿体区域并过滤低质量reads...' | tee -a '$log_file'
        # 统计叶绿体区域总reads数（用于日志记录）
        chloro_total_reads=\$(samtools view -c '$bam_file' '$CHLORO_NAME')
        echo "叶绿体区域原始reads数: \$chloro_total_reads" | tee -a '$log_file'
        # 提取叶绿体区域，过滤MAPQ<10的低质量比对reads
        samtools view -b -q 10 -o '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam' '$bam_file' '$CHLORO_NAME'
        
        # 检查临时BAM是否有效（非空）
        if [ ! -s '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam' ]; then
            echo '警告: BAM中未找到有效叶绿体区域（空文件）' | tee -a '$log_file'
            touch '$result_file'  # 创建空结果文件，避免下游报错
            rm -f '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam'
            exit 0
        fi
        
        # 排序并索引叶绿体BAM
        samtools sort -o '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam' '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam'
        samtools index '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam'
        
        # 步骤2：计算叶绿体区域平均覆盖度
        echo '计算叶绿体区域平均覆盖度...' | tee -a '$log_file'
        avg_depth=\$(samtools depth '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam' | '$tmp_dir/calc_coverage.pl')
        echo "叶绿体区域平均覆盖度: \$avg_depth×" | tee -a '$log_file'
        
        # 步骤3：压缩GATK原始VCF并建立索引
        echo '压缩GATK原始VCF并生成索引...' | tee -a '$log_file'
        # 清理旧的压缩文件（避免残留）
        if [ -f '$gatk_vcf_gz' ] || [ -f '$gatk_vcf_gz.tbi' ]; then
            rm -f '$gatk_vcf_gz' '$gatk_vcf_gz.tbi'
        fi
        bgzip -c '$gatk_vcf' > '$gatk_vcf_gz'
        tabix -p vcf '$gatk_vcf_gz'
        
        # 验证VCF压缩是否成功
        if ! file '$gatk_vcf_gz' | grep -q 'gzip compressed data'; then
            echo '错误: GATK VCF压缩失败（非gzip格式）' | tee -a '$log_file'
            exit 1
        fi
        
        # 步骤4：提取叶绿体变异并筛选（核心修改：加入PASS和AF>0.01）
        echo '提取叶绿体变异并筛选（PASS标签 + AF>0.01）...' | tee -a '$log_file'
        # 第一步：提取叶绿体区域（--targets Pt）并保留PASS过滤标签的变异，压缩输出
        bcftools view \
            --apply-filters PASS \
            --targets '$CHLORO_NAME' \
            -Oz '$gatk_vcf_gz' > '$filtered_temp_vcf.gz'
        # 建立筛选后临时VCF的索引
        tabix -p vcf '$filtered_temp_vcf.gz'
        
        # 第二步：进一步筛选AF>0.01的变异，输出未压缩临时VCF（用于统计变异数）
        bcftools view \
            --include "AF>0.01" \
            '$filtered_temp_vcf.gz' > '$temp_vcf'
        
        # 统计筛选后的变异数（排除VCF头部注释行）
        var_count=\$(grep -v "^#" '$temp_vcf' | wc -l)
        echo "筛选后提取到的叶绿体变异数（PASS + AF>0.01）: \$var_count" | tee -a '$log_file'
        
        # 若无变异，创建空结果文件并退出
        if [ \$var_count -eq 0 ]; then
            echo '警告: 未提取到符合条件的叶绿体变异（PASS + AF>0.01）' | tee -a '$log_file'
            touch '$result_file'
            rm -f '$temp_vcf' '$filtered_temp_vcf.gz' '$filtered_temp_vcf.gz.tbi'  # 清理临时文件
            exit 0
        fi
        
        # 步骤5：清理旧的目标VCF并生成最终压缩VCF
        echo '清理旧文件并生成最终叶绿体VCF...' | tee -a '$log_file'
        if [ -f '$chloro_vcf' ] || [ -f '$chloro_vcf_tbi' ]; then
            echo '清理旧的叶绿体VCF文件和索引...' | tee -a '$log_file'
            rm -f '$chloro_vcf' '$chloro_vcf_tbi'
        fi
        # 压缩筛选后的临时VCF为最终目标文件
        bgzip -c '$temp_vcf' > '$chloro_vcf'
        tabix -p vcf '$chloro_vcf'
        
        # 验证最终VCF和索引是否有效
        if [ ! -f '$chloro_vcf' ] || [ ! -s '$chloro_vcf' ]; then
            echo '错误: 最终叶绿体VCF生成失败（文件不存在或为空）' | tee -a '$log_file'
            exit 1
        fi
        if [ ! -f '$chloro_vcf_tbi' ]; then
            echo '错误: 最终叶绿体VCF索引生成失败' | tee -a '$log_file'
            exit 1
        fi
        echo '最终叶绿体VCF（chloroplast_variants_from_gatk.vcf.gz）生成成功' | tee -a '$log_file'
        
        # 步骤6：计算异质性频率（基于筛选后的VCF）
        echo '计算叶绿体变异的异质性频率...' | tee -a '$log_file'
        # 用bcftools query提取关键字段，管道给Perl脚本计算频率
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t%QUAL\n' '$chloro_vcf' | \
        '$tmp_dir/calc_frequency.pl' > '$result_file'
        
        # 步骤7：记录最终结果统计
        result_count=\$(wc -l '$result_file' | cut -d' ' -f1)
        echo "最终保留的异质性变异数（alt_cnt≥1）: \$result_count" | tee -a '$log_file'
        
        # 步骤8：清理所有临时文件（避免磁盘占用）
        echo '清理临时文件...' | tee -a '$log_file'
        rm -f '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam'
        rm -f '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam'*  # 包含BAM和索引
        rm -f '$gatk_vcf_gz' '$gatk_vcf_gz.tbi'  # 清理GATK压缩VCF
        rm -f '$temp_vcf' '$filtered_temp_vcf.gz' '$filtered_temp_vcf.gz.tbi'  # 清理筛选临时文件
        rm -rf '$tmp_dir'  # 清理临时脚本目录
        
        echo '样本 $sample 处理完成' | tee -a '$log_file'
EOT
    # 注意：EOT必须顶格，且前面无任何字符（包括空格）
    
    echo "样本 $sample 任务已提交（任务名：${sample}_chloro）"
done

echo -e "\n所有样本任务提交完毕！可通过 bjobs 查看任务状态。"
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

# 提取样本名
samples=$(awk -F'\t' 'NR>0 && $1!="" {gsub(/^[ \t\r]+|[ \t\r]+$/, ""); print $1}' "$opts_file")
if [ -z "$samples" ]; then
    echo "错误: 未提取到样本名"
    exit 1
fi

# 循环处理样本
for sample in $samples; do
    echo -e "\n===== 处理样本: $sample ====="
    
    # 定义文件路径
    bam_file="$PROJECT_ROOT/$sample/3_bwa/R.sort.bam"
    bam_index="$PROJECT_ROOT/$sample/3_bwa/R.sort.bai"
    gatk_vcf="$PROJECT_ROOT/$sample/3_gatk/R.raw.vcf"
    gatk_vcf_gz="$PROJECT_ROOT/$sample/3_gatk/R.raw.vcf.gz"
    result_file="$PROJECT_ROOT/$sample/chloroplast_heteroplasmy_freq.txt"
    chloro_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_from_gatk.vcf.gz"
    chloro_vcf_tbi="$chloro_vcf.tbi"
    temp_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_temp.vcf"
    log_file="$PROJECT_ROOT/$sample/chloro_extraction.log"
    
    # 检查必要文件
    if [ ! -f "$bam_file" ] || [ ! -f "$bam_index" ] || [ ! -f "$gatk_vcf" ]; then
        echo "警告: 缺少必要文件，跳过样本 $sample"
        continue
    fi

    # 检查任务是否存在
    if bjobs -w | grep -q "${sample}_chloro"; then
        echo "任务已存在，跳过: $sample"
        continue
    fi
    
    # 创建临时目录和日志
    tmp_dir="$PROJECT_ROOT/$sample/tmp_scripts"
    mkdir -p "$tmp_dir"
    > "$log_file"
    
    # 覆盖度计算脚本
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
    
    # 频率计算脚本
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
    
    # 提交任务（核心修复：修正变量解析和引号使用）
    bsub -q mpi -n 24 -J "${sample}_chloro" <<EOT
        echo '开始处理样本 $sample' | tee -a '$log_file'
        
        # 1. 处理BAM文件
        echo '处理BAM文件...' | tee -a '$log_file'
        samtools view -c '$bam_file' '$CHLORO_NAME' | tee -a '$log_file'
        samtools view -b -q 10 -o '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam' '$bam_file' '$CHLORO_NAME'
        
        if [ ! -s '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam' ]; then
            echo '警告: BAM中未找到有效叶绿体区域' | tee -a '$log_file'
            touch '$result_file'
            rm -f '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam'
            exit 0
        fi
        
        samtools sort -o '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam' '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam'
        samtools index '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam'
        
        # 2. 计算覆盖度
        avg_depth=\$(samtools depth '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam' | '$tmp_dir/calc_coverage.pl')
        echo "平均覆盖度: \$avg_depth×" | tee -a '$log_file'
        
        # 3. 压缩GATK VCF
        echo '压缩GATK VCF文件...' | tee -a '$log_file'
        if [ -f '$gatk_vcf_gz' ] || [ -f '$gatk_vcf_gz.tbi' ]; then
            rm -f '$gatk_vcf_gz' '$gatk_vcf_gz.tbi'
        fi
        bgzip -c '$gatk_vcf' > '$gatk_vcf_gz'
        tabix -p vcf '$gatk_vcf_gz'
        
        if ! file '$gatk_vcf_gz' | grep -q 'gzip compressed data'; then
            echo '错误: GATK VCF压缩失败' | tee -a '$log_file'
            exit 1
        fi
        
        # 4. 提取叶绿体变异（核心修复：正确计算var_count）
        echo '统计并提取叶绿体变异...' | tee -a '$log_file'
        bcftools view -r '$CHLORO_NAME' '$gatk_vcf_gz' > '$temp_vcf'
        
        # 修复变量解析：使用正确的语法计算变异数
        var_count=\$(grep -v "^#" '$temp_vcf' | wc -l)
        echo "提取到的未压缩变异数: \$var_count" | tee -a '$log_file'
        
        if [ \$var_count -eq 0 ]; then
            echo '警告: 未提取到叶绿体变异' | tee -a '$log_file'
            touch '$result_file'
            rm -f '$temp_vcf'
            exit 0
        fi
        
        # 清理旧文件
        if [ -f '$chloro_vcf' ] || [ -f '$chloro_vcf_tbi' ]; then
            echo '清理旧的叶绿体VCF文件...' | tee -a '$log_file'
            rm -f '$chloro_vcf' '$chloro_vcf_tbi'
        fi
        
        # 压缩生成目标文件
        echo '压缩生成叶绿体VCF...' | tee -a '$log_file'
        bgzip -c '$temp_vcf' > '$chloro_vcf'
        tabix -p vcf '$chloro_vcf'
        
        # 验证目标文件
        if [ ! -f '$chloro_vcf' ] || [ ! -s '$chloro_vcf' ]; then
            echo '错误: chloroplast_variants_from_gatk.vcf.gz 生成失败' | tee -a '$log_file'
            exit 1
        fi
        if [ ! -f '$chloro_vcf_tbi' ]; then
            echo '错误: chloroplast_variants_from_gatk.vcf.gz.tbi 索引生成失败' | tee -a '$log_file'
            exit 1
        fi
        echo 'chloroplast_variants_from_gatk.vcf.gz 生成成功' | tee -a '$log_file'
        
        # 5. 计算异质性频率（修复格式符错误）
        echo '计算异质性频率...' | tee -a '$log_file'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t%QUAL\n' '$chloro_vcf' | \
        '$tmp_dir/calc_frequency.pl' > '$result_file'
        
        # 6. 记录结果
        result_count=\$(wc -l '$result_file' | cut -d' ' -f1)
        echo "最终保留的异质性变异数: \$result_count" | tee -a '$log_file'
        
        # 7. 清理临时文件
        rm -f '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam'
        rm -f '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam'*
        rm -f '$gatk_vcf_gz' '$gatk_vcf_gz.tbi'
        rm -f '$temp_vcf'
        rm -rf '$tmp_dir'
        
        echo '处理完成' | tee -a '$log_file'
EOT
    # 注意：上面的EOT必须顶格，且前面不能有任何字符
    
    echo "样本 $sample 任务已提交"
done

echo -e "\n所有样本任务提交完毕"

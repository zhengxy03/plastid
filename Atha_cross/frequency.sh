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
    
    # 定义文件路径
    bam_file="$PROJECT_ROOT/$sample/3_bwa/R.sort.bam"
    bam_index="$PROJECT_ROOT/$sample/3_bwa/R.sort.bai"
    gatk_vcf="$PROJECT_ROOT/$sample/3_gatk/R.filtered.vcf"
    gatk_vcf_gz="$PROJECT_ROOT/$sample/3_gatk/R.filtered.vcf.gz"
    result_file="$PROJECT_ROOT/$sample/chloroplast_heteroplasmy_freq.txt"
    chloro_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_from_gatk.vcf.gz"
    chloro_vcf_tbi="$chloro_vcf.tbi"
    temp_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_temp.vcf"
    filtered_temp_vcf="$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_filtered_temp.vcf"
    log_file="$PROJECT_ROOT/$sample/chloro_extraction.log"
    
    # 检查必要文件
    if [ ! -f "$bam_file" ] || [ ! -f "$bam_index" ] || [ ! -f "$gatk_vcf" ]; then
        echo "警告: 缺少必要文件（BAM/索引/GATK VCF），跳过样本 $sample"
        continue
    fi

    # 检查任务是否已在运行
    if bjobs -w | grep -q "${sample}_chloro"; then
        echo "任务已存在（${sample}_chloro），跳过: $sample"
        continue
    fi
    
    # 创建临时目录
    tmp_dir="$PROJECT_ROOT/$sample/tmp_scripts"
    mkdir -p "$tmp_dir"
    > "$log_file"
    
    # 1️⃣ 覆盖度计算脚本
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
    
    # 2️⃣ 多等位点频率计算脚本（异质性频率 ≥ 1%）
    cat > "$tmp_dir/calc_frequency.pl" << 'EOF'
#!/usr/bin/perl
use strict;
use warnings;

while (<STDIN>) {
    chomp;
    my @f = split(/\t/);
    next unless @f >= 6;

    my ($chrom, $pos, $ref, $alt_str, $ad_str, $qual) = @f[0..5];
    $qual = 0 unless (defined $qual && $qual =~ /^[-+]?\d+(\.\d+)?$/);

    my @alts = split(/,/, $alt_str // "");
    my @ad = split(/,/, $ad_str // "");

    next unless @ad == @alts + 1;   # 确保 AD 字段长度正确

    my $ref_cnt = $ad[0] // 0;
    my $total = $ref_cnt + 0;
    $total += $_ for @ad[1..$#ad];
    next if $total == 0;

    for my $i (0..$#alts) {
        my $alt_cnt = $ad[$i+1] // 0;
        my $freq = $alt_cnt / $total;
        next if $freq < 0.01;   # 过滤 <1%
        printf "%s\t%s\t%s\t%s\t%.6f\t%d\n", $chrom, $pos, $ref, $alts[$i], $freq, $total;
    }
}
EOF
    chmod +x "$tmp_dir/calc_frequency.pl"
    
    # 3️⃣ 提交 LSF 任务
    bsub -q mpi -n 24 -J "${sample}_chloro" <<EOT
        echo '开始处理样本 $sample' | tee -a '$log_file'
        
        # 提取叶绿体区域reads
        samtools view -b -q 10 -o '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam' '$bam_file' '$CHLORO_NAME'
        if [ ! -s '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam' ]; then
            echo '警告: 无有效叶绿体reads' | tee -a '$log_file'
            touch '$result_file'
            exit 0
        fi
        
        samtools sort -o '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam' '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam'
        samtools index '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam'
        
        avg_depth=\$(samtools depth '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam' | '$tmp_dir/calc_coverage.pl')
        echo "叶绿体平均覆盖度: \$avg_depth×" | tee -a '$log_file'
        
        # 压缩并索引 GATK VCF
        bgzip -c '$gatk_vcf' > '$gatk_vcf_gz'
        tabix -p vcf '$gatk_vcf_gz'
        
        # 提取叶绿体变异（PASS 标签）
        echo '提取叶绿体变异（PASS 标签）...' | tee -a '$log_file'
        bcftools view --apply-filters PASS --targets '$CHLORO_NAME' -Oz '$gatk_vcf_gz' > '$filtered_temp_vcf.gz'
        tabix -p vcf '$filtered_temp_vcf.gz'
        
        # 解压供下游 Perl 处理
        bcftools view '$filtered_temp_vcf.gz' > '$temp_vcf'
        
        # 统计变异总数（不含头部）
        var_count=\$(grep -v "^#" '$temp_vcf' | wc -l)
        echo "提取到的叶绿体变异总数（PASS）: \$var_count" | tee -a '$log_file'
        
        if [ \$var_count -eq 0 ]; then
            echo '未检测到任何变异' | tee -a '$log_file'
            touch '$result_file'
            exit 0
        fi
        
        # 压缩最终叶绿体VCF
        bgzip -c '$temp_vcf' > '$chloro_vcf'
        tabix -p vcf '$chloro_vcf'
        echo '生成最终叶绿体VCF完成' | tee -a '$log_file'
        
        # 计算异质性频率，并过滤 <1%
        echo '计算异质性频率并过滤小于1%位点...' | tee -a '$log_file'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t%QUAL\n' '$chloro_vcf' | \
        '$tmp_dir/calc_frequency.pl' > '$result_file'
        
        result_count=\$(wc -l '$result_file' | cut -d' ' -f1)
        echo "保留异质性变异数 (≥1%): \$result_count" | tee -a '$log_file'
        
        # 清理临时文件
        rm -f '$PROJECT_ROOT/$sample/3_bwa/chloroplast_temp.bam'

        rm -f '$gatk_vcf_gz' '$gatk_vcf_gz.tbi'
        rm -f '$temp_vcf' '$filtered_temp_vcf.gz' '$filtered_temp_vcf.gz.tbi'
        rm -rf '$tmp_dir'
        
        echo '样本 $sample 处理完成' | tee -a '$log_file'
EOT

    echo "样本 $sample 任务已提交（任务名：${sample}_chloro）"
done

echo -e "\n所有样本任务提交完毕！可通过 bjobs 查看任务状态。"

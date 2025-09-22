#!/bin/bash
set -euo pipefail

# 手动定义正确的项目根目录（根据你的实际路径）
PROJECT_ROOT="/share/home/wangq/zxy/plastid/Atha_cross"

# 强制验证项目根目录是否存在
if [ ! -d "$PROJECT_ROOT" ]; then
    echo "严重错误: 项目根目录不存在 - $PROJECT_ROOT"
    echo "请检查路径是否正确，当前用户目录为: $(pwd)"
    exit 1
fi

# 验证opts.tsv是否存在
opts_file="$PROJECT_ROOT/opts.tsv"
if [ ! -f "$opts_file" ]; then
    echo "严重错误: opts.tsv文件不存在于项目目录"
    echo "查找路径: $opts_file"
    echo "项目目录下的文件列表:"
    ls -l "$PROJECT_ROOT" | grep "opts.tsv"  # 显示是否有类似文件
    exit 1
fi

# 从opts.tsv提取样本名（严格处理）
samples=$(awk -F'\t' '
    NR > 0 && $1 != "" {          # 只处理非空行
        gsub(/^[ \t\r]+|[ \t\r]+$/, "");  # 清除前后空格和Windows换行符
        print $1
    }
' "$opts_file")

# 检查样本提取结果
if [ -z "$samples" ]; then
    echo "错误: 从opts.tsv中未提取到任何有效样本名"
    echo "opts.tsv的前5行内容:"
    head -n 5 "$opts_file"
    exit 1
fi

# 循环处理样本
for sample in $samples; do
    echo -e "\n===== 处理样本: $sample ====="
    
    # 定义完整路径
    bam_file="$PROJECT_ROOT/$sample/3_bwa/R.sort.bam"
    bam_index="$PROJECT_ROOT/$sample/3_bwa/R.sort.bai"
    gatk_vcf="$PROJECT_ROOT/$sample/3_gatk/R.raw.vcf"
    gatk_vcf_idx="$PROJECT_ROOT/$sample/3_gatk/R.raw.vcf.idx"
    
    # 显示当前检查的路径（关键调试信息）
    echo "当前检查的文件路径:"
    echo "BAM文件: $bam_file"
    echo "BAM索引: $bam_index"
    echo "GATK VCF: $gatk_vcf"
    
    # 检查BAM文件
    if [ ! -f "$bam_file" ]; then
        echo "警告: BAM文件不存在 - $bam_file"
        continue
    fi
    if [ ! -f "$bam_index" ]; then
        echo "警告: BAM索引不存在 - $bam_index"
        continue
    fi
    
    # 检查GATK文件
    if [ ! -f "$gatk_vcf" ]; then
        echo "警告: GATK VCF不存在 - $gatk_vcf"
        continue
    fi
    if [ ! -f "$gatk_vcf_idx" ]; then
        echo "警告: GATK VCF索引不存在 - $gatk_vcf_idx"
        continue
    fi
    
    # 检查结果文件
    result_file="$PROJECT_ROOT/$sample/chloroplast_heteroplasmy_freq.txt"

    
    # 检查任务是否存在
    if bjobs -w | grep -q "${sample}_chloro"; then
        echo "任务已存在，跳过: $sample"
        continue
    fi
    
    # 创建临时脚本目录
    tmp_dir="$PROJECT_ROOT/$sample/tmp_scripts"
    mkdir -p "$tmp_dir"
    
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
    my @ad = split(/,/, $ad_str // "0,0");
    my ($ref_cnt, $alt_cnt) = ($ad[0] // 0, $ad[1] // 0);
    my $total = $ref_cnt + $alt_cnt;
    
    if ($total >= 5 && $qual >= 5 && $alt_cnt >= 1) {
        my $freq = $alt_cnt / $total;
        if ($freq >= 0.001 && $freq <= 0.999) {
            printf "%s\t%s\t%s\t%s\t%.6f\t%d\n", 
                $chrom, $pos, $ref, $alt, $freq, $total;
        }
    }
}
EOF
    chmod +x "$tmp_dir/calc_frequency.pl"
    
    # 提交任务
    bsub -q mpi -n 24 -J "${sample}_chloro" "
        echo '开始处理样本 $sample'
        
        # 提取叶绿体BAM
        samtools view -b -q 10 -o '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt.bam' '$bam_file' 'Pt'
        samtools sort -o '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam' '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt.bam'
        samtools index '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam'
        
        # 计算覆盖度
        avg_depth=\$(samtools depth '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam' | '$tmp_dir/calc_coverage.pl')
        echo \"$sample 叶绿体平均覆盖度: \$avg_depth×\" >> '$PROJECT_ROOT/$sample/chloro_coverage.log'
        
        # 提取叶绿体变异
        bcftools view -r 'Pt' '$gatk_vcf' > '$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_from_gatk.vcf'
        
        # 计算频率
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AD\t%QUAL\n' '$PROJECT_ROOT/$sample/3_bwa/chloroplast_variants_from_gatk.vcf' | '$tmp_dir/calc_frequency.pl' > '$result_file'
        
        # 清理
        rm -f '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt.bam'
        rm -f '$PROJECT_ROOT/$sample/3_bwa/chloroplast_Pt_sorted.bam'*
        rm -rf '$tmp_dir'
        
        echo '样本 $sample 处理完成'
    "
    
    echo "样本 $sample 任务已提交"
done

echo -e "\n所有样本处理完成"

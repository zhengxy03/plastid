#!/usr/bin/perl
use strict;
use warnings;

# 适配TAIR10格式的TRF输出解析（染色体命名如@1、@2...）
my $current_chrom = "";
my $total_entries = 0;
my $written_entries = 0;

print STDERR "开始解析TRF输出文件...\n";

while (my $line = <>) {
    chomp $line;
    $total_entries++;
    
    # 识别染色体名称行（TAIR10格式以@开头，如@1、@2）
    if ($line =~ /^@(\d+)\s/) {  # 只匹配@1、@2、@3、@4、@5
        $current_chrom = $1;
        print STDERR "检测到核染色体: $current_chrom\n";
        next;
    }
    
    # 跳过细胞器染色体（@Mt、@Pt）和其他注释行
    if ($line =~ /^@(Mt|Pt|Un|C|scaffold)/) {
        $current_chrom = "";  # 清空非核染色体
        next;
    }
    
    # 跳过空行
    next if $line =~ /^\s*$/;
    
    # 解析TRF结果行
    my @fields = split(/\s+/, $line);
    
    # 检查字段数量（TRF结果行通常有13-15个字段）
    if (@fields < 13) {
        print STDERR "警告：行 $total_entries 字段数量不足（".@fields."个），跳过\n";
        next;
    }
    
    # 提取坐标信息（TAIR10的TRF输出前两个字段是起始和结束位置）
    my ($start, $end) = @fields[0, 1];
    
    # 验证坐标有效性
    if ($start !~ /^\d+$/ || $end !~ /^\d+$/ || $start >= $end) {
        print STDERR "警告：行 $total_entries 坐标无效（$start-$end），跳过\n";
        next;
    }
    
    # 只处理核染色体（1-5）的重复序列
    if ($current_chrom =~ /^[1-5]$/) {
        # TRF坐标是1-based，转换为BED的0-based
        my $bed_start = $start - 1;
        my $bed_end = $end;
        
        # 提取重复单元信息
        my $period = $fields[2];
        my $copy_number = $fields[3];
        my $consensus = $fields[10] // "N/A";
        
        # 输出BED格式
        print join("\t", 
            $current_chrom,
            $bed_start,
            $bed_end,
            "TRF:".substr($consensus, 0, 20),
            $period,
            sprintf("%.2f", $copy_number)
        )."\n";
        
        $written_entries++;
    }
}

print STDERR "\n解析完成：\n";
print STDERR "总处理行数：$total_entries\n";
print STDERR "写入核染色体串联重复数：$written_entries\n";
if ($written_entries == 0) {
    print STDERR "警告：未找到核染色体的串联重复序列，请检查文件格式\n";
} else {
    print STDERR "结果已成功写入输出文件\n";
}

exit 0;
    

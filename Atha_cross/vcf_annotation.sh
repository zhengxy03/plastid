#!/bin/bash
set -euo pipefail

# ====================== 配置参数 ======================
GFF_GZ="Arabidopsis_thaliana.TAIR10.62.gff3.gz"  # GFF文件路径
HETERO_CSV="chloro_multi_stats.csv"               # 输入CSV
LOCUS_COL=8                                       # 位点列索引（第9列）
SAMPLE_COL=0                                      # 样本名列索引（第1列）
# =======================================================

TEMP_FILE="temp_annotated.csv"
FINAL_OUTPUT="final_annotated_english.csv"

# 检查是否安装Perl CSV模块
if ! perl -MText::CSV -e 1 2>/dev/null; then
    echo "请先安装Perl CSV模块：sudo cpan Text::CSV"
    exit 1
fi

# 步骤1：筛选叶绿体注释
echo "=== Filter chloroplast annotations ==="
gunzip -c "$GFF_GZ" | perl -ne '
    chomp;
    next if /^\s*#/ || /^\s*$/;
    my @cols = split(/\t/);
    print "$_\n" if $cols[0] eq "Pt" && $cols[2] =~ /^gene|CDS|tRNA|rRNA$/;
' > "chloro_Pt.gff3"

# 步骤2：构建基因索引
echo -e "\n=== Build gene index ==="
perl -ne '
    chomp;
    my @cols = split(/\t/);
    my ($type, $start, $end, $strand, $attrs) = @cols[2,3,4,6,8];
    my ($gene_id, $gene_name, $product) = ("Unknown_ID", "Unknown_gene", "No_description");
    $gene_id = $1 if $attrs =~ /ID=([^;]+)/;
    $gene_name = $1 if $attrs =~ /Name=([^;]+)/;
    $product = $1 if $attrs =~ /product=([^;]+)/;
    if ($type eq "gene" || ($type =~ /tRNA|rRNA/ && !exists $seen{$gene_id})) {
        $seen{$gene_id} = 1;
        print join("\t", "$start-$end", $gene_name, $product, $strand, $gene_id) . "\n" if $start <= $end;
    }
' "chloro_Pt.gff3" > "chloro_gene_index.txt"

# 步骤3：注释表格（全英文）
echo -e "\n=== Annotate table ==="
# 处理表头
perl -MText::CSV -e '
    my $csv = Text::CSV->new({ binary => 1, auto_diag => 1 });
    open my $fh, "<", $ARGV[0] or die "$ARGV[0]: $!";
    my $header = $csv->getline($fh);
    # 添加英文注释列
    push @$header, qw(Associated_gene Gene_description Strand Gene_ID Gene_coordinate_range Chloroplast_region);
    print join(",", @$header) . "\n";
    close $fh;
' "$HETERO_CSV" > "$TEMP_FILE"

# 处理数据行
perl -MText::CSV -e '
    my $locus_col = '$LOCUS_COL';
    my $csv = Text::CSV->new({ binary => 1, auto_diag => 1 });
    
    # 加载基因索引
    my %gene;
    open my $gh, "<", "chloro_gene_index.txt" or die "索引文件错误：$!";
    while (<$gh>) {
        chomp;
        my ($range, $name, $def, $strand, $id) = split(/\t/, $_, 5);
        $gene{$range} = [ $name, $def, $strand, $id ];
    }
    close $gh;
    
    # 解析CSV数据行
    open my $fh, "<", $ARGV[0] or die "$ARGV[0]: $!";
    my $header = $csv->getline($fh);  # 跳过表头
    
    while (my $row = $csv->getline($fh)) {
        next if @$row <= $locus_col;
        
        # 提取位点
        my $loc_str = $row->[$locus_col] // "";
        my ($loc_pos) = $loc_str =~ /(\d+)/;
        next unless defined $loc_pos;
        
        # 匹配基因注释
        my ($name, $def, $strand, $id, $range) = ("Non-coding", "No_associated_gene", "None", "None", "None");
        foreach my $r (keys %gene) {
            my ($s, $e) = split(/-/, $r);
            if ($loc_pos >= $s && $loc_pos <= $e) {
                ($name, $def, $strand, $id) = @{$gene{$r}};
                $range = $r;
                last;
            }
        }
        
        # 叶绿体区域（英文）
        my $large_region;
        if ($loc_pos >=1 && $loc_pos <=84079)      { $large_region = "LSC (Large Single Copy)"; }
        elsif ($loc_pos >=84080 && $loc_pos <=110084) { $large_region = "IRb (Inverted Repeat b)"; }
        elsif ($loc_pos >=110085 && $loc_pos <=122218) { $large_region = "SSC (Small Single Copy)"; }
        else                                      { $large_region = "IRa (Inverted Repeat a)"; }
        
        # 添加注释并输出
        push @$row, $name, $def, $strand, $id, $range, $large_region;
        print join(",", @$row) . "\n";
    }
    close $fh;
' "$HETERO_CSV" >> "$TEMP_FILE"

# 步骤4：去重、移除样本名
echo -e "\n=== Process final format ==="
perl -MText::CSV -e '
    my %seen_locus;
    my $locus_col = '$LOCUS_COL';
    my $sample_col = '$SAMPLE_COL';
    my $csv = Text::CSV->new({ binary => 1, auto_diag => 1 });
    
    open my $fh, "<", $ARGV[0] or die "$ARGV[0]: $!";
    my $header = $csv->getline($fh);
    
    # 处理表头（移除样本名）
    splice(@$header, $sample_col, 1);
    unshift @$header, "Highest_frequency_locus";  # 英文表头
    print join(",", @$header) . "\n";
    
    # 处理数据行
    while (my $row = $csv->getline($fh)) {
        my $loc_str = $row->[$locus_col] // "";
        my ($locus) = $loc_str =~ /(\d+)/;
        next unless defined $locus;
        
        if (!$seen_locus{$locus}++) {
            splice(@$row, $sample_col, 1);  # 移除样本名
            unshift @$row, $locus;          # 位点移到第一列
            print join(",", @$row) . "\n";
        }
    }
    close $fh;
' "$TEMP_FILE" > "$FINAL_OUTPUT"

# 清理临时文件
rm "$TEMP_FILE" "chloro_Pt.gff3" "chloro_gene_index.txt"

echo -e "\n✅ Processing completed! Final file: $FINAL_OUTPUT"
head -5 "$FINAL_OUTPUT"

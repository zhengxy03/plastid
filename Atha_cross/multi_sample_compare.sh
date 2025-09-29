#!/bin/bash
set -euo pipefail

# ====================== 1. 配置参数（根据你的文件路径修改） ======================
GFF_GZ="Arabidopsis_thaliana.TAIR10.62.gff3.gz"  # 你的TAIR10 GFF压缩文件
HETERO_CSV="chloro_multi_stats.csv"               # 你的异质性分析表格（CSV格式！）
CHLORO_GFF="chloro_TAIR10.gff3"                  # 临时输出：叶绿体GFF
GENE_INDEX="chloro_gene_index.txt"                # 临时输出：基因索引
ANNOT_OUTPUT="annotated_heteroplasmy.csv"         # 最终输出：带注释的CSV
LOCUS_COL=8                                       # “最高频率位点”列的索引（第9列=8）
# ==============================================================================

# ====================== 2. 步骤1：解压并筛选叶绿体GFF（修复：用NC_000932.1匹配） ======================
echo "=== 步骤1：筛选叶绿体GFF注释 ==="
# 解压GFF → 筛选seqid=NC_000932.1（叶绿体）且类型为gene/CDS/tRNA/rRNA的行
gunzip -c "$GFF_GZ" | perl -ne '
    chomp;
    next if /^\s*#/ || /^\s*$/;  # 跳过注释行和空行
    my @cols = split(/\t/);       # GFF本身是TSV格式，仍按制表符分割
    my $seqid = $cols[0] // "";
    my $type = $cols[2] // "";
    # 关键修复：用NC_000932.1匹配叶绿体（替代原ChrC）
    print "$_\n" if $seqid eq "NC_000932.1" && $type =~ /^gene|CDS|tRNA|rRNA$/;
' > "$CHLORO_GFF"

# 检查筛选结果，若仍为0行则提示
CHLORO_ROWS=$(wc -l < "$CHLORO_GFF")
if [ "$CHLORO_ROWS" -eq 0 ]; then
    echo "警告：叶绿体GFF筛选结果为0行！可能GFF中叶绿体ID不是NC_000932.1，建议用以下命令查看GFF中的seqid："
    echo "gunzip -c $GFF_GZ | grep -v '^#' | cut -f1 | sort -u | head -20"
    exit 1
fi
echo "筛选完成！叶绿体GFF共 $CHLORO_ROWS 行，保存至：$CHLORO_GFF"


# ====================== 3. 步骤2：从叶绿体GFF构建基因索引（不变，GFF是TSV） ======================
echo -e "\n=== 步骤2：构建基因坐标索引 ==="
perl -ne '
    chomp;
    my @cols = split(/\t/);       # GFF是TSV，按制表符分割
    my ($type, $start, $end, $strand, $attrs) = @cols[2,3,4,6,8];
    
    # 解析属性列：提取基因ID、基因名、功能描述
    my ($gene_id, $gene_name, $product) = ("未知ID", "未知基因", "无功能描述");
    $gene_id = $1 if $attrs =~ /ID=([^;]+)/;
    $gene_name = $1 if $attrs =~ /Name=([^;]+)/;
    $product = $1 if $attrs =~ /product=([^;]+)/;
    
    # 优先保留gene类型，避免CDS重复
    if ($type eq "gene" || ($type =~ /tRNA|rRNA/ && !exists $seen{$gene_id})) {
        $seen{$gene_id} = 1;
        print join("\t", "$start-$end", $gene_name, $product, $strand, $gene_id) . "\n"
            if $start <= $end;
    }
' "$CHLORO_GFF" > "$GENE_INDEX"

INDEX_COUNT=$(wc -l < "$GENE_INDEX")
echo "索引构建完成！共 $INDEX_COUNT 个叶绿体基因，保存至：$GENE_INDEX"


# ====================== 4. 步骤3：注释异质性表格（修复：按逗号分割CSV） ======================
echo -e "\n=== 步骤3：注释异质性表格 ==="
# 处理表头：按逗号分割CSV，添加注释列
head -n1 "$HETERO_CSV" | perl -ne '
    chomp;
    my @cols = split(/,/);        # 关键修复：CSV按逗号分割
    # 新增注释列（CSV格式，用逗号分隔）
    push @cols, qw(关联基因名,基因功能描述,基因链方向,基因TAIR ID,基因坐标范围,叶绿体大区域);
    print join(",", @cols) . "\n";
' > "$ANNOT_OUTPUT"

# 处理表格内容：按逗号分割CSV，匹配注释
tail -n+2 "$HETERO_CSV" | perl -sne '
    BEGIN {
        # 加载基因索引（索引是TSV，按制表符分割）
        open(my $fh, "<", $gene_index) or die "无法打开索引：$!";
        while (<$fh>) {
            chomp;
            my ($range, $name, $def, $strand, $id) = split(/\t/, $_, 5);
            $gene{$range} = join(",", $name, $def, $strand, $id);  # 用逗号拼接，便于后续分割
        }
        close($fh);
    }

    chomp;
    my @cols = split(/,/);        # 关键修复：CSV按逗号分割
    my $locus_str = $cols[$locus_col] // "";  # 提取“最高频率位点”列（如111，无需Pt:前缀）
    die "缺少位点列：$_" if $locus_str eq "";

    # 直接使用位点坐标（你的CSV中是纯数字，如111，无需分割Pt:）
    my $loc_pos = $locus_str;
    die "位点格式错误（需为数字）：$locus_str" unless $loc_pos =~ /^\d+$/;

    # 匹配基因信息
    my ($name, $def, $strand, $id, $range) = ("非编码区", "无关联基因", "无", "无", "无");
    foreach my $r (keys %gene) {
        my ($s, $e) = split(/-/, $r);
        if ($loc_pos >= $s && $loc_pos <= $e) {
            ($name, $def, $strand, $id) = split(/,/, $gene{$r});  # 按逗号分割基因信息
            $range = $r;
            last;
        }
    }

    # 判断叶绿体大区域（TAIR10坐标）
    my $large_region;
    if ($loc_pos >=1 && $loc_pos <=84079)      { $large_region = "LSC区（大单拷贝区）"; }
    elsif ($loc_pos >=84080 && $loc_pos <=110084) { $large_region = "IRb区（反向重复区b）"; }
    elsif ($loc_pos >=110085 && $loc_pos <=122218) { $large_region = "SSC区（小单拷贝区）"; }
    else                                      { $large_region = "IRa区（反向重复区a）"; }

    # 合并注释并输出（CSV格式，用逗号分隔）
    push @cols, $name, $def, $strand, $id, $range, $large_region;
    print join(",", @cols) . "\n";
' -- -gene_index="$GENE_INDEX" -locus_col="$LOCUS_COL" >> "$ANNOT_OUTPUT"

echo "注释完成！带注释的CSV表格保存至：$ANNOT_OUTPUT"
echo -e "\n=== 结果预览（前5行） ==="
head -5 "$ANNOT_OUTPUT"
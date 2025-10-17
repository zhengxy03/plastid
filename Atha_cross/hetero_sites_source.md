## calculate the proportion of sources for variant sites
* find parents' GT difference
```
mkdir -p analysis
cd analysis

#提取亲本
cp ../Sample_Col_G/3_gatk/R.filtered.vcf Col.vcf
cp ../Sample_Ler_XL_4/3_gatk/R.filtered.vcf Ler.vcf

bgzip Col.vcf Ler.vcf

bcftools view --apply-filters PASS --max-alleles 2 --targets Pt  --include "AF>0.01" -Oz Col.vcf.gz > Col_Pt.vcf.gz
bcftools view --apply-filters PASS --max-alleles 2 --targets Pt  --include "AF>0.01" -Oz Ler.vcf.gz > Ler_Pt.vcf.gz
bcftools index Col_Pt.vcf.gz
bcftools index Ler_Pt.vcf.gz

#亲本差异位点
bcftools reheader -s <(echo "Col") Col_Pt.vcf.gz -o Col_Pt.renamed.vcf.gz
bcftools reheader -s <(echo "Ler") Ler_Pt.vcf.gz -o Ler_Pt.renamed.vcf.gz
bcftools index Col_Pt.renamed.vcf.gz
bcftools index Ler_Pt.renamed.vcf.gz

bcftools merge -m all Col_Pt.renamed.vcf.gz Ler_Pt.renamed.vcf.gz -Oz -o Col_Ler.merged.vcf.gz
bcftools index Col_Ler.merged.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' Col_Ler.merged.vcf.gz > Col_Ler_GT.tsv
perl -ane 'print if $F[4] ne $F[5]' Col_Ler_GT.tsv > different_sites.tsv
```
* extract progeny GT
```
cd ..
mkdir -p vcf
cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -d "{1}" ] || [ ! -f "{1}/3_gatk/R.filtered.vcf" ]; then
            echo "样本 {1} 文件缺失，跳过"
            continue
        fi

        echo "处理样本: {1}"
        
        # 处理VCF：重命名样本、筛选Pt、过滤PASS、双等位基因、最小AF
        bcftools reheader --samples <(echo {1}) "{1}/3_gatk/R.filtered.vcf" |
            bcftools view \
                --apply-filters PASS \
                --max-alleles 2 \
                --targets Pt \
                -Oz |
            bcftools view --include "AF>0.01" -Oz -o vcf/{1}.vcf.gz

        # 索引VCF
        bcftools index -f vcf/{1}.vcf.gz
    '

bcftools merge --merge all -l <(
        cat opts.tsv |
            cut -f 1 |
            grep -Ev "^Sample_Col_G$|^Sample_Ler_XL_4$" |
            parallel -k -j 1 ' [ -f vcf/{}.vcf.gz ] && echo "vcf/{}.vcf.gz" '
    ) \
    > Atha_cross.vcf

bcftools stats Atha_cross.vcf > Atha_cross.stats
plot-vcfstats  Atha_cross.stats  -p  plots/Atha_cross.stats

cd analysis
bgzip ../Atha_cross.vcf
bcftools index -f ../Atha_cross.vcf.gz

#与父母本比较
bcftools view -R <(cut -f1,2 different_sites.tsv) ../Atha_cross.vcf.gz -Oz -o F2_informative.vcf.gz

#排除Mt\Nc中的相似性位点
cp ~/data/plastid/genome/col0/genome.fa .
faops some genome.fa <(echo 1; echo 2; echo 3; echo 4; echo 5) chr.fa
faops some genome.fa <(echo Pt) Pt.fa
faops some genome.fa <(echo Mt) Mt.fa

bwa index chr.fa
bwa index Mt.fa
bwa index Pt.fa

bwa mem -t 8 chr.fa Pt.fa | samtools sort -o Pt_vs_nuclear.bam
bwa mem -t 8 Mt.fa Pt.fa  | samtools sort -o Pt_vs_mt.bam

bedtools bamtobed -i Pt_vs_nuclear.bam > Pt_vs_nuclear.bed
bedtools bamtobed -i Pt_vs_mt.bam > Pt_vs_mt.bed
cat Pt_vs_nuclear.bed Pt_vs_mt.bed | bedtools sort -i - | bedtools merge -i - > Pt_homology.bed
bcftools view -T ^Pt_homology.bed F2_informative.vcf.gz -Oz -o F2_filtered_homology.vcf.gz
bcftools index -f F2_filtered_homology.vcf.gz
```
* data filtering and processing
```
#稀有位点过滤 (<10% 样本)
bcftools query -f '%CHROM\t%POS[\t%GT]\n' F2_filtered_homology.vcf.gz |
perl -F'\t' -ane '
    $total_samples = '$(cat opts_no_parents.tsv | wc -l)';
    $count_alt = 0;
    $total = 0;
    for ($i = 2; $i < @F; $i++) {
        if ($F[$i] ne "./.") {
            $total++;
            $count_alt++ if $F[$i] =~ /1/;
        }
    }
    if ($total > 0 && $count_alt / $total >= 0.1) {
        print join("\t", @F[0,1]), "\n";
    }
' > common_sites.txt
bcftools view -R common_sites.txt F2_filtered_homology.vcf.gz -Oz -o F2_common.vcf.gz
bcftools index -f F2_common.vcf.gz


#成簇位点过滤 (<50bp)
# 提取并转换位点
bcftools query -f '%CHROM\t%POS\n' F2_common.vcf.gz | \
perl -ane 'print "$F[0]\t" . ($F[1] - 1) . "\t$F[1]\n"' > sites.bed

# 排序和聚类
bedtools sort -i sites.bed | \
bedtools cluster -d 50 -i - > clustered.bed

#处理聚类结果
perl -ane '
    $cluster_count{$F[3]}++;
    push @{$cluster_info{$F[3]}}, [$F[0], $F[2]];
    END {
        foreach $cluster_id (sort { $a <=> $b } keys %cluster_count) {
            if ($cluster_count{$cluster_id} == 1) {
                my $site = $cluster_info{$cluster_id}[0];
                print "$site->[0]\t$site->[1]\n";
            }
        }
    }
' clustered.bed > nonclustered_chrom_pos.txt

bcftools view -R nonclustered_chrom_pos.txt F2_common.vcf.gz -Oz -o F2_noncluster.vcf.gz
bcftools index -f F2_noncluster.vcf.gz


#排除测序或者 PCR 的偏向性错误（单碱基重复±2bp过滤）
grep -v "^>" Pt.fa | tr -d '\n' > Pt.seq.txt

perl -ne '
    while (/(A{5,}|T{5,}|C{5,}|G{5,})/g) {
        $start = $-[0];
        $end   = $+[0]; 
        print "Pt\t$start\t$end\n";
    }
' Pt.seq.txt > homopolymer.bed
samtools faidx Pt.fa
bedtools slop -i homopolymer.bed -g Pt.fa.fai -b 2 > homopolymer_plus2.bed
bcftools view -T ^homopolymer_plus2.bed F2_noncluster.vcf.gz -Oz -o F2_highconf.vcf.gz
bcftools index -f F2_highconf.vcf.gz

# 提取 GT 矩阵（基于高置信度 VCF）
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' F2_highconf.vcf.gz > F2_GT_matrix.tsv

grep -v -E '^Sample_Col_G|^Sample_Ler_XL_4|^Sample_c1c2|^Sample_l2c2|^Sample_l2l3|^Sample_l4c1|^Sample_l4l3' ../opts.tsv > opts_no_parents.tsv
echo -e "CHROM\tPOS\t$(cat opts_no_parents.tsv | cut -f1 | paste -sd '\t' -)" > header.txt
cat header.txt F2_GT_matrix.tsv > F2_GT_matrix_with_header.tsv

bcftools query -f '%CHROM\t%POS\n' F2_highconf.vcf.gz > highconf_sites.list

perl -F'\s+' -ane '
    BEGIN {
        open my $fh, "<", "highconf_sites.list" or die $!;
        while (<$fh>) {
            chomp;
            my @f = split /\s+/;
            $seen{"$f[0]:$f[1]"} = 1;
        }
        close $fh;
    }
    $key = "$F[0]:$F[1]";
    print if $seen{$key};
' different_sites_with_len.tsv > different_sites_with_len_highconf_perl.tsv

perl -F'\t' -ane '
    BEGIN {
        open my $fh, "<", "highconf_sites.list" or die $!;
        while (<$fh>) {
            chomp;
            $seen{$_} = 1;
        }
        close $fh;
    }
    $key = "$F[0]:$F[1]";
    print if $seen{$key};
' F2_sites_ref_alt.tsv > F2_sites_ref_alt_highconf_perl.tsv
```
* calculate the proportion
```
GENOME_SIZE=154478
total_samples=$(cat opts_no_parents.tsv | wc -l)
current=1

> sample_ratio_tmp.tsv
cat opts_no_parents.tsv | while read -r line; do
    sample=$(echo "$line" | cut -f1)
    target_col=$((current + 2))

    echo "[$current/$total_samples] 正在处理样本 $sample，列号 $target_col"

    sample="$sample" target_col="$target_col" GENOME_SIZE="$GENOME_SIZE" \
    perl -F'\t' -ane '
        BEGIN {
            $sample     = $ENV{"sample"};
            $target_col = $ENV{"target_col"};
            $genome_size= $ENV{"GENOME_SIZE"};

            # 父母差异位点长度
            open my $fh1, "<", "different_sites_with_len_highconf.tsv" or die $!;
            while (<$fh1>) {
                chomp;
                my @f = split /\t/;
                my $key = "$f[0]:$f[1]";
                $col_gt{$key}  = $f[4];
                $ler_gt{$key}  = $f[5];
                $site_len{$key}= $f[6];
            }
            close $fh1;

            # F2位点长度
            open my $fh2, "<", "F2_sites_ref_alt_highconf.tsv" or die $!;
            while (<$fh2>) {
                chomp;
                my @f = split /\t/;
                my $key = "$f[0]:$f[1]";
                my $len_ref = length($f[2]);
                my $len_alt = length($f[3]);
                $f2_len{$key} = ($len_ref == $len_alt) ? $len_ref : abs($len_ref - $len_alt);
            }
            close $fh2;

            $count_col = $count_ler = $count_mut = $total = 0;
            $len_col = $len_ler = $len_mut = 0;
        }

        next if $. == 1;   # 跳过header

        my $key = "$F[0]:$F[1]";
        my $gt  = $F[$target_col - 1];
        my $col = $col_gt{$key};
        my $ler = $ler_gt{$key};

        my $this_len = exists $site_len{$key} ? $site_len{$key} : $f2_len{$key};
        next if $gt eq "./.";

        if ($gt eq $col) {
            $count_col++; $len_col += $this_len;
        } elsif ($gt eq $ler) {
            $count_ler++; $len_ler += $this_len;
        } else {
            $count_mut++; $len_mut += $this_len;
        }
        $total++;

        END {
            if ($total > 0 && $genome_size > 0) {
                my $col_ratio = $count_col / $total * 100;
                my $ler_ratio = $count_ler / $total * 100;
                my $mut_ratio = $count_mut / $total * 100;

                my $col_genome = $len_col / $genome_size * 100;
                my $ler_genome = $len_ler / $genome_size * 100;
                my $mut_genome = $len_mut / $genome_size * 100;

                printf("%s\t%.4f\t%.6f\t%.4f\t%.6f\t%.4f\t%.6f\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    $sample, $col_ratio, $col_genome,
                    $ler_ratio, $ler_genome, $mut_ratio, $mut_genome,
                    $count_col, $count_ler, $count_mut,
                    $len_col, $len_ler, $len_mut);
            } else {
                printf("%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n", $sample);
            }
        }
    ' F2_GT_matrix_with_header.tsv >> sample_ratio_tmp.tsv

    current=$((current + 1))
done

echo -e "Sample\tColRatio\tColGenome\tLerRatio\tLerGenome\tMutRatio\tMutGenome\tColCount\tLerCount\tMutCount\tColLength\tLerLength\tMutLength" > sample_ratio_summary.tsv
cat sample_ratio_tmp.tsv >> sample_ratio_summary.tsv

echo -e "| 样本名称 | Col型比例(%) | Col基因组占比(%) | Ler型比例(%) | Ler基因组占比(%) | Mut型比例(%) | Mut基因组占比(%) | Col型位点数 | Ler型位点数 | Mut型位点数 | Col长度 | Ler长度 | Mut长度 |" > sample_ratio_summary.md
echo -e "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |" >> sample_ratio_summary.md

tail -n +2 sample_ratio_summary.tsv | awk -F'\t' '{
    printf("| %s | %s | %s | %s | %s | %s | %s | %d | %d | %d | %d | %d | %d |\n", \
        $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13);
}' >> sample_ratio_summary.md

```
| 样本名称 | Col型比例(%) | Col基因组占比(%) | Ler型比例(%) | Ler基因组占比(%) | Mut型比例(%) | Mut基因组占比(%) | Col型位点数 | Ler型位点数 | Mut型位点数 | Col长度 | Ler长度 | Mut长度 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Sample_14 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_18 | 100.0000 | 0.002589 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 4 | 0 | 0 | 4 | 0 | 0 |
| Sample_19 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_20 | 100.0000 | 0.001295 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 2 | 0 | 0 | 2 | 0 | 0 |
| Sample_21 | 100.0000 | 0.000647 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 1 | 0 | 0 | 1 | 0 | 0 |
| Sample_4 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_5 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Sample_7 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_8 | 100.0000 | 0.001295 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 2 | 0 | 0 | 2 | 0 | 0 |
| Sample_c41 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_c42 | 100.0000 | 0.002589 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 4 | 0 | 0 | 4 | 0 | 0 |
| Sample_c45 | 66.6667 | 0.001295 | 0.0000 | 0.000000 | 33.3333 | 0.000647 | 2 | 0 | 1 | 2 | 0 | 1 |
| Sample_c47 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_c48 | 100.0000 | 0.000647 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 1 | 0 | 0 | 1 | 0 | 0 |
| Sample_c51 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_c52 | 100.0000 | 0.000647 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 1 | 0 | 0 | 1 | 0 | 0 |
| Sample_c54 | 100.0000 | 0.003237 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 5 | 0 | 0 | 5 | 0 | 0 |
| Sample_c57 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_c61 | 100.0000 | 0.001295 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 2 | 0 | 0 | 2 | 0 | 0 |
| Sample_c62 | 100.0000 | 0.002589 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 4 | 0 | 0 | 4 | 0 | 0 |
| Sample_c63 | 100.0000 | 0.000647 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 1 | 0 | 0 | 1 | 0 | 0 |
| Sample_c64 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_c65 | 100.0000 | 0.001295 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 2 | 0 | 0 | 2 | 0 | 0 |
| Sample_c66 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Sample_c73 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_c81 | 75.0000 | 0.001942 | 0.0000 | 0.000000 | 25.0000 | 0.000647 | 3 | 0 | 1 | 3 | 0 | 1 |
| Sample_c82 | 100.0000 | 0.001295 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 2 | 0 | 0 | 2 | 0 | 0 |
| Sample_c83 | 100.0000 | 0.000647 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 1 | 0 | 0 | 1 | 0 | 0 |
| Sample_c84 | 100.0000 | 0.000647 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 1 | 0 | 0 | 1 | 0 | 0 |
| Sample_c85 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_c87 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_c88 | 100.0000 | 0.003237 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 4 | 0 | 0 | 5 | 0 | 0 |
| Sample_c89 | 100.0000 | 0.001942 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 3 | 0 | 0 | 3 | 0 | 0 |
| Sample_c90 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Sample_c91 | 100.0000 | 0.000647 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 1 | 0 | 0 | 1 | 0 | 0 |
| Sample_c92 | 100.0000 | 0.000647 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 1 | 0 | 0 | 1 | 0 | 0 |
| Sample_c93 | 100.0000 | 0.002589 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 4 | 0 | 0 | 4 | 0 | 0 |
| Sample_c94 | 100.0000 | 0.001295 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 2 | 0 | 0 | 2 | 0 | 0 |
| Sample_c95 | 100.0000 | 0.003237 | 0.0000 | 0.000000 | 0.0000 | 0.000000 | 4 | 0 | 0 | 5 | 0 | 0 |
```
library(ggplot2)
library(tidyr)
library(dplyr)

df <- read.table("sample_ratio_summary.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

df_long <- df %>%
  select(Sample, ColGenome, LerGenome) %>%
  pivot_longer(cols = c(ColGenome, LerGenome), names_to = "来源", values_to = "占比") %>%
  filter(!is.na(占比))

df_long$来源 <- recode(df_long$来源,
                        ColGenome = "Col0",
                        LerGenome = "Ler")

p <- ggplot(df_long, aes(x=来源, y=占比)) +
  geom_jitter(width=0.2, size=1, color="black") +
  stat_summary(fun=mean, geom="crossbar", width=0.5, fatten=2, color="red") +
  labs(x="Source", y="Heterozygosity Ratio (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    axis.ticks = element_line(color="black"),
    axis.text = element_text(color="black"),
    axis.title = element_text(color="black"),
    plot.background = element_rect(fill="white", color=NA),
    panel.background = element_rect(fill="white", color=NA)
  )

ggsave("genome_proportion_dotplot.png", p, width=6, height=4, dpi=300)
```
![Heterozygosity Ratio][def]

[def]: ../pic/genome_proportion_dotplot.png
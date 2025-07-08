# Progeny Chloroplast Heterozygosity
Use Atha_cross(*Arabidopsis thaliana* Columbia X Ler) as example
## Basic info
* Genome: GCF_000001735.3, TAIR10, 119.668 Mb<br>
* Chloroplast: NC_000932, Columbia, 154478 bp<br>
* Chloroplast: KX551970, Landsberg erecta, 154515 bp<br> Ler-0，另一种常用野生型，与 Col-0 存在遗传差异
* Mitochondrion: Y08501, 366924 bp<br>
## Project
* https://www.pnas.org/content/109/51/20992.long<br>
* PRJNA178613<br>
## Download
* Reference
```
mkdir -p ~/data/plastid/Atha_cross/genome
cd ~/data/plastid/Atha_cross/genome

for ACCESSION in "NC_000932" "Y08501"; do
    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
    curl $URL -o ${ACCESSION}.fa
done

TAB=$'\t'
cat <<EOF > replace.tsv
NC_000932${TAB}Pt
Y08501${TAB}Mt
EOF

cat NC_000932.fa Y08501.fa |
    faops filter -s stdin stdout |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(echo Pt; echo Mt) genome.fa
```
* illumina
    * Download Metadata from NCBI SRA Run Selector via a web browser
        * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA178613
```
cat SraRunTable.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f Experiment,"Sample\ Name",Bases \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    sed '1 s/^/#/' |
    keep-header -- tsv-sort -k2,2 -k3,3nr |
    tsv-uniq -H -f "Sample\ Name" --max 1 |
    mlr --itsv --ocsv cat \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

rgr md ena_info.tsv --fmt

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"
```
## Symlink
```
cd ~/data/plastid/Atha_cross/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/col0/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.tsv |
    tsv-select -H -f name,srr,bases |
    perl -nla -F'\t' -e '
        BEGIN { our %seen }
        /^name/ and next;
        $seen{$F[0]} and next;
        my $bases = $F[2];
        $bases =~ s/G$//;
        my $cutoff = $bases * 1000 * 1000 * 1000 / $ENV{GENOME_SIZE} * $ENV{FOLD};
        $cutoff = int $cutoff;
        print join qq(\t), ($F[0], $F[1], $cutoff);
        $seen{$F[0]}++;
    ' \
    > opts.tsv

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ ! -f ena/{2}_1.fastq.gz ]; then
            exit;
        fi
        if [ ! -f ena/{2}_2.fastq.gz ]; then
            exit;
        fi

        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        mkdir -p {1}/1_genome
        pushd {1}/1_genome

        cp ../../genome/genome.fa genome.fa
        popd > /dev/null

        mkdir -p {1}/2_illumina
        pushd {1}/2_illumina

        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd > /dev/null
    '
```
## Run anchr
```
cd Atha_cross
cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        if [ ! -d {1} ]; then
            exit;
        fi

        if [ ! -e {1}/2_illumina/R1.fq.gz ]; then
            exit;
        fi
        if [ ! -e {1}/2_illumina/R2.fq.gz ]; then
            exit;
        fi

        if bjobs -w | tr -s " " | cut -d " " -f 7 | grep -w "^{1}$"; then
            echo Job {1} exists
            exit;
        fi

        cd {1}

        echo {1}

        rm *.sh
        anchr template \
            --genome 1000000 \
            --parallel 24 \
            --xmx 80g \
            \
            --fastqc \
            --insertsize \
            --fastk \
            \
            --trim "--dedupe --cutoff {3} --cutk 31" \
            --qual "25" \
            --len "60" \
            --filter "adapter artifact" \
            \
            --bwa Q25L60 \
            --gatk

        bsub -q mpi -n 24 -J "{1}" "
            bash 0_script/2_fastqc.sh
            bash 0_script/2_insert_size.sh
            bash 0_script/2_fastk.sh
            bash 0_script/2_trim.sh
            bash 0_script/9_stat_reads.sh
            bash 0_script/3_bwa.sh
            bash 0_script/3_gatk.sh
        "
    '
```
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
awk '$5 != $6' Col_Ler_GT.tsv > informative_sites.tsv
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
bcftools view -R <(cut -f1,2 informative_sites.tsv) ../Atha_cross.vcf.gz -Oz -o F2_informative.vcf.gz

bcftools index F2_informative.vcf.gz
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' F2_informative.vcf.gz > F2_GT_matrix.tsv
```
* calculate the proportion
```
awk -v OFS='\t' '
    BEGIN {
        while ((getline < "informative_sites.tsv") > 0) {
            key = $1":"$2
            col_gt[key] = $5
            ler_gt[key] = $6
        }

        count_col = count_ler = count_mut = total = 0
    }
    NR==1 {
        header = "CHROM_POS"
        for (i=3; i<=NF; i++) {
            header = header "\t" "SAMPLE" (i-2)
        }
        print header > "classification_matrix.tsv"
        next
    }
    {
        key = $1":"$2
        row = key

        for (i=3; i<=NF; i++) {
            gt = $i
            col = col_gt[key]
            ler = ler_gt[key]

            if (gt == "./.") {
                row = row "\tNA"
                continue
            }

            if (gt == col) {
                row = row "\tCol"
                count_col++
            } else if (gt == ler) {
                row = row "\tLer"
                count_ler++
            } else {
                row = row "\tMut"
                count_mut++
            }
            total++
        }

        print row >> "classification_matrix.tsv"
    }
    END {
        print "来源分类比例：" > "summary.txt"
        if (total > 0) {
            print "Col型\t" count_col "\t" count_col / total * 100 "%" >> "summary.txt"
            print "Ler型\t" count_ler "\t" count_ler / total * 100 "%" >> "summary.txt"
            print "突变型\t" count_mut "\t" count_mut / total * 100 "%" >> "summary.txt"
        } else {
            print "无有效基因型参与分类" >> "summary.txt"
        }
    }
' F2_GT_matrix.tsv
```
每个样本分开统计：
```
total_samples=$(cat ../opts.tsv | grep -Ev "^Sample_Col_G$|^Sample_Ler_XL_4$" | wc -l)
current=1

cat ../opts.tsv | grep -Ev "^Sample_Col_G$|^Sample_Ler_XL_4$" | while read -r line; do
    sample=$(echo "$line" | cut -f1)
    target_col=$((current + 2))  # 列号 = 样本序号 + 2 (跳过CHROM和POS列)
    
    echo "[$current/$total_samples] 正在处理样本 $sample，列号 $target_col"
    
    awk -v OFS='\t' -v target_col="$target_col" -v sample_name="$sample" '
        BEGIN {
            while ((getline < "informative_sites.tsv") > 0) {
                key = $1":"$2
                col_gt[key] = $5
                ler_gt[key] = $6
            }

            count_col = count_ler = count_mut = total = 0
            print "CHROM_POS", sample_name > "classification_matrix_" sample_name ".tsv"
        }

        NR > 1 {
            key = $1":"$2
            gt = $target_col
            col = col_gt[key]
            ler = ler_gt[key]
            
            if (gt == "./.") {
                label = "NA"
                print key, label >> "classification_matrix_" sample_name ".tsv"
                next
            }
            
            if (gt == col) {
                label = "Col"
                count_col++
            } else if (gt == ler) {
                label = "Ler"
                count_ler++
            } else {
                label = "Mut"
                count_mut++
            }
            
            total++
            print key, label >> "classification_matrix_" sample_name ".tsv"
        }

        END {
            print "来源分类比例：" > "summary_" sample_name ".txt"
            if (total > 0) {
                print "Col型\t" count_col "\t" count_col / total * 100 "%" >> "summary_" sample_name ".txt"
                print "Ler型\t" count_ler "\t" count_ler / total * 100 "%" >> "summary_" sample_name ".txt"
                print "突变型\t" count_mut "\t" count_mut / total * 100 "%" >> "summary_" sample_name ".txt"
            } else {
                print "无有效基因型参与分类" >> "summary_" sample_name ".txt"
            }
        }
    ' F2_GT_matrix.tsv
    
    current=$((current + 1))
done

echo -e "样本名称\tCol型比例(%)\tLer型比例(%)\t突变型比例(%)" > sample_ratio_summary.tsv
samples=$(cat ../opts.tsv | cut -f1 | grep -Ev "^Sample_Col_G$|^Sample_Ler_XL_4$")

for sample in $samples; do
    if [ -f "summary_${sample}.txt" ]; then
        col_ratio=$(grep "Col型" "summary_${sample}.txt" | awk '{print $3}')
        ler_ratio=$(grep "Ler型" "summary_${sample}.txt" | awk '{print $3}')
        mut_ratio=$(grep "突变型" "summary_${sample}.txt" | awk '{print $3}')
        
        echo -e "$sample\t$col_ratio\t$ler_ratio\t$mut_ratio" >> sample_ratio_summary.tsv
    else
        echo "警告: 样本 $sample 的统计文件不存在，跳过"
    fi
done

echo -e "| 样本名称 | Col型比例(%) | Ler型比例(%) | 突变型比例(%) |" > sample_ratio_summary.md
echo -e "| --- | --- | --- | --- |" >> sample_ratio_summary.md

tail -n +2 sample_ratio_summary.tsv | awk -F'\t' '{
    printf("| %s | %s | %s | %s |\n", $1, $2, $3, $4);
}' >> sample_ratio_summary.md
```
| 样本名称 | Col型比例(%) | Ler型比例(%) | 突变型比例(%) |
| --- | --- | --- | --- |
| Sample_14 | 100% | 0% | 0% |
| Sample_18 | 100% | 0% | 0% |
| Sample_19 | 50% | 0% | 50% |
| Sample_20 | 100% | 0% | 0% |
| Sample_21 | 66.6667% | 0% | 33.3333% |
| Sample_4 | 100% | 0% | 0% |
| Sample_5 | 100% | 0% | 0% |
| Sample_6 | 0% | 100% | 0% |
| Sample_7 | 75% | 0% | 25% |
| Sample_8 | 100% | 0% | 0% |
| Sample_c1c2 | 66.6667% | 0% | 33.3333% |
| Sample_c41 | 80% | 0% | 20% |
| Sample_c42 | 100% | 0% | 0% |
| Sample_c45 | 80% | 0% | 20% |
| Sample_c47 | 100% | 0% | 0% |
| Sample_c48 | 100% | 0% | 0% |
| Sample_c51 | 100% | 0% | 0% |
| Sample_c52 | 100% | 0% | 0% |
| Sample_c54 | 80% | 0% | 20% |
| Sample_c57 | 100% | 0% | 0% |
| Sample_c61 | 100% | 0% | 0% |
| Sample_c62 | 80% | 0% | 20% |
| Sample_c63 | 0% | 50% | 50% |
| Sample_c64 | 100% | 0% | 0% |
| Sample_c65 | 75% | 0% | 25% |
| Sample_c66 | 100% | 0% | 0% |
| Sample_c73 | 100% | 0% | 0% |
| Sample_c81 | 50% | 0% | 50% |
| Sample_c82 | 100% | 0% | 0% |
| Sample_c83 | 100% | 0% | 0% |
| Sample_c84 | 80% | 0% | 20% |
| Sample_c85 | 100% | 0% | 0% |
| Sample_c87 | 0% | 0% | 100% |
| Sample_c88 | 100% | 0% | 0% |
| Sample_c89 | 100% | 0% | 0% |
| Sample_c90 | 100% | 0% | 0% |
| Sample_c91 | 100% | 0% | 0% |
| Sample_c92 | 100% | 0% | 0% |
| Sample_c93 | 100% | 0% | 0% |
| Sample_c94 | 4% | 88% | 8% |
| Sample_c95 | 5.88235% | 80.3922% | 13.7255% |
| Sample_l2c2 | 2.04082% | 89.7959% | 8.16327% |
| Sample_l2l3 | 7.54717% | 83.0189% | 9.43396% |
| Sample_l4c1 | 1.69492% | 0% | 98.3051% |
| Sample_l4l3 | 1.69492% | 0% | 98.3051% |
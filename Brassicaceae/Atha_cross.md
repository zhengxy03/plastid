# Arabidopsis thaliana Columbia X Ler Genome sequencing
## Basic info
* Genome: GCF_000001735.3, TAIR10, 119.668 Mb<br>
* Chloroplast: NC_000932, Columbia, 154478 bp<br>
* Chloroplast: KX551970, Landsberg erecta, 154515 bp<br> Ler-0，另一种常用野生型，与 Col-0 存在遗传差异
* Mitochondrion: Y08501, 366924 bp<br>
## Project
* https://www.pnas.org/content/109/51/20992.long<br>
* PRJNA178613<br>
* https://www.nature.com/articles/nature14649<br>
* PRJNA243018<br>
* PRJNA232554 - rice<br>
* PRJNA252997 - bee<br>
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
        * https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA243018
```
cat SraRunTable.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f Experiment,"Sample\ Name",Bases \
    > SraRunTable.tsv

cat SraRunTable.2.csv |
    mlr --icsv --otsv cat |
    tsv-filter -H --str-eq Organism:"Arabidopsis thaliana" |
    tsv-select -H -f Experiment,"Sample\ Name",Bases |
    sed '1d' \
    >> SraRunTable.tsv

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
## Run
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
## Pack and clean
```
cd Atha_cross
cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ -f {1}.tar.gz ]; then
            echo "==> {1} .tar.gz"
            exit;
        fi

        if [ ! -f {1}/3_gatk/R.filtered.vcf ]; then
            echo "==> {1} 3_gatk"
            exit;
        fi

#        if [ ! -f {1}/7_merge_anchors/anchor.merge.fasta ]; then
#            echo "==> {1} 7_merge_anchors"
#            exit;
#        fi
#
#        if [ ! -d {1}/9_quast ]; then
#            echo "==> {1} 9_quast"
#            exit;
#        fi

        echo "==> Clean {1}"
        bash {1}/0_script/0_cleanup.sh

        echo "==> Create {1}.tar.gz"

        tar -czvf {1}.tar.gz \
            {1}/1_genome/genome.fa \
            {1}/2_illumina/fastqc \
            {1}/2_illumina/insert_size \
            {1}/2_illumina/fastk \
            {1}/3_bwa/join.tsv \
            {1}/3_bwa/R.dedup.metrics \
            {1}/3_bwa/R.wgs.metrics \
            {1}/3_bwa/R.sort.bam \
            {1}/3_bwa/R.sort.bai \
            {1}/3_gatk \
            {1}/7_merge_anchors/anchor.merge.fasta \
            {1}/7_merge_anchors/others.non-contained.fasta \
            {1}/8_megahit/anchor/anchor.fasta \
            {1}/8_megahit/megahit.non-contained.fasta \
            {1}/8_mr_megahit/anchor/anchor.fasta \
            {1}/8_mr_megahit/megahit.non-contained.fasta \
            {1}/8_spades/anchor/anchor.fasta \
            {1}/8_spades/spades.non-contained.fasta \
            {1}/8_mr_spades/anchor/anchor.fasta \
            {1}/8_mr_spades/spades.non-contained.fasta \
            {1}/8_platanus/anchor/anchor.fasta \
            {1}/8_platanus/platanus.non-contained.fasta \
            {1}/9_quast \
            {1}/*.md

        echo
    '

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            if [ -d {1} ]; then
                echo "==> Remove {1}/"
                rm -fr {1}
            fi
        fi
    '
```
* Remove processed files
```
cd Atha_cross

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ ! -f {1}.tar.gz ]; then
            exit;
        fi

        find ena -type f -name "*{2}*"
    ' |
    xargs rm
```
* Unpack
```
cd Atha_cross

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -d {1} ]; then
            echo "==> {1} exists"
            exit;
        fi

        if [ ! -f {1}.tar.gz ]; then
            echo "==> {1}.tar.gz not exists"
            exit;
        fi

        tar -xzvf {1}.tar.gz
        rm {1}.tar.gz
    '
```
## VCF
```
cd Atha_cross
cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -f {1}.tar.gz ]; then
            echo "==> {1}.tar.gz not exists"
            exit;
        fi

        tar -xOzvf {1}.tar.gz {1}/3_gatk/R.filtered.vcf |
            bcftools reheader --samples <(echo {1}) |
            bcftools view \
                --apply-filters PASS --types snps --max-alleles 2 --targets Pt -Oz |
            bcftools view --include "AF>0.01" -Oz -o vcf/{1}.vcf.gz

        bcftools index -f vcf/{1}.vcf.gz
    '

bcftools merge --merge all -l <(
        cat opts.tsv |
            cut -f 1 |
            parallel -k -j 1 ' [ -f vcf/{}.vcf.gz ] && echo "vcf/{}.vcf.gz" '
    ) \
    > Atha_cross.vcf

rm -fr vcf
```
## stat analysis
```
mkdir -p analysis
cd analysis

#提取亲本
cp ../Sample_Col_G/3_gatk/R.filtered.vcf Col.vcf
cp ../Sample_Ler_XL_4/3_gatk/R.filtered.vcf Ler.vcf

bgzip Col.vcf Ler.vcf

bcftools view --apply-filters PASS --max-alleles 2 --targets Pt  --include "AF>0.01" -Oz Col.vcf.gz > Col_Pt.vcf.gz
bcftools view --apply-filters PASS --max-alleles 2 --targets Pt  --include "AF>0.01" -Oz Ler.vcf.gz > Ler_Pt.vcf.gz



#后代
cd ..
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
bcftools index Col_Pt.vcf.gz
bcftools index Ler_Pt.vcf.gz
bgzip ../Atha_cross.vcf
bcftools index -f ../Atha_cross.vcf.gz


#亲本差异位点
bcftools reheader -s <(echo "Col") Col_Pt.vcf.gz -o Col_Pt.renamed.vcf.gz
bcftools reheader -s <(echo "Ler") Ler_Pt.vcf.gz -o Ler_Pt.renamed.vcf.gz
bcftools index Col_Pt.renamed.vcf.gz
bcftools index Ler_Pt.renamed.vcf.gz

bcftools merge -m all Col_Pt.renamed.vcf.gz Ler_Pt.renamed.vcf.gz -Oz -o Col_Ler.merged.vcf.gz
bcftools index Col_Ler.merged.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' Col_Ler.merged.vcf.gz > Col_Ler_GT.tsv
awk '$5 != $6' Col_Ler_GT.tsv > informative_sites.tsv

#提取F2 GT
bcftools view -R <(cut -f1,2 informative_sites.tsv) ../Atha_cross.vcf.gz -Oz -o F2_informative.vcf.gz

bcftools index F2_informative.vcf.gz
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' F2_informative.vcf.gz > F2_GT_matrix.tsv



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

#单个样本
i=14  # 第14个样本

sample=$(awk -v line=$i 'NR==line { print $1 }' ../opts.tsv)
target_col=$((i + 2))  

echo "正在处理样本 $sample，列号 $target_col"

awk -v OFS='\t' -v target_col="$target_col" -v sample_name="$sample" '
    BEGIN {
        while ((getline < "informative_sites.tsv") > 0) {
            key = $1":"$2
            col_gt[key] = $5
            ler_gt[key] = $6
        }

        count_col = count_ler = count_mut = total = 0
        print "CHROM_POS", sample_name > "classification_matrix_sample.tsv"
    }

    NR > 1 {
        key = $1":"$2
        gt = $target_col
        col = col_gt[key]
        ler = ler_gt[key]

        if (gt == "./.") {
            label = "NA"
            print key, label >> "classification_matrix_sample.tsv"
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
        print key, label >> "classification_matrix_sample.tsv"
    }

    END {
        print "来源分类比例：" > "summary_sample.txt"
        if (total > 0) {
            print "Col型\t" count_col "\t" count_col / total * 100 "%" >> "summary_sample.txt"
            print "Ler型\t" count_ler "\t" count_ler / total * 100 "%" >> "summary_sample.txt"
            print "突变型\t" count_mut "\t" count_mut / total * 100 "%" >> "summary_sample.txt"
        } else {
            print "无有效基因型参与分类" >> "summary_sample.txt"
        }
    }
' F2_GT_matrix.tsv


#汇总
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

```
## rdp4 analysis
```
bcftools consensus -f c44/1_genome/genome.fa c/44/3_gatk/R.filtered.vcf.gz > c44_chrC.fa


cat opts.tsv | cut -f 1 > sample_list

output_dir="consensus_fasta"
mkdir -p ${output_dir}

for sample in $(cat sample_list); do
    echo "[INFO] Processing $sample..."

    bgzip -f "${sample}/3_gatk/R.filtered.vcf"
    vcf="${sample}/3_gatk/R.filtered.vcf.gz"
    faops filter -l 0 ${sample}/1_genome/genome.fa stdout | grep -A 1 '^>Pt' | faops one -l 0 ${sample}/1_genome/genome.fa Pt ${sample}/1_genome/Pt.fa
    ref="${sample}/1_genome/Pt.fa"
    out="consensus_fasta/${sample}_chrC.fa"

    if [[ -f "$vcf" && -f "$ref" ]]; then
        bcftools index -f "$vcf"

        # 过滤出只有 Pt 染色体的VCF
        bcftools view \
            --apply-filters PASS \
            --max-alleles 2 \
            --targets Pt \
            --include "AF>0.01" \
            -Oz -o "${sample}/3_gatk/R.filtered.Pt.vcf.gz" "$vcf"
        bcftools index -f "${sample}/3_gatk/R.filtered.Pt.vcf.gz"

        bcftools consensus -f "$ref" "${sample}/3_gatk/R.filtered.Pt.vcf.gz" > "$out"
        echo "[OK] Wrote $out"
    else
        echo "[WARN] Missing files for $sample"
    fi
done

cat consensus_fasta/*_chrC.fa > all_samples_raw.fasta
mafft --auto all_samples_raw.fasta > all_samples_aligned.fasta
```
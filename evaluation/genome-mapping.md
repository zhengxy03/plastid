# Tools
```
curl -L https://delivery04-mul.dhe.ibm.com/sar/CMA/OSA/0cz8y/0/ibm-aspera-connect_4.2.14.855-HEAD_linux_x86_64.tar.gz | tar xvz
bash ibm-aspera-connect*

cat >> $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh <<EOF
-----BEGIN DSA PRIVATE KEY-----
MIIBuwIBAAKBgQDkKQHD6m4yIxgjsey6Pny46acZXERsJHy54p/BqXIyYkVOAkEp
KgvT3qTTNmykWWw4ovOP1+Di1c/2FpYcllcTphkWcS8lA7j012mUEecXavXjPPG0
i3t5vtB8xLy33kQ3e9v9/Lwh0xcRfua0d5UfFwopBIAXvJAr3B6raps8+QIVALws
yeqsx3EolCaCVXJf+61ceJppAoGAPoPtEP4yzHG2XtcxCfXab4u9zE6wPz4ePJt0
UTn3fUvnQmJT7i0KVCRr3g2H2OZMWF12y0jUq8QBuZ2so3CHee7W1VmAdbN7Fxc+
cyV9nE6zURqAaPyt2bE+rgM1pP6LQUYxgD3xKdv1ZG+kDIDEf6U3onjcKbmA6ckx
T6GavoACgYEAobapDv5p2foH+cG5K07sIFD9r0RD7uKJnlqjYAXzFc8U76wXKgu6
WXup2ac0Co+RnZp7Hsa9G+E+iJ6poI9pOR08XTdPly4yDULNST4PwlfrbSFT9FVh
zkWfpOvAUc8fkQAhZqv/PE6VhFQ8w03Z8GpqXx7b3NvBR+EfIx368KoCFEyfl0vH
Ta7g6mGwIMXrdTQQ8fZs
-----END DSA PRIVATE KEY-----
EOF
chmod 600 $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh

cd $HOME/bin
ln -s ../.aspera/connect/bin/ascp ascp
```
# Reference genomes
* Arabidopsis thaliana Col-0
```
mkdir -p plastid/genome/col0
cd plastid/genome/col0

wget -N https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
#dna_sm:将重复序列和低复杂度区域的碱基转换为小写形式（如 a、t、c、g），而功能基因等非重复区域仍保持大写格式，这是一种 “软屏蔽” 处理方式。这种处理既标记出了重复区域，又保留了重复区域的碱基信息

hnsm order Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz \
    <(for chr in {1,2,3,4,5,Mt,Pt}; do echo $chr; done) \
    -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes
```
* Oryza sativa Japonica Group Cultivar Nipponbare
日本晴
```
mkdir nip && cd nip 

wget -N https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz

hnsm order Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz \
    <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) \
    -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes
```
* Medicago truncatula A17
蒺藜苜蓿<br>

The Ensembl version lacks chromosome naming and does not include chloroplast and mitochondrial genomes, which need to be constructed manually.
```
mkdir a17 && cd a17

aria2c -x 4 -s 2 -c \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/473/485/GCF_003473485.1_MtrunA17r5.0-ANR/GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna.gz

gzip -dc GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna.gz | head -n 2
TAB=$'\t'
cat <<EOF > replace.tsv
NC_053042.1${TAB}1
NC_053043.1${TAB}2
NC_053044.1${TAB}3
NC_053045.1${TAB}4
NC_053046.1${TAB}5
NC_053047.1${TAB}6
NC_053048.1${TAB}7
NC_053049.1${TAB}8
NC_003119.8${TAB}Pt
NC_029641.1${TAB}Mt
EOF

gzip -dcf GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna.gz |
    hnsm replace stdin replace.tsv |
    hnsm order stdin <(for chr in $(seq 1 1 8) Mt Pt; do echo $chr; done) -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes
```
* Solanum lycopersicum Micro-Tom
微型番茄品种，植株矮小紧凑，存在大量单核苷酸多态性（SNPs）和插入 / 缺失（Indels）<br>
```
mkdir microtom && cd microtom

aria2c -x 4 -s 2 -c \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/512/215/GCF_036512215.1_SLM_r2.1/GCF_036512215.1_SLM_r2.1_genomic.fna.gz

TAB=$'\t'
cat <<EOF > replace.tsv
NC_090800.1${TAB}1
NC_090801.1${TAB}2
NC_090802.1${TAB}3
NC_090803.1${TAB}4
NC_090804.1${TAB}5
NC_090805.1${TAB}6
NC_090806.1${TAB}7
NC_090807.1${TAB}8
NC_090808.1${TAB}9
NC_090809.1${TAB}10
NC_090810.1${TAB}11
NC_090811.1${TAB}12
NC_035963.1${TAB}Mt
NC_007898.3${TAB}Pt
EOF

gzip -dcf GCF_036512215*_genomic.fna.gz |
    hnsm replace stdin replace.tsv |
    hnsm order stdin <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes
```
* Solanum lycopersicum Cultivar: Heinz 1706
普通番茄品种<br>
```
mkdir h1706
cd h1706

aria2c -x 4 -s 2 -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.5_SL3.1/GCF_000188115.5_SL3.1_genomic.fna.gz

AB=$'\t'
cat <<EOF > replace.tsv
NC_015438.3${TAB}1
NC_015439.3${TAB}2
NC_015440.3${TAB}3
NC_015441.3${TAB}4
NC_015442.3${TAB}5
NC_015443.3${TAB}6
NC_015444.3${TAB}7
NC_015445.3${TAB}8
NC_015446.3${TAB}9
NC_015447.3${TAB}10
NC_015448.3${TAB}11
NC_015449.3${TAB}12
NC_035963.1${TAB}Mt
NC_007898.3${TAB}Pt
EOF

gzip -dcf GCF*_genomic.fna.gz |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) genome.fa

faops size genome.fa > chr.sizes
```

* Prunus persica PLov2-2N (a double haploid genotype of the peach cv. Lovell)
桃的一种双单倍体<br>
```
mkdir lovell
cd lovell

aria2c -x 4 -s 2 -c \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz

TAB=$'\t'
cat <<EOF > replace.tsv
NC_034009.1${TAB}G1
NC_034010.1${TAB}G2
NC_034011.1${TAB}G3
NC_034012.1${TAB}G4
NC_034013.1${TAB}G5
NC_034014.1${TAB}G6
NC_034015.1${TAB}G7
NC_034016.1${TAB}G8
NC_014697.1${TAB}Pt
EOF

gzip -dcf GCF_000346465*_genomic.fna.gz |
    hnsm replace -s stdin replace.tsv -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes
```
* Glycine max Williams 82
大豆
```
mkdir w82
cd w82

aria2c -x 4 -s 2 -c \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz

TAB=$'\t'
cat <<EOF > replace.tsv
NC_016088.3${TAB}1
NC_016089.3${TAB}2
NC_016090.3${TAB}3
NC_016091.3${TAB}4
NC_038241.1${TAB}5
NC_038242.1${TAB}6
NC_038243.1${TAB}7
NC_038244.1${TAB}8
NC_038245.1${TAB}9
NC_038246.1${TAB}10
NC_038247.1${TAB}11
NC_038248.1${TAB}12
NC_038249.1${TAB}13
NC_038250.1${TAB}14
NC_038251.1${TAB}15
NC_038252.1${TAB}16
NC_038253.1${TAB}17
NC_038254.1${TAB}18
NC_038255.1${TAB}19
NC_038256.1${TAB}20
NC_007942.1${TAB}Pt
NC_020455.1${TAB}Mt
EOF

gzip -dcf GCF*_genomic.fna.gz |
    hnsm replace stdin replace.tsv |
    hnsm order stdin <(for chr in $(seq 1 1 20) Pt Mt; do echo $chr; done) -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes
```
# Download fastq files from ENA
```
mkdir -p ~/data/plastid/ena
cd ~/data/plastid/ena

cat << EOF > source.csv
SRX202246,Atha_Col_0_1,HiSeq 2000 PE100
SRX2527206,Atha_Col_0_2,MiSeq 2000 PE300
SRR616965,Atha_Ler_0,HiSeq 2000 PE100
SRX179254,Osat_Nip,HiSeq 2000 PE100
SRX025260,Osat_Nip_2,Osat_50
SRX673852,Mtru_A17,HiSeq 2000 PE150
SRX150254,Pper_Lovell,Illumina Genome Analyzer IIx PE100
SRX698770,Slyc_H1706,Illumina HiSeq 2000 PE100
SRX7009428,Gmax_W82,HiSeq X Ten
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

rgr md ena_info.tsv --fmt

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"

#aria2c -j 4 -x 4 -s 1 -c -i ena_info.ftp.txt
md5sum --check ena_info.md5.txt

# Bad quality
# SRR1542422
```
# Basic info
* `cutoff = FOLD * DEPTH` <br>
FOLD 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64<br>
* 叶绿体的拷贝数>线粒体>核基因组
# Symlink
```
mkdir -p evaluation
cd evaluation

SRRS=(
    'SRR616966::col0'  # Col-0
    'SRR611086::col0'
    'SRR5216995::col0'
    'SRR616965::col0'  # Ler-0
    'SRR611087::col0'
    'SRR545231::nip'    # Nipponbare
    'SRR063638::nip'
    'SRR1542423::a17'   # A17
    'SRR1572628::h1706' # Heinz 1706
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    GENOME="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        mkdir -p ${BASE_NAME}/1_genome
        pushd ${BASE_NAME}/1_genome

        ln -fs ../../../genome/${GENOME}/genome.fa genome.fa
        cp ../../../genome/${GENOME}/chr.sizes chr.sizes
        popd

        mkdir -p ${BASE_NAME}/2_illumina
        pushd ${BASE_NAME}/2_illumina

        ln -fs ../../../ena/${SRR}_1.fastq.gz R1.fq.gz
        ln -fs ../../../ena/${SRR}_2.fastq.gz R2.fq.gz
        popd
    done
done
```
# Trim, cutoff and mapping
* Rsync to hpcc
```
rsync -avP \
    ~/data/plastid/ \
    wangq@202.119.37.251:data/plastid
```
```
cd evaluation
SRRS=(
    'SRR616966::41'  # Col-0
    'SRR611086::78'
    'SRR5216995::121'
    'SRR616965::42'  # Ler-0
    'SRR611087::79'
    'SRR545231::46'  # Nipponbare
    'SRR063638::15'
    'SRR1542423::23' # A17
    'SRR1572628::5'  # Heinz 1706
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    DEPTH="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")

        echo 1>&2 "==> ${item} ${FOLD}"

        BASE_NAME=${SRR}_${FOLD}
        pushd ${BASE_NAME}

        if [ ! -f 3_bwa/join.tsv ]; then
            rm *.sh
            if [[ "${FOLD}" == "0" ]]; then
                anchr template \
                    --genome $(tsv-summarize 1_genome/chr.sizes --sum 2) \
                    --parallel 24 \
                    --xmx 80g \
                    \
                    --fastqc \
                    --insertsize \
                    --fastk \
                    \
                    --trim "--dedupe" \
                    --qual "25" \
                    --len "60" \
                    --filter "adapter artifact" \
                    \
                    --bwa Q25L60

                bsub -q mpi -n 24 -J "${BASE_NAME}" "
                    bash 0_script/2_fastqc.sh
                    bash 0_script/2_insert_size.sh
                    bash 0_script/_kat.sh
                    bash 0_script/2_trim.sh
                    bash 0_script/9_stat_reads.sh
                    bash 0_script/3_bwa.sh
                "
            else
                anchr template \
                    --genome $(tsv-summarize 1_genome/chr.sizes --sum 2) \
                    --parallel 24 \
                    --xmx 80g \
                    \
                    --trim "--dedupe --cutoff ${CUTOFF} --cutk 31" \
                    --qual "25" \
                    --len "60" \
                    --filter "adapter artifact" \
                    \
                    --bwa Q25L60

                bsub -q mpi -n 24 -J "${BASE_NAME}" "
                    bash 0_script/2_trim.sh
                    bash 0_script/9_stat_reads.sh
                    bash 0_script/3_bwa.sh
                "
            fi

        fi

        popd

    done
done
```
# Combine chromosomes
```
cd evaluation

SRRS=(
    'SRR616966::5'   # Col-0
    'SRR611086::5'
    'SRR5216995::5'
    'SRR616965::5'   # Ler-0
    'SRR611087::5'
    'SRR545231::12'  # Nipponbare
    'SRR063638::12'
    'SRR1542423::8'  # A17
    'SRR1572628::12' # Heinz 1706
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    CHR_NUM="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        echo 1>&2 "==> ${BASE_NAME}"

        pushd ${BASE_NAME}/3_bwa

        cat join.tsv |
            grep -v "^Mt" |
            grep -v "^Pt" |
            tsv-summarize -H --sum chrLength,covLength,bases --min min --max max |
            sed '1d' |
            perl -e '
                my $line = <>;
                chomp $line;
                my ($chrLength, $covLength, $bases, $min, $max, ) = split qq(\t), $line;
                my $covRate = sprintf qq(%.4f), $covLength / $chrLength;
                my $mean = sprintf qq(%.2f), $bases / $chrLength;
                print join qq(\t), (
                    "Nc", $chrLength, $covLength, $covRate, $bases, $mean, $min, $max,
                );
                print qq(\n);
            ' |
            (cat join.tsv | sed "2,$((CHR_NUM+1))d" && cat) \
            > combine.tsv

        popd

    done
done
```
# Merge all results
```
cd evaluation

SRRS=(
    'SRR616966::Col-0'      # Col-0
    'SRR611086::Col-0'
    'SRR5216995::Col-0'
    'SRR616965::Ler-0'      # Ler-0
    'SRR611087::Ler-0'
    'SRR545231::Nipponbare' # Nipponbare
    'SRR063638::NP'
    'SRR1542423::A17'       # A17
    'SRR1572628::Heinz1706' # Heinz 1706
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    STRAIN="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        pushd ${BASE_NAME}/3_bwa > /dev/null

        echo -e "Fold\tchrom\n${FOLD}\tNc\n${FOLD}\tMt\n${FOLD}\tPt" |
            tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

        popd > /dev/null

    done |
        tsv-uniq \
        > ${SRR}_folds.tsv

    echo
    echo "Table: ${STRAIN} ${SRR} Folds"
    echo
    for PART in Nc Mt Pt; do
        cat ${SRR}_folds.tsv |
            tsv-filter -H --str-eq chrom:${PART} |
            mlr --itsv --omd cat
        echo
        echo
    done

done
```
可视化：
```R
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(stringr)
library(scales)

# 样本顺序
ordered_samples <- c(
  'SRR616966::Col-0',
  'SRR611086::Col-0',
  'SRR5216995::Col-0',
  'SRR616965::Ler-0',
  'SRR611087::Ler-0',
  'SRR545231::Nipponbare',
  'SRR063638::NP',
  'SRR1542423::A17',
  'SRR1572628::Heinz1706'
)

sample_mapping <- tibble(
  SRR = str_extract(ordered_samples, "^SRR\\d+"),
  Label = ordered_samples
)

# 颜色和染色体顺序
colors <- c("Nc" = "black", "Pt" = "gray60", "Mt" = "white")
chrom_order <- c("Nc", "Pt", "Mt")
bar_width <- 0.4
dodge_width <- 0.35

# 读取所有 _folds.tsv 文件
all_tsv <- list.files(pattern = "_folds\\.tsv$")
all_sample_plots <- list()

for (srr in sample_mapping$SRR) {
  file <- all_tsv[grepl(srr, all_tsv)]
  if (length(file) == 0) {
    message("WARNING: No file found for ", srr)
    next
  }
  label <- sample_mapping$Label[sample_mapping$SRR == srr]
  
  df <- read_tsv(file, show_col_types = FALSE)
  
  # 强制染色体顺序
  df$chrom <- factor(df$chrom, levels = chrom_order)
  
  # 图1：Mapped bases
  p1 <- ggplot(df, aes(x = factor(Fold), y = bases / 1e6, fill = chrom, group = chrom)) +
    geom_bar(stat = "identity", position = position_dodge(width = dodge_width),
             color = "black", width = bar_width) +
    scale_fill_manual(values = colors, breaks = chrom_order, name = "Genome") +
    labs(
      title = label,
      x = "Fold",
      y = "Mapped bases (Million)"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 14)
    )
  
  # 图3：Proportion of genome，直接用 covRate
  p3 <- ggplot(df, aes(x = factor(Fold), y = covRate, fill = chrom, group = chrom)) +
    geom_bar(stat = "identity", position = position_dodge(width = dodge_width),
             color = "black", width = bar_width) +
    scale_fill_manual(values = colors, breaks = chrom_order, name = "Genome") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
    labs(
      title = label,
      x = "Fold",
      y = "Proportion of genome"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 14)
    )
  
  # 组合两张图，共用图例
  row_plot <- (p1 | p3) + plot_layout(guides = "collect")
  all_sample_plots[[label]] <- row_plot
}

# 按指定顺序排列样本
plots_in_order <- all_sample_plots[ordered_samples[ordered_samples %in% names(all_sample_plots)]]

# 每个样本一行纵向堆叠
big_plot <- wrap_plots(plots_in_order, ncol = 1)

# 保存图片
ggsave("all_samples_combined.png", plot = big_plot,
       width = 16,
       height = 4 * length(plots_in_order),
       limitsize = FALSE)

```
![evaluation](../pic/all_samples_combined.png)

# gene annotation
```
cd ~/zxy/plastid/evaluation/SRR616966_1/1_genome

# 1️⃣ 生成染色体长度文件
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > genome.chrom.sizes



# 核基因组（染色体 1-5）
perl -ne 'chomp; @f=split(/\t/,$_); print "$f[0]\t0\t$f[1]\n" if $f[0]=~/^[1-5]$/;' genome.chrom.sizes > nuclear_regions.bed

# 线粒体
perl -ne 'chomp; @f=split(/\t/,$_); print "$f[0]\t0\t$f[1]\n" if $f[0] eq "Mt";' genome.chrom.sizes > mito_regions.bed

# 叶绿体（可选）
perl -ne 'chomp; @f=split(/\t/,$_); print "$f[0]\t0\t$f[1]\n" if $f[0] eq "Pt";' genome.chrom.sizes > chloroplast_regions.bed

# 3️⃣ 提取落在核/线粒体的残留 reads
cd ~/zxy/plastid/evaluation/SRR616966_1/3_bwa

samtools view -b -L ../1_genome/nuclear_regions.bed R.sort.bam > nuclear_reads.bam
samtools view -b -L ../1_genome/mito_regions.bed R.sort.bam > mito_reads.bam

# 4️⃣ BAM 转 BED
bedtools bamtobed -i nuclear_reads.bam > nuclear_reads.bed
bedtools bamtobed -i mito_reads.bam > mito_reads.bed

# 5️⃣ 可选：提取 FASTA 用于序列分析
bedtools getfasta -fi ../1_genome/genome.fa -bed nuclear_reads.bed -fo nuclear_reads.fa
bedtools getfasta -fi ../1_genome/genome.fa -bed mito_reads.bed -fo mito_reads.fa

# 核基因组残留 reads
samtools view -b -L ../1_genome/nuclear_regions.bed R.sort.bam > nuclear_reads.bam

# 线粒体残留 reads
samtools view -b -L ../1_genome/mito_regions.bed R.sort.bam > mito_reads.bam

# 查看数量
samtools view -c nuclear_reads.bam
samtools view -c mito_reads.bam

bedtools bamtobed -i nuclear_reads.bam > nuclear_reads.bed
bedtools bamtobed -i mito_reads.bam > mito_reads.bed

#统计残留 reads 的覆盖区域特征
# 1. 统计核基因组残留reads的染色体分布
perl -ne 'chomp; @f=split(/\t/); $count{$f[0]}++; END{print "$_\t$count{$_}\n" for sort keys %count}' nuclear_reads.bed

# 2. 计算核基因组残留reads的平均长度
perl -ne 'chomp; @f=split(/\t/); $len=$f[2]-$f[1]; $sum+=$len; $count++; END{printf "平均长度: %.2f\n", $sum/$count}' nuclear_reads.bed

# 3. 分析线粒体残留reads的覆盖密度（输出每1bp的覆盖度）
bedtools genomecov -i mito_reads.bed -g ~/zxy/plastid/evaluation/SRR616966_1/1_genome/genome.chrom.sizes -d | perl -ne 'chomp; @f=split(/\t/); print "$_\n" if $f[0] eq "Mt"' > mito_coverage.txt

#三、筛选多位置比对的 reads）
# 1. 提取核残留reads中存在多位置比对的reads ID
samtools view -f 256 nuclear_reads.bam | perl -ne 'chomp; @f=split(/\t/); print "$f[0]\n"' | sort | uniq > nuclear_multimap_reads.txt

# 2. 从FASTA中提取这些候选嵌合序列
samtools fasta nuclear_reads.bam > nuclear_reads_with_id.fa
perl -ne 'BEGIN {
    open my $fh, "<", "nuclear_multimap_reads.txt" or die "Can\"t open ID file: $!";
    while (<$fh>) {
        chomp;
        $ids{$_} = 1;
    }
    close $fh;
}
if (/^>(\S+)\/\d+$/) {  # 匹配并移除/1或/2后缀
    $keep = exists $ids{$1} ? 1 : 0;
}
print if $keep;' nuclear_reads_with_id.fa > nuclear_chimeric_candidates.fa

#处理重复序列（以转座子为主）
cd ~/zxy/plastid/evaluation/SRR616966_1/1_genome

# 处理TAIR10_Transposable_Elements.txt（根据实际列名调整）
# 列名对应：
# 1: Transposon_Name（转座子ID）
# 2: orientation_is_5prime（方向相关，暂不使用）
# 3: Transposon_min_Start（起始位置）
# 4: Transposon_max_End（终止位置）
# 5: Transposon_Family（家族）
# 6: Transposon_Super_Family（超家族，如LTR/Gypsy）

perl -ne 'chomp;
          next if /^Transposon_Name/;  # 跳过表头行
          @f = split(/\t/);
          # 提取关键信息（处理染色体ID：从转座子名称中提取，如AT1TE...→染色体1）
          $te_id = $f[0];
          if ($te_id =~ /^AT(\d)TE/) {  # 从ID中提取染色体号（AT1TE→1，AT2TE→2）
              $chr = $1;
          } else {
              next;  # 跳过非核染色体的转座子（如线粒体/叶绿体上的）
          }
          $start = $f[2] - 1;  # 转换为BED的0起始格式
          $end = $f[3];
          $family = $f[5];
          $super_family = $f[6];
          $type = "$super_family($family)";  # 合并类型信息
          # 仅保留1-5号染色体
          if ($chr =~ /^[1-5]$/) {
              print "$chr\t$start\t$end\t$type\t$te_id\n";
          }' TAIR10_Transposable_Elements.txt > nuclear_te.bed

# 验证转换结果（查看前5行）
head nuclear_te.bed

# 后续分析步骤与之前一致，直接使用nuclear_te.bed即可
# 1️⃣ 统计转座元件类型分布
cut -f4 nuclear_te.bed | cut -d'(' -f1 | sort | uniq -c | awk '{print $2 "\t" $1}' > te_superfamily_counts.txt
echo "转座元件超家族分布："
cut -f4 nuclear_te.bed | cut -d'(' -f1 | sort | uniq -c | awk '{print $2 "\t" $1}' > te_superfamily_counts.txt

# 2️⃣ 提取核残留reads中位于转座元件区域的reads
cd ~/zxy/plastid/evaluation/SRR616966_1/3_bwa

# 转座元件区域的reads
bedtools intersect -abam nuclear_reads.bam -b ../1_genome/nuclear_te.bed -u > nuclear_te_reads.bam

# 非转座元件区域的reads
bedtools intersect -abam nuclear_reads.bam -b ../1_genome/nuclear_te.bed -v > nuclear_non_te_reads.bam

# 3️⃣ 统计转座元件相关reads的基础数据
echo "核基因组总残留reads：$(samtools view -c nuclear_reads.bam)" > te_stats.txt
echo "转座元件区域reads：$(samtools view -c nuclear_te_reads.bam)" >> te_stats.txt
echo "转座元件区域reads占比：$(echo "scale=2; $(samtools view -c nuclear_te_reads.bam)/$(samtools view -c nuclear_reads.bam)*100" | bc)%" >> te_stats.txt

# 4️⃣ 分析多位置比对序列中来自转座元件的比例
bedtools bamtobed -i nuclear_reads.bam | grep -Ff nuclear_multimap_reads.txt > nuclear_multimap.bed
bedtools intersect -a nuclear_multimap.bed -b ../1_genome/nuclear_te.bed -wa -wb | cut -f4 | sort | uniq > te_multimap_ids.txt

echo "多位置比对序列总数：$(wc -l < nuclear_multimap_reads.txt)" >> te_stats.txt
echo "其中来自转座元件的多比对序列：$(wc -l < te_multimap_ids.txt)" >> te_stats.txt
echo "转座元件导致的多比对占比：$(echo "scale=2; $(wc -l < te_multimap_ids.txt)/$(wc -l < nuclear_multimap_reads.txt)*100" | bc)%" >> te_stats.txt

# 5️⃣ 不同超家族转座元件的reads覆盖深度
cd ~/zxy/plastid/evaluation/SRR616966_1/3_bwa

# 1️⃣ 重新提取超家族并替换名称中的斜杠（避免路径问题）
cut -f4 ../1_genome/nuclear_te.bed | sed 's/^(//; s/)$//' | sort | uniq | sed 's/\//_/g' > te_superfamilies_sanitized.txt

# 显示处理后的超家族名称（不含斜杠）
echo "处理后的超家族名称（已替换斜杠）："
cat te_superfamilies_sanitized.txt

# 2️⃣ 计算不同超家族的平均覆盖深度（修正版）
# 2️⃣ 重新计算超家族覆盖深度（基于覆盖碱基数）
> te_superfamily_depth_final.txt
echo -e "超家族\t转座元件数量\t总长度(bp)\t覆盖碱基数\t平均深度" > te_superfamily_depth_final.txt

while read sanitized_sf; do
    [ -z "$sanitized_sf" ] && continue
    
    # 还原原始超家族名称
    original_sf=$(echo "$sanitized_sf" | sed 's/_/\//g')
    
    # 生成该超家族的BED文件
    grep "($original_sf)" ../1_genome/nuclear_te.bed > ../1_genome/te_${sanitized_sf}.bed
    te_count=$(wc -l < ../1_genome/te_${sanitized_sf}.bed)
    
    if [ $te_count -gt 0 ]; then
        # 计算该超家族所有转座元件的总长度
        total_length=$(awk '{sum += $3 - $2} END {print sum}' ../1_genome/te_${sanitized_sf}.bed)
        
        # 计算被reads覆盖的总碱基数（-d输出每个碱基的覆盖度）
        covered_bases=$(bedtools coverage -a ../1_genome/te_${sanitized_sf}.bed \
                                         -b nuclear_te_reads.bam \
                                         -d | awk '{sum += $7} END {print sum}')
        
        # 计算平均深度（覆盖碱基数 / 总长度）
        if [ $total_length -gt 0 ] && [ $covered_bases -gt 0 ]; then
            avg_depth=$(echo "scale=4; $covered_bases / $total_length" | bc)
        else
            avg_depth=0
        fi
    else
        total_length=0
        covered_bases=0
        avg_depth=0
    fi
    
    # 保留2位小数输出
    printf "%s\t%d\t%d\t%d\t%.2f\n" \
           "$original_sf" $te_count $total_length $covered_bases $avg_depth >> te_superfamily_depth_final.txt
    
    # 清理临时文件
    rm -f ../1_genome/te_${sanitized_sf}.bed
done < te_superfamilies_sanitized.txt

# 3️⃣ 显示最终结果
echo -e "\n最终转座元件超家族深度分析结果："
cat te_superfamily_depth_final.txt




#
# ==============================================
# 第一步：修复BAM文件头信息（解决SQ lines缺失问题）
# ==============================================
cd ~/zxy/plastid/evaluation/SRR616966_1/3_bwa

# 1. 生成参考基因组索引
samtools faidx ../1_genome/genome.fa

# 2. 创建序列字典
awk 'BEGIN{print "@HD\tVN:1.5\tSO:coordinate"} {print "@SQ\tSN:"$1"\tLN:"$2}' ../1_genome/genome.fa.fai > genome.dict

# 3. 下载Picard工具（如已存在可跳过）
wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard-3.1.1.jar -O picard.jar

# 4. 用Picard重建BAM头信息
java -jar picard.jar AddOrReplaceReadGroups \
    I=nuclear_reads.bam \
    O=nuclear_reads_fixed.bam \
    R=../1_genome/genome.fa \
    ID=SRR616966 \
    LB=library1 \
    PL=illumina \
    PU=unit1 \
    SM=sample1

# 5. 为修复后的BAM创建索引
samtools index nuclear_reads_fixed.bam

# 6. 验证修复结果（必须显示@SQ开头的行）
samtools view -H nuclear_reads_fixed.bam | grep '@SQ' | head -3


# ==============================================
# 第二步：构建BLAST数据库
# ==============================================
cd ../1_genome
makeblastdb -in genome.fa -dbtype nucl -out combined_ref_db -parse_seqids
cd ../3_bwa


# ==============================================
# 第三步：筛选候选ID（排除转座元件）
# ==============================================
comm -23 <(sort nuclear_multimap_reads.txt | uniq | sed 's/[[:space:]]//g') \
         <(sort te_multimap_ids.txt | uniq | sed 's/[[:space:]]//g') \
         > non_te_multimap_ids_clean.txt

# 查看候选ID数量
echo "候选ID数量：$(wc -l < non_te_multimap_ids_clean.txt)"


# ==============================================
# 第四步：提取潜在嵌合体序列（兼容低版本samtools）
# ==============================================
# 1. 导出完整SAM文件（含头信息）
samtools view -h nuclear_reads_fixed.bam > full_reads.sam

# 2. 筛选目标ID的reads
awk 'NR==1 || /^@/; BEGIN{while((getline < "non_te_multimap_ids_clean.txt")>0)ids[$1]=1} $1 in ids' full_reads.sam > target_reads.sam

# 3. 转换为FASTA
samtools fasta target_reads.sam > non_te_multimap.fasta

# 4. 清理临时文件
rm -f full_reads.sam target_reads.sam

# 验证提取结果
echo "提取的潜在嵌合体序列数：$(grep -c '^>' non_te_multimap.fasta)"


# ==============================================
# 第五步：BLAST筛选真实嵌合体
# ==============================================
blastn -query non_te_multimap.fasta \
       -db ../1_genome/combined_ref_db \
       -evalue 1e-20 \
       -outfmt "6 qseqid sseqid pident length evalue" \
       -max_target_seqs 2 \
       -num_threads 4 > blast_temp.txt

perl -ne 'chomp; @f=split(/\t/); next if @f<5;
         $id=$f[0]; $chr=$f[1]; $p=$f[2]; $e=$f[4];
         
         # 基因组分类：核(nuc)/叶绿体(cp)/线粒体(mt)
         $g = ($chr=~/^[1-5]$/) ? "nuc" :
              ($chr eq "Pt") ? "cp" :
              ($chr eq "Mt") ? "mt" : "other";
         next if $g eq "other";  # 排除未知基因组
         
         # 保留符合基本质量的比对（可根据数据调整）
         if($p>=70 && $e<=1e-10) {
             push @{$h{$id}}, $g;  # 存储基因组类型
         }
         
         END{
             for my $id (keys %h) {
                 my %unique = map { $_ => 1 } @{$h{$id}};  # 去重基因组类型
                 my @types = keys %unique;
                 
                 # 筛选跨两种及以上基因组的嵌合体
                 if(@types >= 2) {
                     # 三种目标类型：核-叶绿体、核-线粒体、叶绿体-线粒体
                     if( ($types[0] eq "nuc" && $types[1] eq "cp") ||
                         ($types[0] eq "cp" && $types[1] eq "nuc") ) {
                         print "$id\n";  # 核-叶绿体
                     } elsif( ($types[0] eq "nuc" && $types[1] eq "mt") ||
                              ($types[0] eq "mt" && $types[1] eq "nuc") ) {
                         print "$id\n";  # 核-线粒体
                     } elsif( ($types[0] eq "cp" && $types[1] eq "mt") ||
                              ($types[0] eq "mt" && $types[1] eq "cp") ) {
                         print "$id\n";  # 叶绿体-线粒体
                     }
                 }
             }
         }' blast_temp.txt > candidate_chimeric_ids.txt

# 查看所有类型候选嵌合体数量
echo "所有类型嵌合体候选ID数量：$(wc -l < candidate_chimeric_ids.txt)"

# 1. 去除候选ID中的/1和/2后缀，与BAM中的ID格式统一
sed 's/\/[12]$//' candidate_chimeric_ids.txt > candidate_ids_fixed.txt

# 2. 重新提取序列（使用修复后的ID列表）
samtools view -h nuclear_reads_fixed.bam > full_reads2.sam

# 3. 筛选目标ID的reads（此时ID格式已匹配）
awk 'NR==1 || /^@/; BEGIN{while((getline < "candidate_ids_fixed.txt")>0)ids[$1]=1} $1 in ids' full_reads2.sam > target_reads2.sam

# 4. 转换为FASTA
samtools fasta target_reads2.sam > true_chimeric_seqs.fasta

# 5. 清理临时文件
rm -f full_reads2.sam target_reads2.sam candidate_ids_fixed.txt

# 验证结果
echo "修复ID格式后提取的序列数：$(grep -c '^>' true_chimeric_seqs.fasta)"


# 1. BLAST验证嵌合体序列（确保比对准确性）
blastn -query true_chimeric_seqs.fasta \
       -db ../1_genome/combined_ref_db \
       -evalue 1e-10 \
       -outfmt "6 qseqid sseqid length pident evalue" \
       -max_target_seqs 2 \
       -num_threads 4 > chimeric_blast_verify.txt

# 2. 统计三种嵌合体类型（核-叶绿体/核-线粒体/叶绿体-线粒体）及特征
perl -ne 'chomp; @f=split(/\t/); next if @f<5;
         $id=$f[0]; $chr=$f[1]; $len=$f[2]; $p=$f[3];
         
         # 基因组分类
         $g = ($chr=~/^[1-5]$/) ? "nuc" :
              ($chr eq "Pt") ? "cp" :
              ($chr eq "Mt") ? "mt" : "other";
         next if $g eq "other";
         
         # 累积长度和相似性数据
         $len_sum{$id}{$g} += $len;
         $total_len{$id} += $len;
         $p_sum{$id} += $p;
         $p_cnt{$id}++;
         
         END{
             print "嵌合体ID\t类型\t核片段长度\t叶绿体片段长度\t线粒体片段长度\t平均相似性(%)\n";
             my ($nuc_cp, $nuc_mt, $cp_mt) = (0,0,0);
             
             for my $id (keys %len_sum) {
                 my $nuc = $len_sum{$id}{nuc} // 0;
                 my $cp = $len_sum{$id}{cp} // 0;
                 my $mt = $len_sum{$id}{mt} // 0;
                 next if $nuc+$cp+$mt == 0;  # 过滤无效数据
                 
                 my $avg_p = sprintf("%.1f", $p_sum{$id}/$p_cnt{$id});
                 my $type;
                 
                 # 判断嵌合体类型
                 if ($nuc>0 && $cp>0) {
                     $type = "核-叶绿体";
                     $nuc_cp++;
                 } elsif ($nuc>0 && $mt>0) {
                     $type = "核-线粒体";
                     $nuc_mt++;
                 } elsif ($cp>0 && $mt>0) {
                     $type = "叶绿体-线粒体";
                     $cp_mt++;
                 } else {
                     next;  # 仅属于一种基因组的排除
                 }
                 
                 print "$id\t$type\t$nuc\t$cp\t$mt\t$avg_p\n";
             }
             
             # 汇总统计
             my $total = $nuc_cp + $nuc_mt + $cp_mt;
             print "\n===== 嵌合体类型汇总 =====\n";
             print "核-叶绿体：$nuc_cp 条 (".($total>0?sprintf("%.1f",$nuc_cp/$total*100):0)."%)\n";
             print "核-线粒体：$nuc_mt 条 (".($total>0?sprintf("%.1f",$nuc_mt/$total*100):0)."%)\n";
             print "叶绿体-线粒体：$cp_mt 条 (".($total>0?sprintf("%.1f",$cp_mt/$total*100):0)."%)\n";
             print "总计：$total 条\n";
         }' chimeric_blast_verify.txt > chimeric_type_stat.txt

# 3. 查看最终统计结果
echo "===== 嵌合体分析最终结果 ====="
cat chimeric_type_stat.txt

```
# col assembly
```
#SRR611086
cd evaluation
mkdir col0
cd col0
cp ../SRR611086_8/* .
bash 0_script/0_master.sh 

cp 1_genome/genome.fa genome.fa 
faops filter -l 0 genome.fa stdout | grep -A 1 '^>Pt' | faops one -l 0 genome.fa Pt ref_Pt.fa

cp 7_merge_unitigs_bcalm/anchor.non-contained.fasta SRR611086.fa

minimap2 -ax sr ref_Pt.fa SRR611086.fa > SRR611086_aln.sam
samtools view -b SRR611086_aln.sam | samtools sort -o SRR611086_aln.sorted.bam
samtools index SRR611086_aln.sorted.bam
#从已排序的 BAM 文件计算每个位点的测序深度和碱基信息
samtools mpileup -aa -f ref_Pt.fa SRR611086_aln.sorted.bam > SRR611086_aln.pileup

#计算每一个位点的得分，得分标准为：若 contig 序列与参考序列一致，则得分 +1，否则得分为 -1
perl -F'\t' -ane '
BEGIN { print "Position\tScore\n" }
$ref = uc($F[2]);
$pos = $F[1];
$depth = $F[3];
$score = 0;
if ($depth == 0) { print "$pos\t0\n"; next }
$bases = $F[4];
$i = 0;
while ($i < length($bases)) {
    $c = substr($bases,$i,1);
    if ($c eq "." || $c eq ",") {
        if ($i+1 < length($bases) && (substr($bases,$i+1,1) eq "+" || substr($bases,$i+1,1) eq "-")) {
            $score -= 1;
            $i += 2;
            $lenstr = "";
            while ($i < length($bases) && substr($bases,$i,1) =~ /\d/) { $lenstr .= substr($bases,$i,1); $i++ }
            $indel_len = $lenstr+0;
            $i += $indel_len;
        } else { $score += 1; $i++ }
    } elsif ($c eq "^") { $i += 2 }
      elsif ($c eq "\$") { $i++ }  # 这里需要转义$
      elsif (uc($c) =~ /[ACGTN]/) { $score -= 1; $i++ }
      else { $i++ }
}
print "$pos\t$score\n";
' SRR611086_aln.pileup > SRR611086_per_base_scores.txt

perl -ane 'print if $.==1 || $F[1]<0' SRR611086_per_base_scores_perl.txt > positions_below_zero.txt
#28672   -6

#SRR5216995
#28672	-4
```
show in IGV:
![IGV_col0](../pic/IGV_col0.png)
```
scores <- read.table("SRR611086_per_base_scores.txt", header=TRUE, sep="\t")
library(ggplot2)

p <- ggplot(scores, aes(x=Position, y=Score)) +
  geom_point(aes(color = Score < 0), size=0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.3) + 
  scale_color_manual(values = c("black", "red"), guide = "none") +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "SRR611086",
       x = "Genomic Position",
       y = "Score")

ggsave("SRR611086_score.png", plot = p, width = 8, height = 4)
```
![SRR611086-score](../pic/SRR611086_score.png)
![SRR5216995-score](../pic/SRR5216995_score.png)

# Remove intermediate files
```
find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm
```
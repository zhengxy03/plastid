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

# Residual reads analysis
* pretreatment
```
SAMPLE="SRR616966"
FOLDS=(1 2 4 8 16 32 64)  # 与目录名中的fold值值完全匹配
BASE_DIR="$HOME/zxy/plastid/evaluation"
RESULT_DIR="$BASE_DIR/summary_results"
mkdir -p $RESULT_DIR

# 初始化汇总表格（调整统计项，移除高拷贝基因相关）
echo -e "fold\t总残留reads数\t转座元件reads数\t转座元件占比(%)\t多比对序列总数\t转座子多比对数\t转座子多比对占比(%)\t核平均read长度\t核chr1-read数\t核chr2-read数\t核chr3-read数\t核chr4-read数\t核chr5-read数\t线粒体平均覆盖度\t串联重复reads数\t串联重复占比(%)\t其他核残留reads数\t其他核残留占比(%)\t核-叶绿体嵌合体数\t核-线粒体嵌合体数\t叶绿体-线粒体嵌合体数\t总嵌合体数" > $RESULT_DIR/summary_stats.tsv

# 检查所有fold目录及3_bwa子目录中的BAM文件
for fold in "${FOLDS[@]}"; do
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    BWA_DIR="$FOLD_DIR/3_bwa"  # BAM文件所在的子目录
    
    if [ ! -d "$FOLD_DIR" ]; then
        echo "错误：未找到目录 $FOLD_DIR，请检查目录命名是否正确"
        exit 1
    fi
    
    # 检查3_bwa子目录中是否有R.sort.bam
    if [ ! -f "$BWA_DIR/R.sort.bam" ]; then
        echo "错误：$BWA_DIR 中未找到 R.sort.bam 文件"
        exit 1
    fi
done

# 选择fold=1的基因组数据作为共享资源
SHARED_GENOME_DIR="$BASE_DIR/${SAMPLE}_1/1_genome"

# 确保基因组索引存在
if [ ! -f "$SHARED_GENOME_DIR/genome.fa.fai" ]; then
    echo "生成基因组索引..."
    samtools faidx "$SHARED_GENOME_DIR/genome.fa"
fi

# 确保区域BED文件存在
if [ ! -f "$SHARED_GENOME_DIR/nuclear_regions.bed" ]; then
    echo "生成区域BED文件..."
    cut -f1,2 "$SHARED_GENOME_DIR/genome.fa.fai" > "$SHARED_GENOME_DIR/genome.chrom.sizes"
    perl -ne 'chomp; @f=split(/\t/); print "$f[0]\t0\t$f[1]\n" if $f[0]=~/^[1-5]$/' "$SHARED_GENOME_DIR/genome.chrom.sizes" > "$SHARED_GENOME_DIR/nuclear_regions.bed"
    perl -ne 'chomp; @f=split(/\t/); print "$f[0]\t0\t$f[1]\n" if $f[0] eq "Mt"' "$SHARED_GENOME_DIR/genome.chrom.sizes" > "$SHARED_GENOME_DIR/mito_regions.bed"
    perl -ne 'chomp; @f=split(/\t/); print "$f[0]\t0\t$f[1]\n" if $f[0] eq "Pt"' "$SHARED_GENOME_DIR/genome.chrom.sizes" > "$SHARED_GENOME_DIR/chloroplast_regions.bed"
fi

# 确保转座元件和串联重复注释存在
if [ ! -f "$SHARED_GENOME_DIR/nuclear_te.bed" ]; then
    echo "生成转座元件注释..."
    perl -ne 'chomp;
              next if /^Transposon_Name/;
              @f = split(/\t/);
              $te_id = $f[0];
              if ($te_id =~ /^AT(\d)TE/) { $chr = $1; }
              else { next; }
              $start = $f[2] - 1;
              $end = $f[3];
              $family = $f[5];
              $super_family = $f[6];
              $type = "$super_family($family)";
              if ($chr =~ /^[1-5]$/) {
                  print "$chr\t$start\t$end\t$type\t$te_id\n";
              }' "$SHARED_GENOME_DIR/TAIR10_Transposable_Elements.txt" > "$SHARED_GENOME_DIR/nuclear_te.bed"
fi

# 添加串联重复区域注释
if [ ! -f "$SHARED_GENOME_DIR/nuclear_tandem.bed" ]; then
    echo "警告：未找到串联重复注释文件，将使用空文件替代"
    touch "$SHARED_GENOME_DIR/nuclear_tandem.bed"
    # trf "$SHARED_GENOME_DIR/genome.fa" 2 5 7 80 10 50 500 -h -ngs > "$SHARED_GENOME_DIR/tandem_repeats.txt"
    # perl trf.pl "$SHARED_GENOME_DIR/tandem_repeats.txt" > "$SHARED_GENOME_DIR/nuclear_tandem.bed"
fi

if [ ! -f "$SHARED_GENOME_DIR/combined_ref_db.nhr" ]; then
    echo "构建BLAST数据库..."
    makeblastdb -in "$SHARED_GENOME_DIR/genome.fa" -dbtype nucl -out "$SHARED_GENOME_DIR/combined_ref_db" -parse_seqids
fi
```

* 转座子分析
```
for fold in "${FOLDS[@]}"; do
    echo "----- 处理fold $fold（目录：${SAMPLE}_${fold}） -----"
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    WORK_DIR="$FOLD_DIR/3_bwa"  # 工作目录和BAM文件目录一致
    mkdir -p $WORK_DIR
    cd $WORK_DIR

    # 从3_bwa子目录提取残留reads（BAM文件就在当前目录）
    samtools view -b -L "$SHARED_GENOME_DIR/nuclear_regions.bed" "$WORK_DIR/R.sort.bam" > nuclear_reads.bam
    samtools view -b -L "$SHARED_GENOME_DIR/mito_regions.bed" "$WORK_DIR/R.sort.bam" > mito_reads.bam
    samtools index nuclear_reads.bam
    samtools index mito_reads.bam

    # 转成bed文件用于覆盖分析
    bedtools bamtobed -i nuclear_reads.bam > nuclear_reads.bed
    bedtools bamtobed -i mito_reads.bam > mito_reads.bed

    # 统计残留reads的覆盖区域特征
    # 1. 核基因组残留reads的染色体分布
    perl -ne 'chomp; @f=split(/\t/); $count{$f[0]}++; END{
        foreach my $chr (1..5) {
            print "$chr\t".($count{$chr}//0)."\n";
        }
    }' nuclear_reads.bed > $WORK_DIR/nuclear_chrom_dist.txt
    
    chr1=$(grep "^1" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr2=$(grep "^2" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr3=$(grep "^3" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr4=$(grep "^4" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr5=$(grep "^5" $WORK_DIR/nuclear_chrom_dist.txt | cut -f2)
    chr1=${chr1:-0}; chr2=${chr2:-0}; chr3=${chr3:-0}; chr4=${chr4:-0}; chr5=${chr5:-0}
    
    # 2. 核基因组残留reads的平均长度
    avg_length=$(perl -ne 'chomp; @f=split(/\t/); $len=$f[2]-$f[1]; $sum+=$len; $count++; END{printf "%.2f", $sum/$count}' nuclear_reads.bed)
    
    # 3. 线粒体残留reads的平均覆盖度
    bedtools genomecov -i mito_reads.bed -g "$SHARED_GENOME_DIR/genome.chrom.sizes" -d | \
        perl -ne 'chomp; @f=split(/\t/); 
                 if($f[0] eq "Mt") { $sum += $f[2]; $count++ }
                 END{ printf "%.2f", $count>0 ? $sum/$count : 0 }' > $WORK_DIR/mito_avg_coverage.txt
    mito_avg_cov=$(cat $WORK_DIR/mito_avg_coverage.txt)

    # 转座子筛选分析
    samtools view -f 256 nuclear_reads.bam | perl -ne 'chomp; @f=split(/\t/); print "$f[0]\n"' | sort | uniq > nuclear_multimap_reads.txt
    
    bedtools intersect -abam nuclear_reads.bam -b "$SHARED_GENOME_DIR/nuclear_te.bed" -u > nuclear_te_reads.bam
    bedtools intersect -abam nuclear_reads.bam -b "$SHARED_GENOME_DIR/nuclear_te.bed" -v > nuclear_non_te_reads.bam
    
    bedtools bamtobed -i nuclear_reads.bam | grep -Ff nuclear_multimap_reads.txt > nuclear_multimap.bed
    bedtools intersect -a nuclear_multimap.bed -b "$SHARED_GENOME_DIR/nuclear_te.bed" -wa -wb | cut -f4 | sort | uniq > te_multimap_ids.txt
    
    # 计算统计值
    total_reads=$(samtools view -c nuclear_reads.bam)
    te_reads=$(samtools view -c nuclear_te_reads.bam)
    te_ratio=$(echo "scale=4; $te_reads/$total_reads*100" | bc | xargs printf "%.2f")
    multi_total=$(wc -l < nuclear_multimap_reads.txt)
    multi_te=$(wc -l < te_multimap_ids.txt)
    multi_te_ratio=$(echo "scale=4; $multi_te/$multi_total*100" | bc | xargs printf "%.2f")
    
    # 保存到临时文件
    echo -e "$fold\t$total_reads\t$te_reads\t$te_ratio\t$multi_total\t$multi_te\t$multi_te_ratio\t$avg_length\t$chr1\t$chr2\t$chr3\t$chr4\t$chr5\t$mito_avg_cov" > $WORK_DIR/repeat_stats.tmp
done
```
* 串联重复分析
```
for fold in "${FOLDS[@]}"; do
    echo "----- 处理fold $fold（目录：${SAMPLE}_${fold}） -----"
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    WORK_DIR="$FOLD_DIR/3_bwa"
    mkdir -p $WORK_DIR
    cd $WORK_DIR

    # 检查必要的中间文件
    if [ ! -f "nuclear_reads.bam" ] || [ ! -f "nuclear_non_te_reads.bam" ]; then
        echo "错误：未找到重复序列分析的中间文件，请先运行第一个循环"
        exit 1
    fi

    # 1. 串联重复区域的reads数（排除已归类为转座子的reads）
    bedtools intersect -abam nuclear_non_te_reads.bam -b "$SHARED_GENOME_DIR/nuclear_tandem.bed" -u > nuclear_tandem_reads.bam
    tandem_reads=$(samtools view -c nuclear_tandem_reads.bam)
    total_reads=$(samtools view -c nuclear_reads.bam)
    te_reads=$(samtools view -c nuclear_te_reads.bam)
    
    # 计算占比（确保保留两位小数）
    tandem_ratio=$(echo "scale=4; $tandem_reads/$total_reads*100" | bc | xargs printf "%.2f")
    
    # 2. 其他核残留reads数（未被上述类别包含）
    other_reads=$((total_reads - te_reads - tandem_reads))
    
    # 其他核残留占比 = 100% - 转座子占比 - 串联重复占比（确保保留两位小数）
    te_ratio=$(grep -w "$fold" $WORK_DIR/repeat_stats.tmp | cut -f4)
    other_ratio=$(echo "scale=4; 100 - $te_ratio - $tandem_ratio" | bc | xargs printf "%.2f")
    
    # 保存到临时文件
    echo -e "$tandem_reads\t$tandem_ratio\t$other_reads\t$other_ratio" > $WORK_DIR/other_nuclear_stats.tmp
done
```
* 嵌合体序列分析
```
for fold in "${FOLDS[@]}"; do
    echo "----- 处理fold $fold（目录：${SAMPLE}_${fold}） -----"
    FOLD_DIR="$BASE_DIR/${SAMPLE}_${fold}"
    WORK_DIR="$FOLD_DIR/3_bwa"
    mkdir -p $WORK_DIR
    cd $WORK_DIR

    # 检查必要的中间文件（新增nuclear_tandem_reads.bam的检查，确保串联重复分析已完成）
    if [ ! -f "nuclear_multimap_reads.txt" ] || [ ! -f "te_multimap_ids.txt" ] || [ ! -f "nuclear_tandem_reads.bam" ]; then
        echo "错误：未找到转座子/串联重复分析的中间文件，请先运行前两个循环"
        exit 1
    fi

    # 1：提取串联重复相关的多比对reads ID（排除串联重复干扰）
    # 先将串联重复reads转为BED格式，再筛选其中的多比对reads
    bedtools bamtobed -i nuclear_tandem_reads.bam > nuclear_tandem.bed
    # 从总多比对reads中筛选出串联重复相关的ID
    grep -Ff nuclear_multimap_reads.txt nuclear_tandem.bed | cut -f4 | sort | uniq > tandem_multimap_ids.txt
    rm -f nuclear_tandem.bed  # 清理临时文件

    # 2：筛选「非转座元件+非串联重复」的多比对ID
    comm -23 <(sort nuclear_multimap_reads.txt | uniq | sed 's/[[:space:]]//g') \
             <(sort -u <(cat te_multimap_ids.txt tandem_multimap_ids.txt)) \  # 合并转座子+串联重复的多比对ID
             > non_te_tandem_multimap_ids_clean.txt  # 最终用于嵌合体分析的候选ID
    rm -f tandem_multimap_ids.txt  # 清理临时文件

    # 步骤3：提取潜在嵌合体序列
    samtools view -h nuclear_reads.bam > full_reads.sam
    perl -ne '
        BEGIN {
            open(my $fh, "<", "non_te_tandem_multimap_ids_clean.txt") or die $!;
            while (my $line = <$fh>) {
                chomp $line;
                $ids{$line} = 1 if $line ne "";
            }
            close $fh;
        }
        if ($. == 1 || /^@/) { print; next; }
        my @fields = split(/\t/);
        print if exists $ids{$fields[0]};
    ' full_reads.sam > target_reads.sam
    samtools fasta target_reads.sam > non_te_tandem_multimap.fasta  # 修改输出文件名，避免与旧文件混淆
    rm -f full_reads.sam target_reads.sam non_te_tandem_multimap_ids_clean.txt  # 清理临时文件

    # 步骤4：BLAST筛选真实嵌合体
    blastn -query non_te_tandem_multimap.fasta \
           -db "$SHARED_GENOME_DIR/combined_ref_db" \
           -evalue 1e-20 \
           -outfmt "6 qseqid sseqid pident length evalue" \
           -max_target_seqs 3 \
           -num_threads 4 > blast_temp.txt

    # 步骤5：筛选嵌合体ID
    perl -ne 'chomp; @f=split(/\t/); next if @f<5;
             $id=$f[0]; $chr=$f[1]; $p=$f[2]; $e=$f[4];
             $g = ($chr=~/^[1-5]$/) ? "nuc" :
                  ($chr eq "Pt") ? "cp" :
                  ($chr eq "Mt") ? "mt" : "other";
             next if $g eq "other";
             if($p>=70 && $e<=1e-10) { push @{$h{$id}}, $g; }
             END{
                 for my $id (keys %h) {
                     my %unique = map { $_ => 1 } @{$h{$id}};
                     my @types = keys %unique;
                     if(@types >= 2) {
                         if( ($types[0] eq "nuc" && $types[1] eq "cp") || ($types[0] eq "cp" && $types[1] eq "nuc") ) {
                             print "$id\n";
                         } elsif( ($types[0] eq "nuc" && $types[1] eq "mt") || ($types[0] eq "mt" && $types[1] eq "nuc") ) {
                             print "$id\n";
                         } elsif( ($types[0] eq "cp" && $types[1] eq "mt") || ($types[0] eq "mt" && $types[1] eq "cp") ) {
                             print "$id\n";
                         }
                     }
                 }
             }' blast_temp.txt > candidate_chimeric_ids.txt
    rm -f blast_temp.txt  # 清理临时文件

    # 步骤6：提取并验证嵌合体序列
    sed 's/\/[12]$//' candidate_chimeric_ids.txt > candidate_ids_fixed.txt
    samtools view -h nuclear_reads.bam > full_reads2.sam
    perl -ne '
        BEGIN {
            open(my $fh, "<", "candidate_ids_fixed.txt") or die $!;
            while (my $line = <$fh>) {
                chomp $line;
                $ids{$line} = 1 if $line ne "";
            }
            close $fh;
        }
        if ($. == 1 || /^@/) { print; next; }
        my @fields = split(/\t/);
        print if exists $ids{$fields[0]};
    ' full_reads2.sam > target_reads2.sam
    samtools fasta target_reads2.sam > true_chimeric_seqs.fasta
    rm -f full_reads2.sam target_reads2.sam candidate_ids_fixed.txt candidate_chimeric_ids.txt

    # 步骤7：统计嵌合体类型
    blastn -query true_chimeric_seqs.fasta \
           -db "$SHARED_GENOME_DIR/combined_ref_db" \
           -evalue 1e-10 \
           -outfmt "6 qseqid sseqid length pident evalue" \
           -max_target_seqs 3 \
           -num_threads 4 > chimeric_blast_verify.txt

    perl -ne 'chomp; @f=split(/\t/); next if @f<5;
             $id=$f[0]; $chr=$f[1]; $len=$f[2]; $p=$f[3];
             $g = ($chr=~/^[1-5]$/) ? "nuc" :
                  ($chr eq "Pt") ? "cp" :
                  ($chr eq "Mt") ? "mt" : "other";
             next if $g eq "other";
             $len_sum{$id}{$g} += $len;
             $p_sum{$id} += $p;
             $p_cnt{$id}++;
             END{
                 my ($nuc_cp, $nuc_mt, $cp_mt) = (0,0,0);
                 for my $id (keys %len_sum) {
                     my $nuc = $len_sum{$id}{nuc} // 0;
                     my $cp = $len_sum{$id}{cp} // 0;
                     my $mt = $len_sum{$id}{mt} // 0;
                     next if $nuc+$cp+$mt == 0;
                     
                     if ($nuc>0 && $cp>0) {
                         $nuc_cp++;
                     } elsif ($nuc>0 && $mt>0) {
                         $nuc_mt++;
                     } elsif ($cp>0 && $mt>0) {
                         $cp_mt++;
                     }
                 }
                 my $total = $nuc_cp + $nuc_mt + $cp_mt;
                 print "$nuc_cp\t$nuc_mt\t$cp_mt\t$total";
             }' chimeric_blast_verify.txt > $WORK_DIR/chimera_stats.tmp
    rm -f chimeric_blast_verify.txt true_chimeric_seqs.fasta non_te_tandem_multimap.fasta  
done

# 汇总结果
for fold in "${FOLDS[@]}"; do
    WORK_DIR="$BASE_DIR/${SAMPLE}_${fold}/3_bwa"
    repeat_stats=$(cat $WORK_DIR/repeat_stats.tmp)
    other_nuclear_stats=$(cat $WORK_DIR/other_nuclear_stats.tmp)
    chimera_stats=$(cat $WORK_DIR/chimera_stats.tmp)
    echo -e "$repeat_stats\t$other_nuclear_stats\t$chimera_stats" >> $RESULT_DIR/summary_stats.tsv
    rm -f $WORK_DIR/repeat_stats.tmp $WORK_DIR/other_nuclear_stats.tmp $WORK_DIR/chimera_stats.tmp
done


```
* 可视化
```
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(scales)

# 设置中文字体支持
if (.Platform$OS.type == "windows") {
  windowsFonts(
    sans = windowsFont("SimHei"),
    serif = windowsFont("SimSun")
  )
} else if (Sys.info()["sysname"] == "Darwin") {  # macOS
  theme_update(text = element_text(family = "Heiti TC"))
} else {  # Linux
  theme_update(text = element_text(family = "WenQuanYi Micro Hei", size = 12))
}

# 读取数据
result_file <- "~/zxy/plastid/evaluation/summary_results/summary_stats.tsv"
df <- read_tsv(result_file, show_col_types = FALSE)

# 数据处理：转换为长格式并整理标签（注意使用反引号处理特殊字符）
plot_data <- df %>%
  select(
    fold = fold,
    `转座元件占比(%)`,
    `串联重复占比(%)`,
    `其他核残留占比(%)`
  ) %>%
  pivot_longer(
    cols = -fold,
    names_to = "repeat_type",
    values_to = "proportion"
  ) %>%
  mutate(
    fold = factor(fold, levels = unique(fold)),  # 保持fold原始顺序
    repeat_type = factor(repeat_type, 
                        levels = c("转座元件占比(%)", "串联重复占比(%)", "其他核残留占比(%)")) %>%
      # 去除重复的"占比(%)"文字，使标签更简洁
      recode(
        "转座元件占比(%)" = "转座元件",
        "串联重复占比(%)" = "串联重复",
        "其他核残留占比(%)" = "其他核残留"
      )
  )

# 定义颜色和参数
colors <- c(
  "转座元件" = "black",
  "串联重复" = "gray60",
  "其他核残留" = "white"
)
bar_width <- 0.6
dodge_width <- 0.6

# 绘制柱状图
p <- ggplot(plot_data, 
           aes(x = fold, y = proportion, fill = repeat_type, group = repeat_type)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = dodge_width),
           color = "black",
           width = bar_width) +
  # 添加数值标签
  geom_text(aes(label = sprintf("%.2f%%", proportion)),
            position = position_dodge(width = dodge_width),
            vjust = -0.5,
            size = 2.5,
            fontface = "bold") +
  scale_fill_manual(values = colors, name = "核残留reads类型") +
  scale_y_continuous(
    labels = percent_format(scale = 1, accuracy = 0.1),
    limits = c(0, max(plot_data$proportion) * 1.2),
    expand = c(0, 0)
  ) +
  labs(
    x = "Fold",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10, "mm")
  )

# 保存图形
ggsave("repeat_proportion_by_fold.png",
       plot = p,
       width = 12,
       height = 8,
       dpi = 300,
       bg = "white")

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
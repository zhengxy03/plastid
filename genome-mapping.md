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
* `cutoff = FOLD * DEPTH`
FOLD 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64
* 叶绿体的拷贝数>线粒体>核基因组
# Symlink
```
mkdir -p ~/data/plastid/evaluation
cd ~/data/plastid/evaluation

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
    SRR="${item%%::*}" #删除最长匹配的后缀
    GENOME="${item##*::}" #删除最长匹配的前缀

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

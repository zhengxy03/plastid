# Model bacteria
* template
```
Usage: anchr template [OPTIONS]

Options:
      --genome <genome>        Your best guess of the haploid genome size [default: 1000000]
      --se                     Single end mode
      --xmx <xmx>              Set Java memory usage
  -p, --parallel <parallel>    Number of threads [default: 8]
      --queue <queue>          Queue name of the LSF cluster
      --repetitive             Find repetitive regions
      --fastqc                 Run FastQC
      --fastk                  Run FastK
      --insertsize             Calc insert sizes
      --reads <reads>          How many reads to estimate insert sizes [default: 1000000]
      --trim <trim>            Opts for trim [default: --dedupe]
      --sample <sample>        Sampling coverage
      --qual <qual>            Quality threshold [default: "25 30"]
      --len <len>              Filter reads less or equal to this length [default: 60]
      --filter <filter>        Adapter, artifact, or both [default: adapter]
      --quorum                 Run quorum
      --merge                  Run merge reads
      --prefilter <prefilter>  Prefilter=N (1 or 2) for tadpole and bbmerge, 1 use less memories
      --ecphase <ecphase>      Error-correct phases. Phase 2 can be skipped [default: "1 2 3"]
      --bwa <bwa>              Map trimmed reads to the genome
      --gatk                   Calling variants with GATK Mutect2
      --cov <cov>              Down sampling coverages [default: "40 80"]
  -u, --unitigger <unitigger>  Unitigger used: bcalm, bifrost, superreads, or tadpole [default: bcalm]
      --splitp <splitp>        Parts of splitting [default: 20]
      --statp <statp>          Parts of stats [default: 2]
      --readl <readl>          Length of reads [default: 100]
      --uscale <uscale>        The scale factor for upper, (median + k * MAD) * u [default: 2] #定义kmer深度阈值上限
      --lscale <lscale>        The scale factor for upper, (median - k * MAD) / l [default: 3]
      --redo                   Redo anchors when merging anchors
      --extend                 Extend anchors with other contigs
      --gluemin <gluemin>      Min length of overlaps to be glued [default: 30]
      --fillmax <fillmax>      Max length of gaps [default: 100]
      --busco                  Run busco
  -h, --help                   Print help
  -V, --version                Print version


* Info
    * --genome
    * --se
    * --xmx
    * --parallel 8
    * --queue mpi

* Resources
    * --repetitive

* Quality check
    * --fastqc
    * --fastk
    * --insertsize
    * --reads 1000000

* Trimming
    * --trim "--dedupe"
    * --sample "300"
    * --qual "25 30"
    * --len "60"
    * --filter "adapter"

* Post-trimming
    * --quorum
    * --merge
    * --prefilter
    * --ecphase "1 2 3"

* Mapping
    * --bwa
    * --gatk

* Down sampling, unitigs, and anchors
    * --cov "40 80"
    * --unitigger "bcalm"
    * --splitp 20
    * --statp 2
    * --readl 100
    * --uscale 2
    * --lscale 3
    * --redo

* Extend anchors
    * --extend
    * --gluemin 30
    * --fillmax 100

* Validate assemblies
    * --busco
```
## *Mycoplasma genitalium* G37
> g37:生殖支原体，基因组最小的支原体<br>
`g37: reference`
```
mkdir -p g37/1_genome
cd g37/1_genome

curl -O "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/556/925/GCF_040556925.1_ASM4055692v1/GCF_040556925.1_ASM4055692v1_genomic.fna.gz"
gzip -d GCF_040556925.1_ASM4055692v1_genomic.fna.gz
mv GCF_040556925.1_ASM4055692v1_genomic.fna genome.fa
sed -i '1s/.*/>1/' genome.fa
```
`g37: download`

* illumina
```
mkdir -p 2_illumina
cd 2_illumina

ln -s ../ERR486835_1.fastq.gz R1.fq.gz
ln -s ../ERR486835_2.fastq.gz R2.fq.gz
```
`g37: template`
```
anchr template \
    --genome 580076 \
    --parallel 8 \
    --xmx 12g \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "bcalm bifrost superreads tadpole" \
    --statp 2 \
    --readl 125 \
    --uscale 2 \
    --lscale 3 \
    --redo
```
`g37: run`
```
bash 0_script/1_repetitive.sh

bash 0_script/0_master.sh
```
具体过程：
* 质控
* 采样组装
> kunitigs、  tadpole 的组装以及其他组装方法包括基因组组装软件 SPAdes、Megahit、Platanus、SPAdesMR、MegahitMRM(组装基于`De Bruijn`图从头组装)<br>

`De Bruijn`图:
* 定义：
是一种有向图（directed graph），常用来表示所有长度为 k 的序列（k-mer）之间的重叠关系。
    * 节点是长度 k-1 的所有子串
    * 边是这些 k-mers 的滑动窗口连接
    * 在图上找路径 → 得到 contig/scaffold → 组装出完整基因组
* 例子：
```
AGC -> GC
GCT -> CT
CTG -> TG
```
连接成：
```
AG -> GC -> CT -> TG
```
![De Brujin graph](./pic/de%20Brujin%20graph.jpg "De Brujin graph")
* bubble（泡）
    * 定义
> Bubble（泡） 是在 De Bruijn 图里，一对节点之间存在两条或多条并行、可替代的路径。必须被组装器识别和正确处理
    * 主要成因有：
        ⭐ 测序错误 <br>
        ⭐ 真正的变异（SNP、小Indel）<br>
        ⭐ 杂合性（heterozygosity）<br>
        ⭐ 重复序列（repeats）<br>
    * 如果不解决：组装就会分裂出很多碎片；会误把错误或杂合区当成不同的 contig
* unitig:
    * 定义：
    > 在 De Bruijn 图里，没有分支的、可唯一延伸的路径段。
    * `anchr基于unitig进行组装`

## *Escherichia coli* str. K-12 substr. MG1655
> mg1655:野生型菌株，具备完整的野生型基因背景<br>
`mg1655: reference`
```
mkdir -p mg1655/1_genome
cd mg1655/1_genome

curl -O "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
gzip -d GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna genome.fa
sed -i '1s/.*/>1/' genome.fa
```
`mg1655: download`
```
aria2c -x 9 -s 3 -c ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz
aria2c -x 9 -s 3 -c ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz
```
* illumina
```
mkdir -p 2_illumina
cd 2_illumina

ln -s ../MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz R1.fq.gz
ln -s ../MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz R2.fq.gz
```
`mg1655: template`
* rsync to hpcc
```
rsync -avP \
    ~/project/anchr/mg1655/ \
    wangq@202.119.37.251:zxy/mg1655
```
* template
```
anchr template \
    --genome 4641652 \
    --parallel 24 \
    --xmx 80g \
    --queue mpi \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --tile --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --bwa "Q25L60" \
    --gatk \
    \
    --cov "40 80" \
    --unitigger "bcalm bifrost superreads" \
    --statp 2 \
    --readl 151 \
    --uscale 2 \
    --lscale 3 \
    \
    --extend \
    \
    --busco
```
`mg1655: run`
```
bash 0_script/1_repetitive.sh
#bash 0_script/0_master.sh
bsub -q mpi -n 24 -J "mg1655-0_master" "bash 0_script/0_master.sh"
```
## *Escherichia coli* str. K-12 substr. DH5alpha
> dh5alpha:经过诱变的K12菌株，常用于蓝白斑筛选<br>
`dh5alpha: reference`
```
mkdir -p dh5alpha/1_genome
cd dh5alpha/1_genome

curl -O "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/899/475/GCF_002899475.1_ASM289947v1/GCF_002899475.1_ASM289947v1_genomic.fna.gz"
gzip -d GCF_002899475.1_ASM289947v1_genomic.fna.gz
mv GCF_002899475.1_ASM289947v1_genomic.fna genome.fa
sed -i '1s/.*/>1/' genome.fa
```
`dh5alpha: download`
```
cat << EOF > source.csv
SRP251726,dh5alpha,HiSeq 2500 PE125
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

rgr md ena_info.tsv --fmt

aria2c -x 9 -s 3 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt
```
* rsync to hpcc
```
rsync -avP \
    ~/project/anchr/dh5alpha/ \
    wangq@202.119.37.251:zxy/dh5alpha
```
* illumina
```
mkdir -p 2_illumina
cd 2_illumina

ln -s ../SRR11245239_1.fastq.gz R1.fq.gz
ln -s ../SRR11245239_2.fastq.gz R2.fq.gz
```

`mg1655: template`
```
anchr template \
    --genome 4583637 \
    --parallel 24 \
    --xmx 80g \
    --queue mpi \
    \
    --repetitive \
    \
    --fastqc \
    --insertsize \
    --fastk \
    \
    --trim "--dedupe --cutoff 30 --cutk 31" \
    --qual "25 30" \
    --len "60" \
    --filter "adapter artifact" \
    \
    --quorum \
    --merge \
    --ecphase "1 2 3" \
    \
    --cov "40 80" \
    --unitigger "bcalm bifrost superreads" \
    --statp 2 \
    --readl 125 \
    --uscale 2 \
    --lscale 3 \
    --redo \
    \
    --extend \
    \
    --busco
```
`dh5alpha: run`
```
bash 0_script/1_repetitive.sh
bsub -q mpi -n 24 -J "dh5alpha-0_master" "bash 0_script/0_master.sh"
```
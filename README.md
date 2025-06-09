# anchr

# download anchr and dependencies
## anchr
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

cargo install --git https://github.com/wang-q/anchr --branch main
```
## runtime dependencies
```
#cbp
curl -LO https://github.com/wang-q/cbp/releases/latest/download/cbp.linux
chmod +x cbp.linux
./cbp.linux init
source ~/.bashrc

cbp install openjdk jq parallel pigz
cbp install fastqc sickle bwa samtools picard
cbp install spoa
cbp install bcalm bifrost
cbp install tsv-utils faops intspan
cbp install dazzdb daligner

# linux
cbp install mosdepth

brew tap brewsci/bio
brew tap wang-q/tap


FILES=$(curl -fsSL https://api.github.com/repos/wang-q/builds/git/trees/master?recursive=1 |
        jq -r '.tree[] | select( .path | startswith("tar/") ) | .path' |
        sed 's|^tar/||')
cd $HOME/bin/

for file in $FILES; do
    echo "download: $file"
    
    curl -fsSL --retry 5 --retry-delay 2 \
        "https://github.com/wang-q/builds/raw/master/tar/$file" \
        -o "$file"

done

for file in *.linux.tar.gz; do
    echo "tar: $file"
    tar xvzf "$file"
done

curl -LO https://github.com/brentp/mosdepth/releases/download/v0.3.10/mosdepth
chmod +x mosdepth


brew install wang-q/tap/tsv-utils

anchr dep install | bash
anchr dep check | bash
#dazz
#poa

brew install binutils
brew link binutils --force
parallel -j 1 -k --line-buffer '
    Rscript -e '\'' if (!requireNamespace("{}", quietly = FALSE)) { install.packages("{}", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN") } '\''
    ' ::: \
        argparse minpack.lm \
        ggplot2 scales viridis

# Optional
# quast - assembly quality assessment
curl -LO https://github.com/ablab/quast/releases/download/quast_5.3.0/quast-5.3.0.tar.gz
tar xvfz quast-5.3.0.tar.gz
cd quast-5.3.0
python3 ./setup.py install

# Optional: leading assemblers
brew install brewsci/bio/megahit

brew install spades
```
## examples
### Individual subcommands
```
mkdir -p anchr/test
cd anchr/test
#lambda
for F in R1.fq.gz R2.fq.gz; do
    1>&2 echo ${F}
    curl -fsSLO "https://raw.githubusercontent.com/wang-q/anchr/main/tests/Lambda/${F}"
done
```
修剪
* trim
```
Usage: anchr trim [OPTIONS] <infiles>...

Arguments:
  <infiles>...  Sets the input file to use

Options:
  -q, --qual <qual>          Quality threshold [default: 25]
  -l, --len <len>            Filter reads less or equal to this length [default: 60]
      --filter <filter>      Adapter, artifact, or both [default: adapter]
      --trimq <trimq>        Quality score for 3' end [default: 15]
      --trimk <trimk>        Kmer for 5' adapter trimming [default: 23]
      --matchk <matchk>      Kmer for decontamination [default: 27]
      --cutk <cutk>          Kmer for cutoff [default: 31]
      --adapter <adapter>    The adapter file
      --artifact <artifact>  The artifact file
      --prefix <prefix>      Prefix of trimmed reads [default: R]
      --dedupe               Do the dedupe step
      --tile                 With normal Illumina names, do tile-based filtering
      --cutoff <cutoff>      Min kmer depth cutoff
      --sample <sample>      The sampling step
      --xmx <xmx>            Set Java memory usage
  -p, --parallel <parallel>  Number of threads [default: 8]
  -o, --outfile <outfile>    Output filename. [stdout] for screen [default: trim.sh]
  -h, --help                 Print help
  -V, --version              Print version


<file1> [file2]

Fastq files can be gzipped
```
```
mkdir -p trim
pushd trim #切换目录并打印当前目录

anchr trim \
    ../R1.fq.gz ../R2.fq.gz \
    -q 25 -l 60 \
    -o stdout |
    bash #将 anchr trim 的输出作为命令传递给 bash 执行
popd #返回上一级目录并打印跳转后的目录
```
合并
* mergeread
```
Usage: anchr mergeread [OPTIONS] <infiles>...

Arguments:
  <infiles>...  Sets the input file to use

Options:
  -l, --len <len>              Filter reads less or equal to this length [default: 60]
  -q, --qual <qual>            Quality score for 3' end [default: 15]
      --prefilter <prefilter>  Prefilter=N (1 or 2) for tadpole and bbmerge
      --ecphase <ecphase>      Error-correct phases. Phase 2 can be skipped [default: "1 2 3"] #ecphase：指定纠错阶段（Error Correction Phase）
                                    #阶段 1：基于 k-mer 频率的纠错（高丰度 k-mer 更可信）
                                    #阶段 2：利用配对末端的互补信息相互验证
                                    #阶段 3：参考基因组辅助纠错（如果提供）
      --prefixm <prefixm>      Prefix of merged reads [default: M]
      --prefixu <prefixu>      Prefix of unmerged reads [default: U]
      --xmx <xmx>              Set Java memory usage
  -p, --parallel <parallel>    Number of threads [default: 8]
  -o, --outfile <outfile>      Output filename. [stdout] for screen [default: merge.sh]
  -h, --help                   Print help
  -V, --version                Print version


<R1> [R2] [Rs]

Fastq files can be gzipped
R1 and R2 are paired; or R1 is interleaved
Rs is single
```
```
mkdir -p merge
pushd merge

anchr mergeread \
    ../trim/R1.fq.gz ../trim/R2.fq.gz ../trim/Rs.fq.gz #双端数据中无法配对的序列被单独保存为 Rs
    --ecphase "1 2 3" \
    --parrallel 4 \
    -o stdout |
    bash
popd
```
丢弃
* quorum
```
Run quorum to discard bad reads

Usage: anchr quorum [OPTIONS] <infiles>...

Arguments:
  <infiles>...  Sets the input file to use

Options:
      --coverage <jf>        Jellyfish hash size [default: 500000000]
      --estsize <estsize>    Estimated genome size [default: auto]
      --prefix <prefix>      Prefix of .cor.fa.gz [default: pe]
  -p, --parallel <parallel>  Number of threads [default: 8]
  -o, --outfile <outfile>    Output filename. [stdout] for screen [default: quorum.sh]
  -h, --help                 Print help
  -V, --version              Print version


<PE file1> [PE file2] [SE file]

Fastq files can be gzipped
```
```
pushd trim
anchr quorum \
    R1.fq.gz R.fq.gz \
    -o stdout |
    bash
popd

pushd trim/Q25L60
anchr quorum \
    R1.fq.gz R2.fq.gz Rs.fq.gz \
    -o stdout |
    bash
popd
```
构建单元序列unitigs
* unitigs
```
Options:
  -u, --unitigger <unitigger>  Which unitig constructor to use: bcalm, bifrost, superreads, or tadpole [default: superreads]
                                - bcalm     # 基于Bloom Filter的算法, 快速判断一个 k-mer 是否存在于测序数据中
                                - bifrost   # 图形化构建工具, 处理多样本数据，可识别样本间的变异
                                - superreads # 基于k-mer的超读长算法, 将重叠的短读长合并为更长的连续序列，减少后续组装的复杂度
                                - tadpole   # BBtools中的快速组装工具, 贪心延伸：从种子 k-mer 开始，每次选择最可能的延伸路径，直到无法继续
      --estsize <estsize>      Estimated genome size [default: auto]
      --kmer <kmer>            K-mer size to be used [default: 31]
      --min <min>              Minimal length of unitigs [default: 1000]
      --merge                  Merge unitigs from all k-mers
  -p, --parallel <parallel>    Number of threads [default: 8]
  -o, --outfile <outfile>      Output filename. [stdout] for screen [default: unitigs.sh]
  -h, --help                   Print help
  -V, --version                Print version

#成对末端测序 reads 
#存储流程运行所需的环境参数或配置信息
<pe.cor.fa> <env.json>

Fasta files can't be gzipped
```
bcalm
```
mkdir -p bcalm
pushd bcalm

anchr unitigs \
    ../trim/pe.cor.fa ../trim/env.json \
    -u bcalm \
    --kmer "31 41 51 61 71 81" \
    --parrallel 4 \
    -o unitigs.sh
bash unitigs.sh
popd
```
锚定
* anchors
```
#用于从 单元型序列（unitigs） 和 测序读长（reads） 中识别锚定序列（Anchors）
Usage: anchr anchors [OPTIONS] <infiles>...

Arguments:
  <infiles>...  Sets the input file to use

Options:
      --min <min>            Minimal length of anchors [default: 1000]
      --mincov <mincov>      Minimal coverage of reads [default: 5]
      --readl <readl>        Length of reads [default: 100]
      --mscale <mscale>      The scale factor for MAD, median +/- k * MAD [default: 3]
      --lscale <lscale>      The scale factor for lower, (median - k * MAD) / l [default: 3] #异常值检测参数（基于中位数绝对偏差 MAD）
      --uscale <uscale>      The scale factor for upper, (median + k * MAD) * u [default: 2]
      --fill <fill>          Fill holes short than or equal to this [default: 1]
      --ratio <ratio>        Fill large holes (opt.fill * 10) when covered ratio larger than this [default: 0.98]
      --longest              Only keep the longest proper region
      --keepedge             Keep edges of anchors
  -p, --parallel <parallel>  Number of threads [default: 8]
  -o, --outfile <outfile>    Output filename. [stdout] for screen [default: anchors.sh]
  -h, --help                 Print help
  -V, --version              Print version


<contig.fasta> <pe.cor.fa> [more reads]

Fasta files can‘t be gzipped

To get single-copy regions, set --uscale to 1.5
```
```
mkdir bcalm/anchors
pushd bcalm/anchors

anchr anchors \
    ../unitigs.fasta \
    ../pe.cor.fa \
    --readl 150 \
    --keepedge \
    -p 4 \
    -o anchors.sh
bash anchors.sh
popd
```

unitigs -superreads
```
gzip -dcf trim/pe.cor.fa.gz > trim/pe.cor.fa
mkdir -p superreads
pushd superreads

anchr unitigs \
    ../trim/pe.cor.fa ../trim/env.json \
    --kmer "31 41 51 61 71 81" \
    --parallel 4 \
    -o unitigs.sh
bash unitigs.sh
popd
```
unitigs -tadpole
```
mkdir -p tadpole
pushd tadpole

anchr unitigs \
    ../trim/pe.cor.fa ../trim/env.json \
    -u tadpole \
    --kmer "31 41 51 61 71 81" \
    --parallel 4 \
    -o unitigs.sh
bash unitigs.sh
popd
```
### fetch data
* ena
```
Usage: anchr ena [OPTIONS] <infile>

Arguments:
  <infile>  Sets the input file to use

Options:
  -o, --outfile <outfile>  Output filename. [stdout] for screen [default: stdout]
  -h, --help               Print help
  -V, --version            Print version


* info - Grab information from ENA
* prep - Create downloading scripts
```
example:
```
cat << EOF > source.csv
ERX452667,G37,MiSeq
EOF

anchr ena info | perl - -v source.csv > ena_info.yml #配置文件/信息交互
anchr ena prep | perl - ena_info.yml

rgr md ena_info.tsv --fmt #markdown

aria2c -j 4 -x 4 -s 2 -c --file-allocation=none -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt
```

### template

### Overlaps-standalone
```
mkdir ovlpr
pushd ovlpr
curl -s https://api.github.com/repos/wang-q/anchr/contents/tests/ovlpr \
  > files.json

FILENAMES=$(cat files.json | jq -r '.[] | select(.type == "file") | .name')
for F in $FILENAMES; do
    echo "下载文件: $F"
    curl -fsSLO "https://raw.githubusercontent.com/wang-q/anchr/main/tests/ovlpr/$F"
done
```

* dazzname
```
Usage: anchr dazzname [OPTIONS] <infiles>...

Arguments:
  <infiles>...  Set the input files to use

Options:
      --prefix <prefix>    Prefix of record names [default: read]
      --start <start>      Starting index [default: 1]
      --no-replace         Do not write a .replace.tsv
  -o, --outfile <outfile>  Output filename. [stdout] for screen [default: stdout]
  -h, --help               Print help
  -V, --version            Print version
```
* show2ovlp
```
Usage: anchr show2ovlp [OPTIONS] <show.txt> <replace.tsv>

Arguments:
  <show.txt>     Set the input file to use
  <replace.tsv>  Set the input file to use

Options:
      --orig               Original names of sequences
  -o, --outfile <outfile>  Output filename. [stdout] for screen [default: stdout]
  -h, --help               Print help
  -V, --version            Print version


f_id and g_id are integers, --orig convert them to the original ones
```
* paf2ovlp

example
```
anchr dazzname 1_4.anchor.fasta -o stdout
anchr show2ovlp 1_4.show.txt 1_4.replace.tsv --orig
anchr paf2ovlp 1_4.pac.paf
```
* covered
```
Usage: anchr covered [OPTIONS] <infiles>...

Arguments:
  <infiles>...  Set the input files to use

Options:
  -c, --coverage <coverage>  minimal coverage [default: 3]
  -l, --len <len>            minimal length of overlaps [default: 1000]
  -i, --idt <idt>            minimal identities of overlaps [default: 0.0]
      --paf                  PAF as input format
      --longest              only keep the longest span
      --base                 per base coverage
      --mean                 mean coverage
  -o, --outfile <outfile>    Output filename. [stdout] for screen [default: stdout]
  -h, --help                 Print help
  -V, --version              Print version
```
```
echo "1_4.anchor.fasta;1_4.pac.fasta" |
    parallel --colsep ";" -j 1 "
        minimap2 -cx asm20 {1} {2} |
            anchr paf2ovlp stdin |
            tsv-sort
        minimap2 -cx asm20 {2} {1} |
            anchr paf2ovlp stdin |
            tsv-sort
    " |
    anchr covered stdin --mean

anchr covered 1_4.pac.paf.ovlp.tsv
anchr covered 11_2.long.paf --paf
anchr covered 1_4.pac.paf.ovlp.tsv --base
anchr covered 1_4.pac.paf.ovlp.tsv --mean
```
组装
```
anchr restrict 1_4.ovlp.tsv 1_4.restrict.tsv
```
### Overlap-dalinger
* overlap
```
Usage: anchr overlap [OPTIONS] <infiles>...

Arguments:
  <infiles>...  Set the input files to use

Options:
  -l, --len <len>            minimal length of overlaps [default: 500]
  -i, --idt <idt>            minimal identities of overlaps [default: 0.7]
      --serial               Serials instead of original names in outputs  #在输出中使用序列号代替原始序列名称
      --all                  All overlaps instead of proper ones #输出所有重叠，而非仅保留 “适当重叠”
  -p, --parallel <parallel>  Number of threads [default: 8]
  -o, --outfile <outfile>    Output filename. [stdout] for screen [default: stdout]
  -h, --help                 Print help
  -V, --version              Print version


* This command is for small files.
* All operations are running in a tempdir and no intermediate files are kept.
```
```
anchr overlap 1_4.pac.fasta
anchr overlap 1_4.pac.fasta --idt 0.8 --len 2500 --serial
```
* orient
```
Usage: anchr orient [OPTIONS] <infiles>...

Arguments:
  <infiles>...  Set the input files to use

Options:
  -l, --len <len>            minimal length of overlaps [default: 1000]
  -i, --idt <idt>            minimal identities of overlaps [default: 0.85]
  -r, --restrict <restrict>  Restrict to known pairs
  -p, --parallel <parallel>  Number of threads [default: 8]
  -o, --outfile <outfile>    Output filename. [stdout] for screen [default: stdout]
  -h, --help                 Print help
  -V, --version              Print version


* This command is for small files.
* All operations are running in a tempdir and no intermediate files are kept.
```
```
anchr orient 1_4.anchor.fasta 1_4.pac.fasta
anchr orient 1_4.anchor.fasta 1_4.pac.fasta -r 1_4.2.restrict.tsv
```
* contained
```
anchr contained contained.fasta
```
* merge
```
anchr merge merge.fasta -o test.fasta
```
* overlap2 #必须是两个FASTA 文件
```
Usage: anchr overlap2 [OPTIONS] <infile> <infile2>

Arguments:
  <infile>   Sets the first input file
  <infile2>  Sets the second input file

Options:
  -d, --dir <dir>            Change working directory [default: .]
      --p1 <p1>              Prefix of the first [default: anchor]
      --p2 <p2>              Prefix of the second [default: long]
      --pd <pd>              Prefix of the result files
  -b, --block <block>        Block size in Mbp [default: 20]
  -l, --len <len>            minimal length of overlaps [default: 500]
  -i, --idt <idt>            minimal identities of overlaps [default: 0.7]
      --all                  All overlaps instead of proper ones
  -p, --parallel <parallel>  Number of threads [default: 8]
  -h, --help                 Print help
  -V, --version              Print version


* Since `daligner` cannot perform overlap comparisons between two databases,
  the current strategy is to place the two FASTA files into a single database
  and then split that database into multiple chunks. The chunk containing the
  first sequence will be compared against the other chunks. This approach can
  still significantly reduce the computational workload.

* The inputs can't be stdin
* All intermediate files (.fasta, .replace.tsv, .db, .las, .show.txt, .ovlp.tsv)
  are kept in the working directory.

* .ovlp.tsv records the serials
```
```
anchr overlap2 1_4.anchor.fasta 1_4.pac.fasta -d tmp
```
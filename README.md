# PolyA-miner<img src = "img/PolyA-miner_logo.png" width="180"> 
Accurate assessment of differential alternative poly-adenylation from 3'Seq data using vector projections and non-negative matrix factorization -Yalamanchili H.K. et al.

## Usage
```
python3 PolyA-miner.py <options>
  ```

## Arguments  
  ```
 Required arguments:
  -mode {bam,fastq}     Run mode options: 'bam' to start from mapped data, 'fastq' to start from raw data (default: bam)
  -c1 C1 [C1 ...]       Comma-separated list of condition1 files. Full path for BAMs (index files are also expected) or Just file names for fastq (default: None)
  -c2 C2 [C2 ...]       Comma-separated list of condition2 files. Full path for BAMs (index files are also expected) or Just file names for fastq (default: None)
  -fasta FASTA          Reference fasta sequence (default: None)
  -bed BED              Reference genes bed file (default: None)
  -pa PA                PolyA annotations file standard 6 column bed format (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -d D                  Base directory of input fastq files. Valid for -mode fastq (default: None)
  -o O                  Out put directory (default: PolyAminer_OUT)
  -s {0,1,2}            Strand information 0: un-stranded 1: fwd-strand 2:rev-strand. (default: 0)
  -index INDEX          Reference genome bowtie2 index. Valid for -mode fastq (default: None)
  -umi UMI              Length of UMIs, 0 if not used (default: 0)
  -apaBlock APABLOCK    Window size for annotated polyA sites (default: 30)
  -mdapa MDAPA          Cluster distance for annotated polyA sites: Merge polyA sites with in this distance. (default: 0)
  -md MD                Cluster distance for de-novo polyA sites: Merge polyA sites with in this distance (default: 0)
  -anchor ANCHOR        Overlap in "bp" for mapping de-novo polyA sites to annotated polyA sites (default: 1)
  -expNovel {1,0}       Explore novel APA sites 0: only annotated sites 1: de-novo (default: 0)
  -novel_d NOVEL_D      Distance from annotated TES to map novel pA sites (default: 1000)
  -p P                  No. of processors to use (default: 4)
  -ip IP                Internal priming window (default: 50)
  -a A                  Internal priming polyA fraction (default: 0.65)
  -pa_p PA_P            pOverA filter: P (default: 0.6)
  -pa_a PA_A            pOverA filter: A (default: 5)
  -pa_m PA_M            pOverA filter: M (default: 2)
  -gene_min GENE_MIN    Min counts per Gene (default: 10)
  -apa_min APA_MIN      Min. proportion per APA (default: 0.05)
  -t {BB,iNMF}          Statistical Test- BB: for beta-binomial or iNMF: for iterative NMF. For small sample size BB is recommended (default: BB)
  -i I                  No. of NMF iterations. Valid only for -t iNMF (default: 100)
  -outPrefix OUTPREFIX  Output file/s prefix (default: PolyAminer_Out)
  ```
  
## Dependency
1) Data Processing: fastp, bowtie2, samtools, featureCounts/Subread, umi-tools, bedtools  
2) Python modules: pandas, cython, pybedtools, scipy, sklearn, statsmodels, rpy2
3) R libraries: countdata 

Test environment: fastp v0.20.0, bedtools v2.27.1, bowtie2 v2.3.5.1, samtools v1.10, featureCounts v2.03

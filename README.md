# PolyA-miner
Accurate assessment of differential alternative poly-adenylation from 3'Seq data using vector projections and non-negative matrix factorization -Yalamanchili H.K. et al.

## Usage
```
PolyA-miner.py [-h] -d D -o O -pa PA -index INDEX -fasta FASTA -bed BED
                      -c1 C1 [C1 ...] -c2 C2 [C2 ...] [-i I] [-p P]
                      [-expNovel] [-ip IP] [-pa_p PA_P] [-pa_a PA_A]
                      [-pa_m PA_M] [-apa_min APA_MIN] [-gene_min GENE_MIN]
                      [-mdapa MDAPA] [-md MD] [-apaBlock APABLOCK]
                      [-anchor ANCHOR] [-key KEY]
  ```

PolyA-miner.py: the following arguments are required: -d, -o, -pa, -index, -fasta, -bed, -c1, -c2   

## Arguments  
  ```
  -h, --help          show this help message and exit
  -d D                Base directory of input fastq files
  -o O                Out put directory
  -pa PA              PolyA annotations file
  -index INDEX        Reference genome bowtie2 index
  -fasta FASTA        Reference fasta sequence
  -bed BED            Reference genes bed file
  -c1 C1 [C1 ...]     Comma-separated list of condition1 fastq files
  -c2 C2 [C2 ...]     Comma-separated list of condition2 fastq files
  -i I                No. of NMF iterations
  -p P                No. of processors to use
  -expNovel           Explore novel APA sites
  -ip IP              Internal priming window
  -pa_p PA_P          pOverA filter: P
  -pa_a PA_A          pOverA filter: A
  -pa_m PA_M          pOverA filter: M
  -apa_min APA_MIN    Min. proportion per APA
  -gene_min GENE_MIN  Min counts per Gene
  -mdapa MDAPA        Cluster distance for annotated polyA sites
  -md MD              Cluster distance for de-novo polyA sites
  -apaBlock APABLOCK  Block size for annotated polyA sites
  -anchor ANCHOR      Anchor length for de-novo polyA sites
  -key KEY            Output file key
  ```
  
## Dependency
1) Data Processing: fastp, bowtie2, samtools    
2) Python libraries: pandas, cython, pybedtools, scipy, sklearn, statsmodels    

PolyA-miner will check for key requied packages and install them. If any of the installations fail, try installing python3.X-dev library ```sudo apt-get install python3.X-dev```  
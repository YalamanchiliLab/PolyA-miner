# BAM Mode #
# Beta Binomial # # Annotated #
python3 PolyA-miner/PolyA-miner.py -mode bam -o TestOutPut_bam_BB_BAM_DB -pa ReferenceFiles/Test.RefPolyA.bed -index /mnt/local_storage/ylab/Index_Files/Human_hg38v33/bowtie2/GRCh38v33 -fasta ReferenceFiles/Test.Genome.fa -bed ReferenceFiles/Test.Genes.bed -p 20 -pa_p 0.6 -pa_a 3 -pa_m 1 -ip 30 -outPrefix BamBBAno -c1 TestBAM_Data/C1.bam,TestBAM_Data/C2.bam,TestBAM_Data/C3.bam -c2 TestBAM_Data/T1.bam,TestBAM_Data/T2.bam,TestBAM_Data/T3.bam -expNovel 0 -t BB

# Fastq Mode #
# Beta Binomial # # Annotated #
#python3 PolyA-miner/PolyA-miner.py -mode fastq -d TestFastq_Data -o TestOutPut_fastq_BB_DB -pa ReferenceFiles/Test.RefPolyA.bed -index <bowtie2 index> -fasta ReferenceFiles/Test.Genome.fa -bed ReferenceFiles/Test.Genes.bed -p 20 -pa_p 0.6 -pa_a 3 -pa_m 1 -ip 30 -outPrefix FastqBBAno -c1 C1.fastq.gz,C2.fastq.gz,C3.fastq.gz -c2 T1.fastq.gz,T2.fastq.gz,T3.fastq.gz -expNovel 0 -umi 0 -t BB

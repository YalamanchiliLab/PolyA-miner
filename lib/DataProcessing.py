
import concurrent.futures as cf, os, glob, time

def process_rawfastq(d, s, mc):
    if int(mc) > 16:
        mc = '16'
    cmd = 'fastp -i ' + s + ' -a AGATCGGAAGAGC -f 6 -g -l 40 -h ' + d + '/' + d.split('/')[-1] + '_trim-report.html -w ' + mc + ' -Q -o ' + d + '/' + d.split('/')[-1] + '_trim.fastq'
    print(cmd)
    os.system(cmd)
    cmd = 'pigz -p ' + str(mc) + ' ' + d + '/' + d.split('/')[-1] + '_trim.fastq'
    print(cmd)
    os.system(cmd)


def mapping_bowtie2(d, mc, ref_genome):
    cmd = 'bowtie2 -p ' + str(mc) + ' --very-sensitive-local -x ' + ref_genome + ' -U ' + d + '/' + d.split('/')[-1] + '_trim.fastq.gz -S ' + d + '/' + d.split('/')[-1] + '.sam'
    print(cmd)
    os.system(cmd)


def sam2bam(f):
    cmd = 'samtools view -bS ' + f + ' -o ' + f.replace('.sam', '.bam')
    print(cmd)
    os.system(cmd)


def sortbam(f, mc):
    cmd = 'samtools sort -@ ' + str(mc) + ' ' + f + ' -o' + f.replace('.bam', '.sorted.bam')
    print(cmd)
    os.system(cmd)
    #os.system('rm ' + f)


def indexbam(f):
    cmd = 'samtools index ' + f
    print(cmd)
    os.system(cmd)

    
def DataProcessing(baseDir, outDir, np, ref_genome, controls, treated, fkey, lf):
    logfile = open(lf,'a')
    nc = len(controls)
    nt = len(treated)
    samples = controls + treated
    for s in samples:
        s = baseDir + '/' + s
        d = outDir + s.split('/')[-1].replace('.fastq.gz', '')
        os.system('mkdir ' + d)
        process_rawfastq(d, s, np)
        localdate = time.strftime('%a %m/%d/%Y')
        localtime = time.strftime('%H:%M:%S')
        logfile.write('# Finished fastp: ' + localdate + ' at: ' + localtime + 'for ' + d + ' \n')
        os.system('rm *fastp.json')
        mapping_bowtie2(d, np, ref_genome)
        localdate = time.strftime('%a %m/%d/%Y')
        localtime = time.strftime('%H:%M:%S')
        logfile.write('# Finished mapping on ' + localdate + ' at: ' + localtime + 'for ' + d + ' \n')

    files = glob.glob(outDir + '*/' + '*.sam')
    with cf.ProcessPoolExecutor(max_workers=int(np)) as (executor):
        result = list(executor.map(sam2bam, files))
    os.system('rm ' + outDir + '*/' + '*.sam')
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished converting sam 2 bam on ' + localdate + ' at: ' + localtime + ' ##\n')
    files = glob.glob(outDir + '*/' + '*.bam')
    for f in files:
        sortbam(f, np)

    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished sorting bams on ' + localdate + ' at: ' + localtime + ' ##\n')
    files = glob.glob(outDir + '*/' + '*.sorted.bam')
    with cf.ProcessPoolExecutor(max_workers=int(np)) as (executor):
        result = list(executor.map(indexbam, files))
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished indexing bams on ' + localdate + ' at: ' + localtime + ' ##\n')
    os.system('rm ' + outDir + '*/' + '*_trim.fastq.gz')
    os.system('multiqc ' + outDir + ' -o ' + outDir + fkey + '_MultiQCReport')
    logfile.close()
    return (1)


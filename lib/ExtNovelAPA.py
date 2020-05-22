
import numpy as np, pandas as pd
from pybedtools import BedTool
import os, sys, glob, time, uuid, concurrent.futures as cf
from numpy import array, dot, arccos, clip
from numpy.linalg import norm

def chk_strand(row, APAblock):
    if row['strand'] == '+':
        row['start'] = row['start'] - APAblock
    if row['strand'] == '-':
        row['end'] = row['end'] + APAblock
    return (row)

def make_bg(file):
    tfp = str(uuid.uuid4())
    bamf = BedTool(file)
    bgf = bamf.genome_coverage(bg=True, strand='+')
    bgf.saveas(tfp)
    df = pd.read_csv(tfp, sep='\t', header=None, names=['Chr', 'Start', 'End', 'Strand'], index_col=None)
    df['Strand'] = '+'
    df.to_csv(tfp, sep='\t', index=False, header=None)
    tfn = str(uuid.uuid4())
    bgf = bamf.genome_coverage(bg=True, strand='-')
    bgf.saveas(tfn)
    df = pd.read_csv(tfn, sep='\t', header=None, names=['Chr', 'Start', 'End', 'Strand'], index_col=None)
    df['Strand'] = '-'
    df.to_csv(tfn, sep='\t', index=False, header=None)
    os.system('cat  ' + tfp + ' ' + tfn + ' > ' + file.replace('.sorted.bam', '.bg'))
    temp = BedTool(file.replace('.sorted.bam', '.bg'))
    temp = temp.sort()
    temp.saveas(file.replace('.sorted.bam', '.bg'))
    os.system('rm ' + tfp + ' ' + tfn)


def mergeBG(f):
    pars = f.split(',')
    with open(pars[0], 'r') as f:
        lines = f.readlines()
    f.close()
    fw = open(pars[0].replace('.bg', '.bed'), 'w')
    for line in lines:
        data = line.strip().split('\t')
        if 'M' not in data[0] and '_' not in data[0]:
           fw.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3] + '\tN\t' + data[3] + '\n')
    fw.close()
    bed = BedTool(pars[0].replace('.bg', '.bed'))
    bed = bed.sort()
    merged = bed.merge(s=True)
    merged.saveas(pars[0].replace('.bg', '.bed'))


def add_gene_name(a, b):
    os.system('bedtools sort -i ' + a + ' > ' + a.replace('.bed', '.sorted.bed'))
    cmd = 'bedtools closest -a ' + a.replace('.bed', '.sorted.bed') + ' -b ' + b + ' -s -id -D a -t first -k 1 > AltPA.temp.gene.bed'
    os.system(cmd)
    tempdf = pd.read_csv('AltPA.temp.gene.bed', sep='\t', index_col=None, header=None)
    tempdf = tempdf.iloc[:, [0, 1, 2, 5, 9, 12]]
    tempdf.columns = ['Chr', 'Start', 'End', 'Strand', 'Gene', 'Distance']
    tempdf = tempdf[['Chr', 'Start', 'End', 'Gene', 'Distance', 'Strand']]
    os.system('rm AltPA.temp.gene.bed ' + a.replace('.bed', '.sorted.bed'))
    return (tempdf)


def computeA(string, size):
    string = string.strip().upper()
    windows = [string[i:i + size] for i in range(len(string) - (size - 1))]
    for w in windows:
        if w.count('A') >= 12:
            return (0)

    return (1)


def UTR3_counts(f):
    file = f.split(',')[0]
    saf = f.split(',')[1]
    cmd = 'featureCounts -a ' + saf + ' -F SAF -f -O --readExtension5 0 --readExtension3 0 -M -s 0 -T 5 -o ' + file.replace('.sorted.bam', '') + '.UTR3.Counts.txt ' + file
    os.system(cmd)


def ExtNovelAPA(outDir, fkey, ref_polyA, ref_bed, ref_fasta, APAblock, mddb, md, anchor, iplen, np, lf):
    logfile = open(lf,'a')
    df = pd.read_csv(ref_polyA, sep='\t', index_col=None, header=None, names=['chr', 'start', 'end', 'gene', 'PAtype', 'strand'])
    df = df[df['gene'] != 'na']
    df = df.apply(lambda row: chk_strand(row, APAblock), axis=1)
    df.to_csv(outDir + fkey + '.APSitesDB.bed', sep='\t', index=None, header=None)
    bed = BedTool(outDir + fkey + '.APSitesDB.bed')
    bed = bed.sort()
    merged = bed.merge(s=True, c=[4, 5], o='distinct', d=mddb)
    merged.saveas(outDir + fkey + '.APSitesDB.bed')
    df = pd.read_csv(outDir + fkey + '.APSitesDB.bed', sep='\t', index_col=None, header=None, names=['chr', 'start', 'end', 'strand', 'gene', 'PAtype'])
    df = df[['chr', 'start', 'end', 'gene', 'PAtype', 'strand']]
    df = df[~df['chr'].str.contains('_')]
    df = df[~df['gene'].str.contains(',')]
    df.to_csv(outDir + fkey + '.APSitesDB.bed', sep='\t', index=None, header=None)

    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished extracting AP sites from PolyA_DB: ' + localdate + ' at: ' + localtime + ' \n')
    
    # make bged grapg files #
    files = glob.glob(outDir + '*/' + '*.sorted.bam')
    with cf.ProcessPoolExecutor(max_workers=int(np)) as (executor):
        result = list(executor.map(make_bg, files))
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finifhed generating bed graph tracks: ' + localdate + ' at: ' + localtime + ' \n')
    
    # sort and merge bed graph files #
    files = glob.glob(outDir + '*/' + '*.bg')
    passlist = []
    for file in files:
        passlist.append(file)
    with cf.ProcessPoolExecutor(max_workers=int(np)) as (executor):
        result = list(executor.map(mergeBG, passlist))
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished generating bed: ' + localdate + ' at: ' + localtime + ' \n')
    os.system("rm "+outDir + '*/' + '*.bg')
    
    # extract features #
    files = glob.glob(outDir + '*/' + '*.bed')
    file = (' ').join(files)
    os.system('cat ' + file + ' > ' + outDir + fkey + '_Jumbo.bed')
    os.system("rm "+outDir + '*/' + '*.bed')
    df = pd.read_csv(outDir + fkey + '_Jumbo.bed', sep='\t', header=None, names=['Chr', 'Start', 'End', 'Strand'], index_col=None)
    df['N1'] = 'N1'
    df['N2'] = 'N2'
    df = df[['Chr', 'Start', 'End', 'N1', 'N2', 'Strand']]
    df.to_csv(outDir + fkey + '_Jumbo.bed', sep='\t', index=False, header=None)
    pooled = BedTool(outDir + fkey + '_Jumbo.bed')
    pooled = pooled.sort()
    pooled.saveas(outDir + fkey + '_Jumbo.bed')
    merged_bg = pooled.merge(s=True, d=md)
    merged_bg.saveas(outDir + fkey + '_merged.features.bed')
    df = pd.read_csv(outDir + fkey + '_merged.features.bed', sep='\t', header=None, names=['Chr', 'Start', 'End', 'Strand'], index_col=None)
    nrow = df.shape[0]
    names = []
    for n in range(1, nrow + 1):
        names.append('3UTR_' + str(n))
    df['GeneID'] = names
    df = df[['Chr', 'Start', 'End', 'GeneID', 'GeneID', 'Strand']]
    df.to_csv(outDir + fkey + '_denovoAPAsites.bed', sep='\t', index=False, header=None)
    os.system('rm ' + outDir + fkey + '_Jumbo.bed ' + outDir + fkey + '_merged.features.bed')

    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finifhed generating features bed file on ' + localdate + ' at: ' + localtime + ' ##\n')
    
    # add geen names #
    df = add_gene_name(outDir + fkey + '_denovoAPAsites.bed', ref_bed)
    df = df[df.Distance >= -16000]
    df.to_csv(outDir + fkey + '_denovoAPAsites.bed', sep='\t', index=False, header=None)
    
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished mapping features to Genes ' + localdate + ' at: ' + localtime + ' ##\n')
    
    # skip annotated sites from miss-priming screen #
    os.system('bedtools intersect -a ' + outDir + fkey + '_denovoAPAsites.bed -b ' + outDir + fkey + '.APSitesDB.bed -s -wo > ' + outDir + fkey + '_AnnotatedAPA.bed')
    with open(outDir + fkey + '_AnnotatedAPA.bed', 'r') as fr:
        lines = fr.readlines()
    fr.close()
    fw = open(outDir + fkey + '_AnnotatedAPA.bed', 'w')
    for line in lines:
        data = line.strip().split('\t')
        if int(data[12]) >= anchor:
            fw.write(('\t').join(data[0:6]) + '\n')
    fw.close()
    df = pd.read_csv(outDir + fkey + '_AnnotatedAPA.bed', sep='\t', index_col=None, names=['Chr', 'Start', 'End', 'GeneID', 'TID', 'Strand'])
    df = df.drop_duplicates(subset=['Chr', 'Start', 'End'], keep='first', inplace=False)
    df.to_csv(outDir + fkey + '_AnnotatedAPA.bed', sep='\t', index=False, header=None)
    
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished compiling annotated features : ' + localdate + ' at: ' + localtime + ' \n')
    
    # mark novel putative APA sites #
    os.system('bedtools intersect -a ' + outDir + fkey + '_denovoAPAsites.bed -b ' + outDir + fkey + '_AnnotatedAPA.bed -s -v > ' + outDir + fkey + '_novelAPA.bed')
    with open(outDir + fkey + '_novelAPA.bed', 'r') as fr:
        lines = fr.readlines()
    fr.close()
    fw = open(outDir + fkey + ('_novelAPA.bed').replace('.bed', '.seq.bed'), 'w')
    for line in lines:
        data = line.strip().split('\t')
        if data[5] == '+' and data[0] != 'chrM':
            fw.write(data[0] + '\t' + str(int(data[2]) - 10) + '\t' + str(int(data[2]) + iplen) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')
        if data[5] == '-' and data[0] != 'chrM':
            fw.write(data[0] + '\t' + str(int(data[1]) - iplen) + '\t' + str(int(data[1]) + 10) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')

    fw.close()
    os.system('bedtools getfasta -fo ' + outDir + fkey + '_novelAPA.bed'.replace('.bed', '.seq.fasta') + ' -name -s -fi ' + ref_fasta + ' -bed ' + outDir + fkey + '_novelAPA.bed'.replace('.bed', '.seq.bed'))
    with open(outDir + fkey + '_novelAPA.bed'.replace('.bed', '.seq.fasta'), 'r') as fr:
        lines = fr.readlines()
    fr.close()
    fw = open(outDir + fkey + '_novelAPA.bed'.replace('.bed', '.cleaned.bed'), 'w')
    for i in range(0, len(lines), 2):
        try:
            status = computeA(lines[i + 1], 20)
            if status == 1:
                fw.write(lines[i].replace('>', '').replace('_', '\t'))
        except:
            pass
    fw.close()

    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished cleaning Internally primed non annotated features : ' + localdate + ' at: ' + localtime + ' \n')
   
    # combine annotateed and novel APA sites #
    os.system('cat ' + outDir + fkey + '_AnnotatedAPA.bed ' + outDir + fkey + '_novelAPA.bed'.replace('.bed', '.cleaned.bed') + ' >' + outDir + fkey + '_denovoAPAsites.bed')
    pooled = BedTool(outDir + fkey + '_denovoAPAsites.bed')
    pooled = pooled.sort()
    pooled.saveas(outDir + fkey + '_denovoAPAsites.bed')
    os.system('rm ' + outDir + fkey + '_novelAPA.* ' + outDir + fkey + '_AnnotatedAPA.bed ')
    
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished compiling annotated and un-annotated AP sites: ' + localdate + ' at: ' + localtime + ' \n')
    
    df = pd.read_csv(outDir + fkey + '_denovoAPAsites.bed', sep='\t', index_col=None, header=None, names=['Chr', 'Start', 'End', 'gene_id', 'Distance', 'Strand'])
    df = df.drop_duplicates(subset=['Chr', 'Start', 'End', 'Strand'])
    df['GeneID'] = df['gene_id'] + '@' + df['Chr'] + '_' + df['Start'].apply(str) + '_' + df['End'].apply(str) + '_' + df['Strand']
    df = df[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
    df.to_csv(outDir + fkey + '_denovoAPAsites.saf', sep='\t', index=False, header=None)
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished making SAF file: ' + localdate + ' at: ' + localtime + ' \n')
    logfile.close()
    return (1)


import pandas as pd, time
from pybedtools import BedTool

def chk_strand(row, APAblock):
    if row['strand'] == '+':
        row['start'] = row['start'] - APAblock
    if row['strand'] == '-':
        row['end'] = row['end'] + APAblock
    return row


def ExtAnnotatedAPA(outDir, fkey, ref_polyA, APAblock, mddb, lf):
    logfile=open(lf,"a")
    df = pd.read_csv(ref_polyA, sep='\t', index_col=None, header=None, names=['chr', 'start', 'end', 'gene', 'PAtype', 'strand'])
    df = df[df['gene'] != 'na']
    df = df[df['PAtype'] != 'Intron']
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
    df = pd.read_csv(outDir + fkey + '.APSitesDB.bed', sep='\t', index_col=None, names=['Chr', 'Start', 'End', 'Gene', 'APA', 'Strand'])
    df['GeneID'] = df['Gene'] + '@' + df['Chr'] + '_' + df['Start'].apply(str) + '_' + df['End'].apply(str) + '_' + df['Strand']
    df = df[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
    df.to_csv(outDir + fkey + '_APA.saf', sep='\t', index=False, header=None)
    logfile.close()
    return (1)


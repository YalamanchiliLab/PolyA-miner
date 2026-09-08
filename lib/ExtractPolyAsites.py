# Extract PolyA sites from DB and de-novo #

import os, sys, glob, time, uuid, math
import concurrent.futures as cf
import pandas as pd,subprocess
from pybedtools import BedTool
import pybedtools as pb


def chk_strand(row, APAblock):
	if row['strand'] == '+':
		row['start'] = row['start'] - APAblock
	if row['strand'] == '-':
		row['end'] = row['end'] + APAblock
	return (row)


def addAttribute(df,alist):
	for element in alist:	
		df[element] = df['attributes'].str.extract(rf'{element}\s+"(.*?)"')
	# Drop the original attributes column
	df.drop(columns=['attributes'], inplace=True)
	return(df)


def makeGeneBed(args_gtf,args_fkey,args_o):
	args_o=args_o.strip("/")+"/"+args_fkey
	gtf = pb.BedTool(args_gtf)
	genes = gtf.filter(lambda x: x[2] == 'gene' and 'chr' in x[0]).saveas()
	df = genes.to_dataframe()
	#print(df)
	# Extract gene_id and gene_name from attributes column
	df=addAttribute(df,['gene_id','gene_name'])
	df['start'] = df['start'] - 1
	df[['seqname', 'start', 'end', 'gene_id', 'gene_name', 'strand']].to_csv(args_o+'.Genes.bed', sep='\t', index=False,header=False)
	return(args_o+'.Genes.bed')


def get_intron_regions(genes_df, cds_df, utr5_df, utr3_df):
	# Convert df_cds, df_utr5, and df_utr3 to BedTool objects
	genes_bed = pb.BedTool.from_dataframe(genes_df)
	cds_bed = pb.BedTool.from_dataframe(cds_df)
	utr5_bed = pb.BedTool.from_dataframe(utr5_df)
	utr3_bed = pb.BedTool.from_dataframe(utr3_df)
	intron = pb.BedTool(genes_bed).subtract(cds_bed, s=True).subtract(utr5_bed, s=True).subtract(utr3_bed, s=True)
	return (intron)


def getTables(gtf):
	cds = gtf.filter(lambda x: x[2] == 'CDS').saveas()
	utr = gtf.filter(lambda x: x[2] == 'UTR').saveas()
	stop_codons = gtf.filter(lambda x: x[2] == 'stop_codon' and 'chr' in x[0]).saveas()
	start_codons = gtf.filter(lambda x: x[2] == 'start_codon' and 'chr' in x[0]).saveas()
	
	# Label UTRs as 5'UTR and 3'UTR
	df_cds = cds.to_dataframe()
	df_cds=addAttribute(df_cds,['gene_id'])
	df_cds['exon_id'] = "CDS"
	df_cds=df_cds[['seqname', 'start', 'end', 'gene_id', 'exon_id', 'strand']]


	df_utr = utr.to_dataframe()
	df_utr=addAttribute(df_utr,['gene_id','transcript_id'])
	df_utr=df_utr[['seqname', 'start', 'end', 'transcript_id','gene_id','strand']]


	df_stop = stop_codons.to_dataframe()
	df_stop=addAttribute(df_stop,['gene_id','transcript_id'])
	df_stop=df_stop[['seqname', 'start', 'end', 'transcript_id','gene_id','strand']]


	df_start = start_codons.to_dataframe()
	df_start=addAttribute(df_start,['gene_id','transcript_id'])
	df_start=df_start[['seqname', 'start', 'end', 'transcript_id','gene_id','strand']]

	# Merge UTR dataframe with both start and stop codons
	merged_utr = pd.merge(df_utr, df_stop, on=['transcript_id', 'gene_id', 'seqname', 'strand'], how='left', suffixes=('', '_stop'))
	merged_utr = pd.merge(merged_utr, df_start, on=['transcript_id', 'gene_id', 'seqname', 'strand'], how='left', suffixes=('', '_start'))

	#merged_utr.to_csv("UTR_Merged.bed", sep='\t', index=False, header=True)
	

	def label_utr(row):
		if row['strand'] == '+':
			if pd.notna(row['start_start']) and row['end'] <= row['start_start']:
				return '5UTR'
			if pd.notna(row['start_stop']) and row['start'] >= row['start_stop']:
				return '3UTR'
		elif row['strand'] == '-':
			if pd.notna(row['end_stop']) and row['end'] <= row['end_stop']:
				return '3UTR'
			if pd.notna(row['start_start']) and row['start'] >= row['start_start']:
				return '5UTR'
		return 'unknown'

	# Label the UTRs
	merged_utr['UTR_type'] = merged_utr.apply(label_utr, axis=1)



	# Separate the UTRs based on their type
	utr5_df = merged_utr[merged_utr['UTR_type'] == '5UTR'].copy()[['seqname', 'start', 'end', 'transcript_id', 'gene_id', 'strand']]
	utr3_df = merged_utr[merged_utr['UTR_type'] == '3UTR'].copy()[['seqname', 'start', 'end', 'transcript_id', 'gene_id', 'strand']]


	#df_cds = cds.to_dataframe()

	utr5_bed=pb.BedTool.from_dataframe(utr5_df)
	utr5_bed = utr5_bed.sort()
	utr5_merge=utr5_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
	utr5_df=utr5_merge.to_dataframe(disable_auto_names=True,header=None)
	utr5_df.columns=['chr','start','end','GeneID','Feature','strand']
	utr5_df['Feature']="UTR5"

	utr3_bed=pb.BedTool.from_dataframe(utr3_df)
	utr3_bed = utr3_bed.sort()
	utr3_merge=utr3_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
	utr3_df=utr3_merge.to_dataframe(disable_auto_names=True,header=None)
	utr3_df.columns=['chr','start','end','GeneID','Feature','strand']
	utr3_df['Feature']="UTR3"

	
	cds_bed = pb.BedTool.from_dataframe(df_cds)
	cds_bed = cds_bed.sort()
	cds_merge=cds_bed.merge(s=True, c=[4,5,6], o='distinct', d=1)
	#cds_merge=cds_merge.subtract(utr5_merge, s=True)
	#cds_merge=cds_merge.subtract(utr3_merge, s=True)
	cds_df=cds_merge.to_dataframe(disable_auto_names=True,header=None)
	cds_df.columns=['chr','start','end','GeneID','Feature','strand']
	cds_df['Feature']="CDS"


	utr5_merge=utr5_merge.subtract(cds_merge, s=True)
	utr5_df=utr5_merge.to_dataframe(disable_auto_names=True,header=None)
	utr5_df.columns=['chr','start','end','GeneID','Feature','strand']
	utr5_df['Feature']="UTR5"

	utr3_merge=utr3_merge.subtract(cds_merge, s=True)
	utr3_df=utr3_merge.to_dataframe(disable_auto_names=True,header=None)
	utr3_df.columns=['chr','start','end','GeneID','Feature','strand']
	utr3_df['Feature']="UTR3"

	
	return (cds_df, utr5_df, utr3_df)


def mapAPA2Features(APAfile,Features):
	PA_df=pd.read_csv(APAfile,sep="\t",header=None,index_col=None)
	PA_df.columns=['Chr', 'Start', 'End', 'GeneID', 'feature', 'Strand']
	PA=pb.BedTool.from_dataframe(PA_df)

	FT = pb.BedTool(Features)
	PA_FT_bed = PA.intersect(FT, nonamecheck=True, s=True, wao=True)
	PA_FT_df=PA_FT_bed.to_dataframe(disable_auto_names=True, header=None)
	PA_FT_df.columns=['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand','Chr_2', 'Start_2', 'End_2', 'GeneID_2', 'APAID_2', 'Strand_2','overlap']	
	# get no-map df #
	PA_FT_df_nomap=PA_FT_df[PA_FT_df['Chr_2']=="."].copy()

	PA_FT_df=PA_FT_df[(((PA_FT_df['Strand'] =="+") & (PA_FT_df['APAID_2'] !="UTR3")) & (PA_FT_df['End'] >= PA_FT_df['Start_2']) & (PA_FT_df['End'] <= PA_FT_df['End_2'])) | (((PA_FT_df['Strand'] =="+") & (PA_FT_df['APAID_2'] =="UTR3")) & (PA_FT_df['End'] >= PA_FT_df['Start_2'])) | (((PA_FT_df['Strand'] =="-") & (PA_FT_df['APAID_2'] !="UTR3")) & (PA_FT_df['Start'] <= PA_FT_df['End_2']) & (PA_FT_df['Start'] >= PA_FT_df['Start_2'])) | (((PA_FT_df['Strand'] =="-") & (PA_FT_df['APAID_2'] =="UTR3")) & (PA_FT_df['Start'] <= PA_FT_df['End_2']))]
	PA_FT_df['APAID']=PA_FT_df['APAID_2']
	PA_FT_df=PA_FT_df[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
		
	PA_FT_bed=pb.BedTool.from_dataframe(PA_FT_df)
	PA_FT_bed = PA_FT_bed.sort()
	PA_FT_bed=PA_FT_bed.merge(c=[4,5,6],s=True,o="distinct")
	PA_FT_df=PA_FT_bed.to_dataframe(disable_auto_names=True, header=None)
	PA_FT_df.columns=['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']

	# process nomap df#
	PA_FT_df_nomap=PA_FT_df_nomap[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
	PA_FT_df_nomap['APAID']="UN"
		
	PA_nomap_FT=pb.BedTool.from_dataframe(PA_FT_df_nomap).intersect(FT, nonamecheck=True, s=False, wao=True)
	PA_nomap_FT_df=PA_nomap_FT.to_dataframe(disable_auto_names=True, header=None)
	if len(PA_nomap_FT_df) >= 1:
		PA_nomap_FT_df.columns=['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand','Chr_2', 'Start_2', 'End_2', 'GeneID_2', 'APAID_2', 'Strand_2','overlap']
		UTR3_EX_df=PA_nomap_FT_df[PA_nomap_FT_df['Start_2']==-1].copy()
		UTR3_EX_df=UTR3_EX_df[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
		UTR3_EX_df['APAID']="UTR3-EX"
		PAmap_df=pd.concat([PA_FT_df,UTR3_EX_df], axis=0,ignore_index=True)
	else:
		PAmap_df=PA_FT_df.copy()
		#print(PA_FT_df)

	PA_df=PA_df.merge(PAmap_df,on=['Chr', 'Start', 'End', 'GeneID', 'Strand'],how='left').fillna("UN")
	PA_df=PA_df.drop(columns=['feature'])
	PA_df=PA_df[['Chr', 'Start', 'End', 'GeneID', 'APAID', 'Strand']]
	PA_df.to_csv(APAfile, sep='\t', index=False, header=False)
	
def MarkFeatures(args_gtf,args_fkey,args_o,APAfile):
	args_o=args_o.strip("/")+"/"+args_fkey
	gtf = pb.BedTool(args_gtf)
	gtf = gtf.filter(lambda x: 'chr' in x[0]).saveas()
	
	# Get Gene bed #
	genes = gtf.filter(lambda x: x[2] == 'gene').saveas()
	gene_df = genes.to_dataframe()

	# Extract gene_id and gene_name from attributes column
	gene_df['gene_id'] = gene_df['attributes'].str.extract(r'gene_id\s+"(.*?)"')
	gene_df['gene_name'] = gene_df['attributes'].str.extract(r'gene_name\s+"(.*?)"')
	
	# Drop the original attributes column
	gene_df.drop(columns=['attributes'], inplace=True)
	gene_df['start'] = gene_df['start'] - 1
	gene_df[['seqname', 'start', 'end', 'gene_id', 'gene_name', 'strand']].to_csv(args_o+'.Genes.bed', sep='\t', index=False,header=False)
	
	# Get features #
	cds_df,utr5_df,utr3_df=getTables(gtf)

	#cds_df.to_csv(args_o+".CDS.bed",sep="\t",header=False,index=False)
	#utr5_df.to_csv(args_o+".5UTR.bed",sep="\t",header=False,index=False)
	#utr3_df.to_csv(args_o+".3UTR.bed",sep="\t",header=False,index=False)

	intron_regions = get_intron_regions(gene_df[['seqname', 'start', 'end', 'gene_id', 'gene_name', 'strand']], cds_df,utr5_df,utr3_df)
	intr_df = intron_regions.to_dataframe()
	intr_df.columns=['chr','start','end','GeneID','Feature','strand']
	intr_df['Feature']='Intron'
	#intr_df.to_csv(args_o+".Introns.bed",sep="\t",header=False,index=False)

	feature_df=pd.concat([intr_df,cds_df,utr5_df,utr3_df], axis=0,ignore_index=True)
	feature_df.to_csv(args_o+".Features.bed",sep="\t",header=False,index=False)

	mapAPA2Features(APAfile,args_o+".Features.bed")

	return(1)


def ExtAPAfromPolyA_DB(outDir, fkey, ref_polyA, APAblock, mddb,logfile):
	df = pd.read_csv(ref_polyA, sep='\t', index_col=None, header=None, names=['chr', 'start', 'end', 'gene', 'PAtype', 'strand'])
	df = df[df['gene'] != 'na']
	#df = df[df['PAtype'] != 'Intron']
	#df = df.apply(lambda row: chk_strand(row, APAblock), axis=1) # Slow 

	# Adjust 'start' for rows where 'strand' is '+'
	df.loc[df['strand'] == '+', 'start'] = df['start'] - APAblock
	# Adjust 'end' for rows where 'strand' is '-'
	df.loc[df['strand'] == '-', 'end'] = df['end'] + APAblock

	df.to_csv(outDir + fkey + '.APSitesDB.bed', sep='\t', index=None, header=None)
	bed = BedTool(outDir + fkey + '.APSitesDB.bed')
	bed = bed.sort()
	#merged = bed.merge(s=True, c=[4, 5], o='distinct', d=mddb) # Worked with bedtools 2.26 but failed with 2.29
	merged = bed.merge(s=True, c=[4, 5, 6], o='distinct', d=mddb) # Worked with bedtools 2.29
	merged.saveas(outDir + fkey + '.APSitesDB.bed')
	#df = pd.read_csv(outDir + fkey + '.APSitesDB.bed', sep='\t', index_col=None, header=None, names=['chr', 'start', 'end', 'strand', 'gene', 'PAtype']) # Worked with bedtools 2.26 but failed with 2.29
	df = pd.read_csv(outDir + fkey + '.APSitesDB.bed', sep='\t', index_col=None, header=None, names=['chr', 'start', 'end', 'gene', 'PAtype', 'strand']) # Worked with bedtools 2.29
	df = df[['chr', 'start', 'end', 'gene', 'PAtype', 'strand']]
	print(df["chr"].dtype)
	print(df["chr"].head())
	print(df["start"].dtype)
	print(df["start"].head())
	
	df = df[~df['chr'].str.contains('_')]
	df = df[~df['gene'].str.contains(',')]
	df.to_csv(outDir + fkey + '.APSitesDB.bed', sep='\t', index=None, header=None)
	
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished extracting AP sites from PolyA_DB: ' + localdate + ' at: ' + localtime + ' \n')
	return(1)

def makeSAF(outDir, fkey):
	df = pd.read_csv(outDir + fkey + '.APSitesDB.bed', sep='\t', index_col=None, names=['Chr', 'Start', 'End', 'Gene', 'APA', 'Strand'])
	df['GeneID'] = df['Gene'] + '@' + df['Chr'] + '_' + df['Start'].apply(str) + '_' + df['End'].apply(str) + '_' + df['Strand']
	df = df[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
	df.to_csv(outDir + fkey + '_APA.saf', sep='\t', index=False, header=None)
	return (1)

def make_bg(files):
	tfp = files[1]+str(uuid.uuid4())
	bamf = BedTool(files[0])
	bgf = bamf.genome_coverage(bg=True, strand='+')
	#bgf.saveas(tfp)
	#df = pd.read_csv(tfp, sep='\t', header=None, names=['Chr', 'Start', 'End', 'Strand'], index_col=None)
	#df['Strand'] = '+'
	#df.to_csv(tfp, sep='\t', index=False, header=None)
	bgf_df=bgf.to_dataframe()
	bgf_df.columns=['Chr', 'Start', 'End', 'Strand']
	bgf_df=bgf_df[bgf_df['Strand']>3]
	bgf_df['Strand'] = '+'
	bgf_df.to_csv(tfp, sep='\t', index=False, header=None)

	tfn = files[1]+str(uuid.uuid4())
	bgf = bamf.genome_coverage(bg=True, strand='-')
	#bgf.saveas(tfn)
	#df = pd.read_csv(tfn, sep='\t', header=None, names=['Chr', 'Start', 'End', 'Strand'], index_col=None)
	#df['Strand'] = '-'
	#df.to_csv(tfn, sep='\t', index=False, header=None)
	bgf_df=bgf.to_dataframe()
	bgf_df.columns=['Chr', 'Start', 'End', 'Strand']
	bgf_df=bgf_df[bgf_df['Strand']>3]
	bgf_df['Strand'] = '-'
	bgf_df.to_csv(tfn, sep='\t', index=False, header=None)

	# move temp files to strand specific BG #
	os.system('mv  ' + tfp +' '+files[1]+files[0].split("/")[-1].replace('.bam', '_+_.bg'))
	temp_p = BedTool(files[1]+files[0].split("/")[-1].replace('.bam', '_+_.bg'))
	temp_p = temp_p.sort()
	temp_p.saveas(files[1]+files[0].split("/")[-1].replace('.bam', '_+_.bg'))
	
	os.system('mv  ' + tfn +' '+files[1]+files[0].split("/")[-1].replace('.bam', '_-_.bg'))
	temp_n = BedTool(files[1]+files[0].split("/")[-1].replace('.bam', '_-_.bg'))
	temp_n = temp_n.sort()
	temp_n.saveas(files[1]+files[0].split("/")[-1].replace('.bam', '_-_.bg'))
	return(1)
	
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
	#merged = bed.merge(s=True) # Worked with bedtools 2.26 but failed with 2.29
	merged = bed.merge(s=True,c=[6], o='distinct') # For bedtools 2.29
	merged.saveas(pars[0].replace('.bg', '.bed'))
	return(1)

def add_gene_name(a, b):

	a_bed = BedTool(a)
	# Sort the BedTool object
	sorted_bed = a_bed.sort()
	sorted_bed.saveas(a)

	b_bed = BedTool(b)
	# Sort the BedTool object
	sorted_bed = b_bed.sort()
	sorted_bed.saveas(b)

	# Perform the closest operation with the specified options
	closest_result = a_bed.closest(b_bed, nonamecheck=True, s=True, id=True, D='a', t='first', k=1)
	closest_result.saveas('AltPA.temp.gene.bed')

	tempdf = pd.read_csv('AltPA.temp.gene.bed', sep='\t', index_col=None, header=None)
	tempdf = tempdf.iloc[:, [0, 1, 2, 5, 9, 12]]
	tempdf.columns = ['Chr', 'Start', 'End', 'Strand', 'Gene', 'Distance']
	tempdf = tempdf[['Chr', 'Start', 'End', 'Gene', 'Distance', 'Strand']]

	#os.system('rm AltPA.temp.gene.bed ' + a.replace('.bed', '.sorted.bed'))
	os.system('rm AltPA.temp.gene.bed')
	return (tempdf)

def computeA(string, size,prpA):
	string = string.strip().upper()
	windows = [string[i:i + size] for i in range(len(string) - (size - 1))]
	for w in windows:
		if w.count('A')/float(size) >= prpA:
			return (0)
	return (1)

def ExtNovelAPA(outDir, fkey, ref_bed, ref_fasta, md, anchor, iplen, iplen_up,novelD, prpA, npc, mode, c1,c2,logfile,gtf,args_ignore):	
	# Make bged grapg files #
	if mode =="fastq":
		files = glob.glob(outDir + '*.bam')
		for i in range(0,len(files)):
			files[i]=[files[i],outDir]
	samples=c1+c2
	beds_p=[]
	beds_n=[]
	if mode == "bam":
		for i in range(0,len(samples)):
			beds_p.append(outDir+samples[i].rsplit("/",1)[1].replace(".bam","_+_.bed"))
			beds_n.append(outDir+samples[i].rsplit("/",1)[1].replace(".bam","_-_.bed"))
			samples[i]=[samples[i],outDir]
		files=samples

	with cf.ProcessPoolExecutor(max_workers=int(npc)) as (executor):
		result = list(executor.map(make_bg, files))
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished generating bed graph tracks: ' + localdate + ' at: ' + localtime + ' \n')
	
	# Sort and merge bed graph files #
	files = glob.glob(outDir + '*.bg')
	passlist = []
	for file in files:
		passlist.append(file)
	with cf.ProcessPoolExecutor(max_workers=int(npc)) as (executor):
		result = list(executor.map(mergeBG, passlist))
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished generating bed: ' + localdate + ' at: ' + localtime + ' \n')

	# Condition 1 samples Extract features #
	x_p=BedTool()
	C1df_p = x_p.multi_intersect(i=beds_p[0:len(c1)],header=True).to_dataframe(disable_auto_names=True)	
	#print(beds_p)
	#print(beds_p[0:len(c1)])
	#print(C1df_p.columns)
	#print(C1df_p.head)

	C1df_p=C1df_p[C1df_p['num']>=math.ceil(0.6*len(c1))]
	C1df_p=C1df_p[['chrom','start','end']]
	x_n=BedTool()
	C1df_n = x_n.multi_intersect(i=beds_n[0:len(c1)],header=True).to_dataframe(disable_auto_names=True)	
	C1df_n=C1df_n[C1df_n['num']>=math.ceil(0.6*len(c1))]
	C1df_n=C1df_n[['chrom','start','end']]
	
	# Condition 2 samples Extract features #
	y_p=BedTool()
	C2df_p = y_p.multi_intersect(i=beds_p[len(c1):len(c1+c2)],header=True).to_dataframe(disable_auto_names=True)	
	#print(C2df_p.head)
	#print(C2df_p.columns)
	#print(beds_p[len(c1):len(c1+c2)])
	#C2df_p.to_csv("OnlyPCATest.txt",sep="\t",header=True)
	C2df_p=C2df_p[C2df_p['num']>=math.ceil(0.6*len(c2))]
	C2df_p=C2df_p[['chrom','start','end']]
	y_n=BedTool()
	C2df_n = y_n.multi_intersect(i=beds_n[len(c1):len(c1+c2)],header=True).to_dataframe(disable_auto_names=True)
	C2df_n=C2df_n[C2df_n['num']>=math.ceil(0.6*len(c2))]
	C2df_n=C2df_n[['chrom','start','end']]

	# merge C1 and C1 #
	#C1df_p=C1df_p.append(C2df_p, ignore_index=True)
	C1df_p = pd.concat([C1df_p, C2df_p], ignore_index=True)
	C1df_p = C1df_p[~C1df_p['chrom'].str.contains('_|\\.')]
	C1df_p['strand']='+'
	
	#C1df_n=C1df_n.append(C2df_n, ignore_index=True)
	C1df_n = pd.concat([C1df_n, C2df_n], ignore_index=True)
	C1df_n = C1df_n[~C1df_n['chrom'].str.contains('_|\\.')]
	C1df_n['strand']='-'

	#Trim_Features_df=C1df_p.append(C1df_n, ignore_index=True)
	Trim_Features_df = pd.concat([C1df_p, C1df_n], ignore_index=True)

	Trim_Features_df['N1']='N1'
	Trim_Features_df['N2']='N2'
	Trim_Features_df=Trim_Features_df[['chrom','start','end','N1','N2','strand']]
	Trim_Features=BedTool.from_dataframe(Trim_Features_df)
	Trim_Features = Trim_Features.sort()
	Trim_Features.saveas(outDir + fkey + '_Jumbo.bed')

	# Merge #
	#merged_bg = pooled.merge(s=True, d=md) # Worked with bedtools 2.26 but failed with 2.29
	merged_bg = Trim_Features.merge(s=True,c=[6], o='distinct') # For bedtools 2.29
	merged_bg.saveas(outDir + fkey + '_merged.features.bed')

	df = pd.read_csv(outDir + fkey + '_merged.features.bed', sep='\t', header=None, names=['Chr', 'Start', 'End', 'Strand'], index_col=None)
	nrow = df.shape[0]
	names = []
	for n in range(1, nrow + 1):
		names.append('3UTR_' + str(n))
	df['GeneID'] = names
	df['width'] = df['End'] - df['Start']
	df=df[df['width']>30]
	df = df[['Chr', 'Start', 'End', 'GeneID', 'GeneID', 'Strand']]
	df.to_csv(outDir + fkey + '_denovoAPAsites.bed', sep='\t', index=False, header=None)
	os.system('rm ' + outDir + fkey + '_Jumbo.bed ' + outDir + fkey + '_merged.features.bed')
	os.system("rm "+" ".join(beds_p+beds_n))
	os.system("rm "+outDir+"*.bg")

	
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished generating features bed file on ' + localdate + ' at: ' + localtime + ' ##\n')

	# add gene names #
	df = add_gene_name(outDir + fkey + '_denovoAPAsites.bed', ref_bed)
	#df = df[df.Distance >= -16000]
	df = df[df.Distance >= (novelD*-1)]
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

	
	########## skip this if you donot want to filter annotated sites #
	'''
	with open(outDir + fkey + '_AnnotatedAPA.bed', 'r') as fr:
		lines = fr.readlines()
	fr.close()
	fw = open(outDir + fkey + ('_AnnotatedAPA.bed').replace('.bed', '.seq.bed'), 'w')
	for line in lines:
		data = line.strip().split('\t')
		if data[5] == '+' and data[0] != 'chrM':
			if ((int(data[2]) - 25) > 0):
				fw.write(data[0] + '\t' + str(int(data[2]) - 25) + '\t' + str(int(data[2]) + 25) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')
			else:
				fw.write(data[0] + '\t' + str(int(data[2])) + '\t' + str(int(data[2]) + 25) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')
		if data[5] == '-' and data[0] != 'chrM':
			if ((int(data[1]) - 25) > 0):
				fw.write(data[0] + '\t' + str(int(data[1]) - 25) + '\t' + str(int(data[1]) + 25) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')
			else:
				fw.write(data[0] + '\t' + str(int(data[1])) + '\t' + str(int(data[1]) + 25) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')
	fw.close()
	
	os.system('bedtools getfasta -fo ' + outDir + fkey + '_AnnotatedAPA.bed'.replace('.bed', '.seq.fasta') + ' -nameOnly -s -fi ' + ref_fasta + ' -bed ' + outDir + fkey + '_AnnotatedAPA.bed'.replace('.bed', '.seq.bed'))
	# nameOnly update for recent bedtools #
	with open(outDir + fkey + '_AnnotatedAPA.bed'.replace('.bed', '.seq.fasta'), 'r') as fr:
		lines = fr.readlines()
	fr.close()
	fw = open(outDir + fkey + '_AnnotatedAPA.bed'.replace('.bed', '.cleaned.bed'), 'w')
	for i in range(0, len(lines), 2):
		try:
			status = computeA(lines[i + 1], 20, 0.80)
			if status == 1:
				# fw.write(lines[i].replace('>', '').replace('_', '\t')) 
				fw.write(lines[i].rsplit("(",1)[0].replace('>', '').replace('_', '\t')+"\n") # for extra strand info bedtools v.2.27
		except:
			pass
	fw.close()

	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished cleaning Internally primed annotated features : ' + localdate + ' at: ' + localtime + ' \n')
	'''
	##########

	os.system("cp "+outDir + fkey + '_AnnotatedAPA.bed ' +outDir + fkey + ('_AnnotatedAPA.bed').replace('.bed', '.cleaned.bed'))
	

	# Mark novel putative APA sites #
	os.system('bedtools intersect -a ' + outDir + fkey + '_denovoAPAsites.bed -b ' + outDir + fkey + '_AnnotatedAPA.bed -s -v > ' + outDir + fkey + '_novelAPA.bed')
	with open(outDir + fkey + '_novelAPA.bed', 'r') as fr:
		lines = fr.readlines()
	fr.close()
	fw = open(outDir + fkey + ('_novelAPA.bed').replace('.bed', '.seq.bed'), 'w')
	for line in lines:
		data = line.strip().split('\t')
		if data[5] == '+' and data[0] != 'chrM':
			if ((int(data[2]) - iplen_up) > 0):
				fw.write(data[0] + '\t' + str(int(data[2]) - iplen_up) + '\t' + str(int(data[2]) + iplen) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')
			else:
				fw.write(data[0] + '\t' + str(int(data[2])) + '\t' + str(int(data[2]) + iplen) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')
		if data[5] == '-' and data[0] != 'chrM':
			if ((int(data[1]) - iplen) > 0):
				fw.write(data[0] + '\t' + str(int(data[1]) - iplen) + '\t' + str(int(data[1]) + iplen_up) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')
			else:
				fw.write(data[0] + '\t' + str(int(data[1])) + '\t' + str(int(data[1]) + iplen_up) + '\t' + data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3].replace('_', '-') + '_' + data[4].replace('_', '-') + '_' + data[5] + '\tN\t' + data[5] + '\n')

	fw.close()
	os.system('bedtools getfasta -fo ' + outDir + fkey + '_novelAPA.bed'.replace('.bed', '.seq.fasta') + ' -nameOnly -s -fi ' + ref_fasta + ' -bed ' + outDir + fkey + '_novelAPA.bed'.replace('.bed', '.seq.bed'))
	with open(outDir + fkey + '_novelAPA.bed'.replace('.bed', '.seq.fasta'), 'r') as fr:
		lines = fr.readlines()
	fr.close()
	fw = open(outDir + fkey + '_novelAPA.bed'.replace('.bed', '.cleaned.bed'), 'w')
	for i in range(0, len(lines), 2):
		try:
			status = computeA(lines[i + 1], 20, prpA)
			if status == 1:
				# fw.write(lines[i].replace('>', '').replace('_', '\t')) 
				fw.write(lines[i].rsplit("(",1)[0].replace('>', '').replace('_', '\t')+"\n") # for extra strand info bedtools v.2.27
		except:
			pass
	fw.close()

	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished cleaning Internally primed non annotated features : ' + localdate + ' at: ' + localtime + ' \n')
   
	# Combine annotateed and novel APA sites #
	os.system('cat ' + outDir + fkey + '_AnnotatedAPA.bed'.replace('.bed', '.cleaned.bed')+' '+ outDir + fkey + '_novelAPA.bed'.replace('.bed', '.cleaned.bed') + ' >' + outDir + fkey + '_denovoAPAsites.bed')
	pooled = BedTool(outDir + fkey + '_denovoAPAsites.bed')
	pooled = pooled.sort()
	pooled.saveas(outDir + fkey + '_denovoAPAsites.bed')
	os.system('rm ' + outDir + fkey + '_novelAPA.* ' + outDir + fkey + '_AnnotatedAPA.* ')
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished compiling annotated and un-annotated AP sites: ' + localdate + ' at: ' + localtime + ' \n')
	

	###################################
	#Map PA sites #
	###################################
	if MarkFeatures(gtf,fkey,outDir,outDir + fkey + '_denovoAPAsites.bed'):	
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Completed mapping PolyA sites to gene features : '+localdate+' at: ' + localtime+' \n')
		pass
	else:
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write("\nError in mapping PolyA sites to gene features ...\n")
		print ("\nError in mapping PolyA sites to gene features ...\n")
		exit()
	
	####################################
	# Skip selected PA sites #
	####################################
	ignore_features="".join(args_ignore).replace(" ","").split(",")
	if len(ignore_features) >=1:
		PA_FT_df=pd.read_csv(outDir + fkey + '_denovoAPAsites.bed',sep="\t",header=None,index_col=None)
		PA_FT_df.columns=['chr','start','end','GeneID','Feature','strand']
		for ig in ignore_features:
			#print(ig)
			PA_FT_df=PA_FT_df[~PA_FT_df['Feature'].str.contains(ig)]
	PA_FT_df.to_csv(outDir + fkey + '_denovoAPAsites.bed',sep="\t",header=False,index=False)
	
	df = pd.read_csv(outDir + fkey + '_denovoAPAsites.bed', sep='\t', index_col=None, header=None, names=['Chr', 'Start', 'End', 'gene_id', 'Feature', 'Strand'])
	df = df.drop_duplicates(subset=['Chr', 'Start', 'End', 'Strand'])
	df['GeneID'] = df['gene_id'] + '@' + df['Chr'] + '_' + df['Start'].apply(str) + '_' + df['End'].apply(str) + '_' + df['Strand']+ '@' + df['Feature']
	df = df[['GeneID', 'Chr', 'Start', 'End', 'Strand']]
	df.to_csv(outDir + fkey + '_denovoAPAsites.saf', sep='\t', index=False, header=None)
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	os.system("rm "+outDir + fkey +".Features.bed ")
	logfile.write('# Finished making SAF file: ' + localdate + ' at: ' + localtime + ' \n')
	return (1)

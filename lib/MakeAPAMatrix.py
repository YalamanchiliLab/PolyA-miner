# Make APA count matrix #

import pandas as pd
import numpy as np
import concurrent.futures as cf 
import os, sys, glob, time, subprocess

def GeneFilter(pars):
	gene=pars[0]; tdf=pars[1]; controls=pars[2]; treated=pars[3]; mge=float(pars[4])
	if tdf.shape[0] <= 1:
		return ('0')
	ct = list(tdf[controls].mean(axis=1))
	#ctn = len(ct)
	cts = sum(ct)
	tr = list(tdf[treated].mean(axis=1))
	#trn = len(tr)
	trs = sum(tr)
	if cts < mge or trs < mge:
		return ('0')
	else:
		return (gene)

def count2proplist(lst):
	tot = sum(lst)
	if tot > 0:
		for i in range(0, len(lst)):
			lst[i] = lst[i] / tot

		return lst
	temp = []
	for i in range(0, len(lst)):
		temp.append(0)

	return (temp)

def PropFilter(gene, df, controls, treated, mip):
	tdf = df.copy()
	tdf = tdf[tdf['gene_id'] == gene]
	if tdf.shape[0] <= 1:
		return
	else:
		ct = count2proplist(list(tdf[controls].mean(axis=1)))
		tr = count2proplist(list(tdf[treated].mean(axis=1)))
		passlist = []
		for i in range(0, len(ct)):
			if ct[i] <= mip and tr[i] <= mip:
				continue
			passlist.append(i)
		return (tdf.iloc[passlist, :])

def check_PoverA(row, P, A, M):
	count = 0.0
	nrow = list(row)[2:]
	for i in nrow:
		if i >= float(A):
			count = count + 1.0
		if i < float(M):
			return (0)
	if count >= float(P) * len(nrow):
		return (1)
	else:
		return(0)

def PoverA(pars):
	pars[0][pars[4]] = pars[0].apply(lambda row: check_PoverA(row, pars[1], pars[2], pars[3]), axis=1)
	return (pars[0])

def APA_counts(f):
	file = f.split(',')[0]
	saf = f.split(',')[1]
	if len(f.split(',')) == 3:
		outfile=f.split(',')[2]+file.split("/")[-1].replace('.bam','')
	else:
		outfile=file.replace('.sorted.bam','')
	log = subprocess.run(['featureCounts','-a',saf,'-F','SAF','-f','-O','--readExtension5','0','--readExtension3','0','-M','-s','0','-T','5','-o',outfile+'.APA.Counts.txt',file],stderr=subprocess.PIPE,stdout=subprocess.PIPE,shell=False)
	output = log.stderr.decode(encoding='utf-8').split("\n")
	for oline in output:
		if "Total reads :" in oline or "Total alignments :" in oline:
			n=oline.strip().split(": ")[1].split(" ")[0]
			return(n)

def MakeMatrix(outDir, npc, fkey, PA_P, PA_A, M, controls, treated, mip, mge, mode,samples,logfile):
	if os.path.isfile(outDir + fkey + '_denovoAPAsites.saf'):
		saf=outDir + fkey + '_denovoAPAsites.saf'
	else:
		saf=outDir + fkey + '_APA.saf'

	if mode =="fastq":
		files = glob.glob(outDir + '*.bam')
		for i in range(0, len(files)):
			files[i] = files[i] + ',' + saf
		for i in range(0, len(controls)):
			controls[i] = controls[i].replace('.fastq.gz', '')
		for i in range(0, len(treated)):
			treated[i] = treated[i].replace('.fastq.gz', '')
	
	if mode == "bam":
		files=samples
		for i in range(0, len(files)):
			files[i] = files[i] + ',' + saf+ ',' + outDir
		for i in range(0, len(controls)):
			controls[i] = controls[i].split("/")[-1].replace('.bam', '')
		for i in range(0, len(treated)):
			treated[i] = treated[i].split("/")[-1].replace('.bam', '')
	
	with cf.ProcessPoolExecutor(max_workers=npc) as (executor):
		result = list(executor.map(APA_counts, files))
	
	# Lib Size #
	fwls=open(outDir+"LibSize.txt","w")
	for i in range(0, len(files)):
		fwls.write(files[i].split(",")[0].split("/")[-1].split('.'+mode)[0])
		if i < len(files)-1:
			fwls.write("\t")
	fwls.write("\n"+"\t".join(result))
	fwls.close()

	os.system('rm ' + outDir + fkey + "*.saf " + outDir + '*.txt.summary')
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished feature counts: ' + localdate + ' at: ' + localtime + ' \n')

	samples = glob.glob(outDir +'*.APA.Counts.txt')
	df = pd.read_csv(samples[0], comment='#', sep='\t', index_col=None)
	df.columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']+[list(df.columns)[-1].split("/")[-1].split('.'+mode)[0]]
	for s in samples[1:]:
		dftemp = pd.read_csv(s, comment='#', sep='\t', index_col=None)
		df[list(dftemp.columns)[-1].split("/")[-1].split('.'+mode)[0]] = dftemp[list(dftemp.columns)[-1]].copy()

	df['feature_id'] = df['Chr'] + '_' + df['Start'].astype(str) + '_' + df['End'].astype(str) + '_' + df['Strand']
	df[['gene_id', 'isoform']] = df.Geneid.str.split('@', expand=True)	
	design = ['feature_id', 'gene_id'] + controls + treated
	df = df[design]
	df.to_csv(outDir + fkey + '_APA.CountMatrix.txt', sep='\t', index=False)
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished making data matrix : ' + localdate + ' at: ' + localtime + ' \n')
	
	# Clean files #
	log = subprocess.run(["rm"]+samples,stderr=subprocess.DEVNULL,shell=False)
	
	# Low proportion APA filtering #
	genes = list(set(list(df['gene_id'])))
	fdf = pd.DataFrame()
	for gene in genes:
		t = PropFilter(gene, df, controls, treated, mip)
		try:
			if t.shape[0] > 0:
				fdf = fdf.append(t, ignore_index=True)
		except:
			pass
	fdf.to_csv(outDir + fkey + '_APA.CountMatrix.PR.txt', sep='\t', index=False)
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished filterin low prop APA sites : ' + localdate + ' at: ' + localtime + ' \n')

	# PoverA filter.#
	cdf = fdf[['feature_id', 'gene_id'] + controls]
	tdf = fdf[['feature_id', 'gene_id'] + treated]
	passlist = [[cdf, str(PA_P), str(PA_A), str(M), 'C'], [tdf, str(PA_P), str(PA_A), str(M), 'T']]

	with cf.ProcessPoolExecutor(max_workers=npc) as (executor):
		result = list(executor.map(PoverA, passlist))
	df = pd.DataFrame()
	df = pd.merge(result[0], result[1], on=['feature_id', 'gene_id'], how='outer')
	df = df[design + ['C'] + ['T']]
	df = df.fillna(0)
	df = df[(df['C'] > 0) | (df['T'] > 0)]
	df = df.drop(columns=['C', 'T'])
	df.to_csv(outDir + fkey + '_APA.CountMatrix.PA.PR.txt', sep='\t', index=False)
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished PoverA filtering : ' + localdate + ' at: ' + localtime + ' \n')
	
	# gene exp filter #
	genes = list(set(list(df['gene_id'])))
	passlist=[]
	#genes = ['MECP2', 'VMA21', 'LAMC1', 'PAK1', 'PAK2']
	for gene in genes:
		passlist.append([gene, df[df['gene_id'] == gene].copy(), controls, treated, mge])
	with cf.ProcessPoolExecutor(max_workers=npc) as (executor):
		result = list(executor.map(GeneFilter, passlist))
	
	result=list(set(result))
	if '0' in result:
		result.remove('0')
	fdf = pd.DataFrame()
	fdf=df[df['gene_id'].isin(result)] 
	fdf.to_csv(outDir + fkey + '_APA.CountMatrix.GFil.PA.PR.txt', sep='\t', index=False)
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finished filtering out low exp genes : ' + localdate + ' at: ' + localtime + ' \n')
	
	os.system("rm "+outDir + fkey + '_APA.CountMatrix.PR.txt '+outDir + fkey + '_APA.CountMatrix.PA.PR.txt ' +outDir + fkey + '_APA.CountMatrix.txt')
	return (1)

# Testing -BetaBinomial and iNMF #

import os, sys, glob, time, math, warnings
import numpy as np, pandas as pd, statistics
import concurrent.futures as cf
import statsmodels.stats.multitest as multi
from scipy.stats import beta, chi2
from sklearn.decomposition import NMF
from sklearn import preprocessing
from numpy import array, dot, arccos, clip
from numpy.linalg import norm
from scipy.stats import combine_pvalues


def geneLevelBB(gdf):
	Gene=list(gdf['gene_id'])[0]
	PolASites_No=str(len(gdf['feature_id']))
	PolyASites=",".join(list(gdf['feature_id']))
	DeltaU= ",".join(list(map(str,list(gdf['DeltaU']))))
	if gdf['DeltaU'].max() >= abs(gdf['DeltaU'].min()):
		i=gdf['DeltaU'].argmax()
	else:
		i=gdf['DeltaU'].argmin()
	maxPolyA=gdf['feature_id'].iloc[i]
	maxDeltaU=gdf['DeltaU'].iloc[i]
	Stat,Gpval= combine_pvalues(list(gdf["P-value"]), "stouffer", weights = np.absolute(list(gdf["DeltaU"])))
	return([Gene,PolASites_No,PolyASites,DeltaU,maxPolyA,maxDeltaU,Gpval])

def runBBtest(apa_matrix, nc, nt, out, key, npc,logfile):
	try:
		from rpy2.robjects.packages import importr
		import rpy2.robjects as ro
		from rpy2.robjects import pandas2ri
		pandas2ri.activate()
		ro.r['options'](warn=-1)
		countdata=importr('countdata')
		apa_df=pd.read_csv(apa_matrix,sep="\t",header=0,index_col=None)
		temp = apa_df.drop(columns = ["feature_id"]).groupby(["gene_id"]).sum().add_suffix("_G")
		apa_df = apa_df.merge(temp, on="gene_id", how = "inner")
		cols=apa_df.columns
		# Adjust counts #
		apa_mat=apa_df.values
		for row in apa_mat:
			ck=row[2:2+nc]; tk=row[2+nc:2+nc+nt]
			cn=row[2+nc+nt:2+nc+nt+nc]; tn=row[2+nc+nt+nc:2+nc+nt+nc+nt]
			if (ck==cn).all():
				cn=cn+1
				row[2+nc+nt:2+nc+nt+nc]=cn
			if (tk==tn).all():
				tn=tn+1
				row[2+nc+nt+nc:2+nc+nt+nc+nt]=tn
		apa_df2 = pd.DataFrame(apa_mat, columns = cols)
		
		ns=nc+nt
		K=cols[2:(ns+2)];N=cols[(ns+2):((ns+2)*2)]
		Kdf=apa_df2[K]+1;Ndf=apa_df2[N]+1
		Kdf_matrix=Kdf.values; Ndf_matrix=Ndf.values
		Ndf = pd.DataFrame(Ndf_matrix, columns=N)
		grp=[]
		for i in range(0,nc):
			grp.append("C1")
		for i in range(0,nt):
			grp.append("C2")
		Kdf_r=ro.conversion.py2rpy(Kdf)
		Ndf_r=ro.conversion.py2rpy(Ndf)
		n_threads=ro.conversion.py2rpy(npc)
		alternative = ro.conversion.py2rpy("two.sided")
		verbose = ro.vectors.BoolVector([False])
		results=countdata.bb_test(Kdf_r,Ndf_r,ro.StrVector(grp),alternative,n_threads,verbose)

		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Finished BB test on APAsites: ' + localdate + ' at: ' + localtime + ' \n')
		apa_df['P-value']=results[results.names.index('p.value')][:].tolist()
		apa_df["P-value"] = pd.to_numeric(apa_df["P-value"])
		#apa_df = apa_df[apa_df['P-value'].notna()]
		cN=Ndf[N[:nc]]
		cN.columns = cN.columns.str.rstrip('_G')
		cK=Kdf[K[:nc]]
		tN=Ndf[N[nc:]]
		tN.columns = tN.columns.str.rstrip('_G')
		tK=Kdf[K[nc:]]
		normC = (cK/cN).mean(axis=1)
		normT= (tK/tN).mean(axis=1)
		deltaU = normT - normC
		apa_df["DeltaU"] = deltaU
		apa_df.to_csv(out.rstrip("/")+"/"+key+"_APA_BBstats.txt",sep="\t",header=True,index=False)
		# Gene level #
		passlist=[]
		for gene in list(set(apa_df['gene_id'])):
			passlist.append(apa_df[apa_df['gene_id'] == gene])
		with cf.ProcessPoolExecutor(max_workers=npc) as (executor):
			result = list(executor.map(geneLevelBB, passlist))
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Finished Gene level stats: ' + localdate + ' at: ' + localtime + ' \n')
		gene_df=pd.DataFrame(result, columns = ["Gene","No.PAsites","PolyASites","DeltaU","MaxPASite","MaxDeltaU","G-Pval"])
		gene_df['AdjG-Pval'] = (multi.multipletests(gene_df['G-Pval'], method='fdr_bh', is_sorted=False, returnsorted=False))[1].tolist()		
		gene_df.to_csv(out.rstrip("/")+"/"+key+'_Gene_Stats.txt',sep="\t",header=True,index=False)
		return(1)
	except:
		return(0)



def estBetaParams(data):
	mu = sum(data) / float(len(data))
	var = float(np.var(data))
	alpha = ((1 - mu) / var - 1 / mu) * (mu * mu)
	beta = alpha * (1 / mu - 1)
	return (alpha, beta)


def ComputeLRT(adjmatrix, ncr, ntr, nruns):
	#print(adjmatrix)
	adjmatrix[adjmatrix == 0] = [1]
	adjmatrix[adjmatrix == 100] = [99]
	#print(adjmatrix)
	intra_cluster = []
	for i in range(0, ncr):
		for j in range(i + 1, ncr):
			intra_cluster.append(adjmatrix[(i, j)])

	for i in range(ncr, ntr + ncr):
		for j in range(i + 1, ncr + ntr):
			intra_cluster.append(adjmatrix[(i, j)])

	inter_cluster = []
	for i in range(0, ncr):
		for j in range(ncr, ncr + ntr):
			inter_cluster.append(adjmatrix[(i, j)])

	intra_cluster = np.array(intra_cluster) / nruns
	inter_cluster = np.array(inter_cluster) / nruns
	Stability = np.sum(intra_cluster) / (np.sum(intra_cluster) + np.sum(inter_cluster))
	intra_cluster = [0.01 if x == 0 else x for x in intra_cluster]
	intra_cluster = [0.99 if x == 1 else x for x in intra_cluster]
	inter_cluster = [0.01 if x == 0 else x for x in inter_cluster]
	inter_cluster = [0.01 if x == 1 else x for x in inter_cluster]

	#print(np.var(intra_cluster),np.var(inter_cluster))
	if float(np.var(intra_cluster)) <= 0.0000001 :
		#print(intra_cluster)
		for i in range(0,len(intra_cluster)):
			intra_cluster[i] = intra_cluster[i] + (i+1)/1000
		#print(intra_cluster)
	if float(np.var(inter_cluster)) <= 0.0000001 :
		for i in range(0,len(inter_cluster)):
			inter_cluster[i] = inter_cluster[i] + (i+1)/1000

	a1, b1 = estBetaParams(intra_cluster)
	l1 = beta.logpdf(intra_cluster, a1, b1)
	a2, b2 = estBetaParams(inter_cluster)
	l2 = beta.logpdf(inter_cluster, a2, b2)
	a3, b3 = estBetaParams(np.append(intra_cluster, inter_cluster))
	l0 = beta.logpdf(np.append(intra_cluster, inter_cluster), a3, b3)
	LR = 2 * ((np.sum(l1) + np.sum(l2))-np.sum(l0))
	if math.isnan(LR):
		exit()
	p = chi2.sf(LR, 2)
	#print(a1,b1,a2,b2,a3,b3)
	return (LR, p, Stability)


def count2prop(sdf):
	X = sdf.values
	X_sums = np.sum(X, axis=0)
	for i in range(X.shape[0]):
		for j in range(X.shape[1]):
			X[(i, j)] = X[(i, j)] / X_sums[j]
	return (X)

def RunNMF(x):
	sdf = x[1][x[1]['gene_id'] == x[0]]
	x[1] = 'NULL'
	if sdf.shape[0] > 1:
		APAsites = list(sdf['feature_id'])
		strand = APAsites[0].split('_')[-1]
		Samples = list(sdf.columns[2:])
		sdf = sdf[sdf.columns[2:]]
		osdf = sdf.copy()
		sdf = count2prop(sdf)
		sdf[np.isnan(sdf)] = 0.0
		cvector = np.zeros(shape=x[2] + x[3])
		adjmatrix = np.zeros(shape=(x[2] + x[3], x[2] + x[3]))
		H_g = np.zeros(shape=(2, x[2] + x[3]))
		W_g = np.zeros(shape=(len(APAsites), 2))
		for run in range(0, x[4]):
			model = NMF(n_components=2, init='random', max_iter=1000, solver='cd', beta_loss='frobenius', random_state=None, alpha=0.0, shuffle=True)
			W = model.fit_transform(sdf)
			H = model.components_
			cols = H.shape[1]
			for j in range(0, cols):
				if float(H[(0, j)]) >= float(H[(1, j)]):
					cvector[j] = 0
				else:
					cvector[j] = 1

			for p in range(0, x[2] + x[3]):
				for q in range(p + 1, x[2] + x[3]):
					if cvector[p] == cvector[q]:
						adjmatrix[(p, q)] = adjmatrix[(p, q)] + 1
		#print(x[0])
		LR, p, stability = ComputeLRT(adjmatrix, x[2], x[3], x[4])	
		APAscores = list(abs((np.mean(sdf[:, 0:x[2]], axis=1)) - (np.mean(sdf[:, x[2]:x[2] + x[3]], axis=1))))
		return ([x[0], str(len(APAsites)), (',').join(APAsites), (',').join(map(str, APAscores)), APAsites[APAscores.index(max(APAscores))],str(APAscores[APAscores.index(max(APAscores))]),str(stability), str(LR), str(p)])
	   
def iNMFtest(outDir, fkey, nc, nt, ni, npc, matrix,logfile):	
	try:
		warnings.filterwarnings("ignore")
		df = pd.read_csv(matrix, sep='\t', index_col=None)
		genes = list(set(list(df['gene_id'])))
		dflb = pd.read_csv(outDir.rstrip("/")+"/LibSize.txt", sep='\t', index_col=None)
		cols = df.columns[2:]
		for i in range(0, len(cols)):
			df[cols[i]] = df[cols[i]] / float(dflb[cols[i]][0]) * 1000000.0

		passlist = []
		#genes = ['MECP2', 'VMA21', 'LAMC1', 'PAK1', 'PAK2']
		for gene in genes:
			passlist.append([gene, df, nc, nt, ni, outDir])

		with cf.ProcessPoolExecutor(max_workers=npc) as (executor):
			result = list(executor.map(RunNMF, passlist))
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Finished iNMF test on : ' + localdate + ' at: ' + localtime + ' \n')
		gene_df=pd.DataFrame(result, columns = ["Gene","No.PAsites","PolyASites","APAscores","MaxAPAsite","MaxAPAscore","Stability","LRstat","G-Pval"])
		gene_df['G-Pval'] = gene_df['G-Pval'].astype(float)
		gene_df['AdjG-Pval'] = (multi.multipletests(gene_df['G-Pval'], method='fdr_bh', is_sorted=False, returnsorted=False))[1].tolist()		
		gene_df.to_csv(outDir.rstrip("/")+"/"+fkey+'_Gene_Stats.txt',sep="\t",header=True,index=False)
		return(1)
	except:
		return(0)
		




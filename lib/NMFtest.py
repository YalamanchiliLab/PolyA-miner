
from scipy.stats import beta, chi2
from sklearn.decomposition import NMF
from sklearn import preprocessing
from numpy import array, dot, arccos, clip
from numpy.linalg import norm
import os, sys, glob, time, uuid, math, warnings
from fpdf import FPDF
import numpy as np, pandas as pd, statistics, matplotlib, warnings, matplotlib.gridspec as gridspec, concurrent.futures as cf, matplotlib.pyplot as plt, matplotlib.image as mpimg, statsmodels.stats.multitest as multi


def scale(x, out_range=(0, 1)):
    domain = (
     np.min(x), np.max(x))
    y = (x - (domain[1] + domain[0]) / 2) / (domain[1] - domain[0])
    return (y * (out_range[1] - out_range[0]) + (out_range[1] + out_range[0]) / 2)

def ComputeLengtheningScore(W_g, strand):
    vd = np.asarray((np.zeros(shape=(1, W_g.shape[0])))[0])
    vs = np.asarray((np.zeros(shape=(1, W_g.shape[0])))[0])
    if strand == '+':
        vs[0] = 1
        vs[-1] = 0
        vd[0] = 0
        vd[-1] = 1
    if strand == '-':
        vs[0] = 0
        vs[-1] = 1
        vd[0] = 1
        vd[-1] = 0
    try:
        a = count2proplist(list(W_g[:, 0].flat))
        b = count2proplist(list(W_g[:, 1].flat))
    except:
        exit()

    cos_theta_cp = dot(vs, a) / (norm(vs) * norm(a))
    cos_theta_cd = dot(vd, a) / (norm(vd) * norm(a))
    cos_theta_tp = dot(vs, b) / (norm(vs) * norm(b))
    cos_theta_td = dot(vd, b) / (norm(vd) * norm(b))
    tp = norm(b) * cos_theta_tp
    cp = norm(a) * cos_theta_cp
    td = norm(b) * cos_theta_td
    cd = norm(a) * cos_theta_cd
    return (cp, tp, cd, td, (td - tp) / (cd - cp))


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


def count2proplist(lst):
    tot = sum(lst)
    if tot > 0:
        for i in range(0, len(lst)):
            lst[i] = lst[i] / tot

        return (lst)
    temp = []
    for i in range(0, len(lst)):
        temp.append(0)

    return (temp)


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
        APAmean = np.hstack((sdf[:, 0:x[2]].sum(axis=1, keepdims=1), sdf[:, x[2]:x[2] + x[3]].sum(axis=1, keepdims=1)))
        CP, TP, CD, TD, FC = ComputeLengtheningScore(APAmean, strand)
        APAscores = list(abs((np.mean(sdf[:, 0:x[2]], axis=1)) - (np.mean(sdf[:, x[2]:x[2] + x[3]], axis=1))))
        if x[6]:
            return ([x[0], str(len(APAsites)), ('; ').join(APAsites), ('; ').join(map(str, APAscores)), str(APAscores.index(max(APAscores)) + 1), str(CP), str(TP), str(CD), str(TD), str(FC), str(stability), str(p), adjmatrix])
        else: 
            return ([x[0], str(len(APAsites)), ('; ').join(APAsites), ('; ').join(map(str, APAscores)), str(APAscores.index(max(APAscores)) + 1), str(CP), str(TP), str(CD), str(TD), str(FC), str(stability), str(p)])
       
def NMFtest(outDir, fkey, nc, nt, ni, np, matrix, LS, plots, lf):
    warnings.filterwarnings("ignore")
    logfile=open(lf,"a")
    df = pd.read_csv(matrix, sep='\t', index_col=None)
    genes = list(set(list(df['gene_id'])))
    dflb = pd.read_csv(LS, sep='\t', index_col=None)
    cols = df.columns[2:]
    for i in range(0, len(cols)):
        df[cols[i]] = df[cols[i]] / float(dflb[cols[i]][0]) * 1000000.0

    passlist = []
    #genes = ['MECP2', 'VMA21', 'LAMC1', 'PAK1', 'PAK2']
    for gene in genes:
        passlist.append([gene, df, nc, nt, ni, outDir, plots])

    with cf.ProcessPoolExecutor(max_workers=np) as (executor):
        result = list(executor.map(RunNMF, passlist))
    result=list(result) 
    adjmatrix_dict={} 
    if plots:
        for i in range(0,len(result)):
            adjmatrix_dict[result[i][0]]=result[i].pop()
    fw = open(outDir + '/' + fkey + '.PolyA-miner.Results.txt', 'w')
    fw.write('Gene\tNo. APAsites\tAPAsites\tAPAsitesScores\tMaxAPASwitch\tCP\tTP\tCD\tTD\tFC\tStability\tPvalue\n')
    for p in result:
        try:
            fw.write(('\t').join(p) + '\n')
        except:
            pass
    fw.close()
    df = pd.read_csv(outDir + '/' + fkey + '.PolyA-miner.Results.txt', sep='\t', index_col=None)
    df['Adj-Pvalue'] = (multi.multipletests(df['Pvalue'], method='fdr_bh', is_sorted=False, returnsorted=False))[1].tolist()
    df.to_csv(outDir + '/' + fkey + '.PolyA-miner.Results.txt', sep='\t', index=False)
    localdate = time.strftime('%a %m/%d/%Y')
    localtime = time.strftime('%H:%M:%S')
    logfile.write('# Finished test on : ' + localdate + ' at: ' + localtime + ' \n')
    logfile.close()
    if plots:
        return (1,adjmatrix_dict)
    else:
        return (1)



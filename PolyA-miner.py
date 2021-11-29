# PolyA-miner.py -help # 
# Hari Krishna Y et al., last update 10/28/2021 #

import os, sys
import time,argparse,subprocess, pandas as pd
sys.path.append("/".join(os.path.abspath(sys.argv[0]).split("/")[0:-1])+"/lib")
import CheckDependency, DataProcessing, ExtractPolyAsites
import MakeAPAMatrix, GenePolyAIndex, STest


def check_files (checkfiles,logfile):
	for checkfile in checkfiles:
		if os.path.exists(checkfile):
			pass
		else:
			logfile.write("\nError cannot locate "+checkfile+"\n")
			exit()
	return(1)

def main():
	parser = argparse.ArgumentParser(description='''PolyA-miner v1.1: Inferring alternative poly-adenylation changes
	from 3'Seq data  - Yalamanchili H.K. et al. \n''',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	optional = parser._action_groups.pop()
	required = parser.add_argument_group('Required arguments')
	parser._action_groups.append(optional)

	required.add_argument('-mode',help='Run mode options: \'bam\' to start from mapped data, \'fastq\' to start from raw data',choices=['bam','fastq'],default='bam',required='True',type=str)
	optional.add_argument('-d',help='Base directory of input fastq files. Valid for -mode fastq ',type=str)
	optional.add_argument('-o',help='Out put directory',type=str,default='PolyAminer_OUT')
	required.add_argument('-c1',help='Comma-separated list of condition1 files. Full path for BAMs (index files are also expected) or Just file names for fastq', nargs='+',required='True',type=str)
	required.add_argument('-c2',help='Comma-separated list of condition2 files. Full path for BAMs (index files are also expected) or Just file names for fastq ', nargs='+',required='True',type=str)
	parser.add_argument('-s',help='Strand information 0: un-stranded 1: fwd-strand 2:rev-strand. ',choices=[0,1,2],type=int,default=0)
	
	# Ref. files
	optional.add_argument('-index',help='Reference genome bowtie2 index. Valid for -mode fastq',type=str)
	required.add_argument('-fasta',help='Reference fasta sequence',required='True',type=str)
	required.add_argument('-bed',help='Reference genes bed file',required='True',type=str)
	required.add_argument('-pa',help='PolyA annotations file standard 6 column bed format',type=str)
	
	# Optional #
	optional.add_argument('-umi',help='Length of UMIs, 0 if not used', type=int,default=0)
	optional.add_argument('-apaBlock',help='Window size for annotated polyA sites',type=int, default=30)
	optional.add_argument('-mdapa',help='Cluster distance for annotated polyA sites: Merge polyA sites with in this distance. ',type=int, default=0)
	optional.add_argument('-md',help='Cluster distance for de-novo polyA sites: Merge polyA sites with in this distance',type=int, default=0)
	optional.add_argument('-anchor',help='Overlap in "bp" for mapping de-novo polyA sites to annotated polyA sites ',type=int, default=1)
	# Tuning #
	optional.add_argument('-expNovel',help='Explore novel APA sites 0: only annotated sites 1: de-novo',choices=[1,0],type=int,default=0)
	optional.add_argument('-p',help='No. of processors to use',type=int,default=4)
	optional.add_argument('-ip',help='Internal priming window',type=int, default=50)
	optional.add_argument('-a',help='Internal priming polyA fraction',type=float, default=0.6)
	optional.add_argument('-pa_p',help='pOverA filter: P ',type=float, default=0.6)
	optional.add_argument('-pa_a',help='pOverA filter: A ',type=int, default=5)
	optional.add_argument('-pa_m',help='pOverA filter: M ',type=int, default=2)
	optional.add_argument('-gene_min',help='Min counts per Gene',type=int, default=10)
	optional.add_argument('-apa_min',help='Min. proportion per APA',type=float, default=0.05)
	#
	optional.add_argument('-t',help='Statistical Test- BB: for beta-binomail or iNMF: for iterative NMF ',type=str,default="DM")
	optional.add_argument('-i',help='No. of NMF iterations. Valid only for -t iNMF',type=int,default=100)
	optional.add_argument('-outPrefix',help='Output file/s prefix', default="PolyAminer_Out",type=str)
	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)
	args=parser.parse_args()
	args_dict = vars(args)
	
	if args.mode =='bam' and args.umi >1:
		print("\nUMI option is only valid in fastq mode ...\n")
		exit()
	if args.mode =='bam' and len(args.index) >= 1:
		print("\nindex option is only valid in fastq mode. Skipping ...\n")

	# Check Outdir #
	args.o=args.o.rstrip("/")
	if (os.path.isdir(args.o)):
		log=subprocess.run(['rm','-R',args.o],stderr=subprocess.DEVNULL,shell=False)
	log=subprocess.run(['mkdir',args.o],stderr=subprocess.DEVNULL,shell=False)
	args.o=args.o+"/"
	
	# Logging #
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile = open(args.o+args.outPrefix+'.PolyA-miner.log.txt','w')
	lf=args.o+args.outPrefix+'.PolyA-miner.log.txt'
	logfile.write('# Starting PolyA-miner : '+localdate+' at: ' + localtime+' \n')
	logfile.write('# *********** Arguments ************** \n')
	for k in args_dict.keys():
		if k in ["c1","c2"]:
			logfile.write("\t"+k+' : ' + "\t".join(args_dict[k])+'\n')
		else:
			logfile.write("\t"+k+' : ' + str(args_dict[k])+'\n')
	logfile.write('# *********** ********* ************** \n')

	# CheckDependency # 
	if CheckDependency.checkDep(logfile):
		pass
	else:
		logfile.write('Fix dependencies and run ...\n')
		print('Fix dependencies and run ...\n')
		exit()
		
	# if no ref. polyA use dummy #
	localdate = time.strftime('%a_%m_%d_%Y')
	localtime = time.strftime('%H_%M_%S')
	if args.pa is None:
		if args.expNovel == 1:
			dummy_refPA=args.o+"DummyRefPA_"+localdate+"_"+localtime+".bed"
			fw=open(dummy_refPA,"w")
			fw.write("chr1	634823	634823	D000001	DummyGene	+\n")
			fw.close()
			args.pa=dummy_refPA
		else:
			logfile.write("Err: Reference PolyA file is needed if expNovel is set to 0 ...\n")
			print("Err: Reference PolyA file is needed if expNovel is set to 0 ...\n")
			exit()
	###################################################
	# Module 1A: Arguments check and Data preprocessing 
	###################################################
	# Number of samples #
	controls="".join(args.c1).replace(" ","").split(",")
	nc=len(controls)
	treated="".join(args.c2).replace(" ","").split(",")
	nt=len(treated)

	if args.mode=='bam':
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Run mode: from preprocessed bam files: '+localdate+' at: ' + localtime+' \n')
		# check files #
		checkfiles=controls+treated+[args.fasta]+[args.bed]+[args.pa]
		if check_files(checkfiles,logfile):
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Arguments checked : '+localdate+' at: ' + localtime+' \n')
		
	if args.mode=='fastq':
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Run mode: from raw fastq data : '+localdate+' at: ' + localtime+' \n')
		# check basedir #:
		if os.path.isdir(args.d):
			pass
		else:
			logfile.write("\nError cannot locate "+args.d+"\n")
			exit()
		
		# Check Ref #:
		rg="/".join(args.index.split('/')[:-1])
		if os.path.isdir(rg):
			pass
		else:
			logfile.write("\nError cannot locate "+rg+"\n")
			exit()
		
		# Check files #
		cck=[]
		for c in controls:
			cck.append(args.d+"/"+c)
		tck=[]
		for t in treated:
			tck.append(args.d+"/"+t)
		checkfiles=cck+tck+[args.fasta]+[args.bed]+[args.pa]
		if check_files(checkfiles,logfile):
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Arguments checked : '+localdate+' at: ' + localtime+' \n')
		try:
			samples = controls + treated
			nc = len(controls)
			nt = len(treated)
			baseDir =args.d; outDir=args.o; np=args.p; ref_genome=args.index; fkey=args.outPrefix;
			for s in samples:
				DataProcessing.process_rawfastq(args.d,args.o, s, args.umi,args.p)
				localdate = time.strftime('%a %m/%d/%Y')
				localtime = time.strftime('%H:%M:%S')
				logfile.write('# Finished fastp/umi: ' + localdate + ' at: ' + localtime + 'for ' + s + ' \n')
				DataProcessing.mapping_bowtie2(args.o, s, args.p, args.index,args.umi)
				localdate = time.strftime('%a %m/%d/%Y')
				localtime = time.strftime('%H:%M:%S')
				logfile.write('# Finished mapping on ' + localdate + ' at: ' + localtime + 'for ' + s + ' \n')
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Completed data processing : '+localdate+' at: ' + localtime+' \n')
			pass

		except:
			logfile = open(args.o+args.outPrefix+'.PolyA-miner.log.txt','a')
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Error in data processing module ... \n')
			logfile.close()
			print ("\nError in data processing module ...\n")
			exit()

	################################
	# Module2: Extract PolyA sites #
	################################
	if ExtractPolyAsites.ExtAPAfromPolyA_DB(args.o, args.outPrefix, args.pa, args.apaBlock, args.mdapa,logfile) == 1:
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Completed extracting annotated polyadenylation sites : '+localdate+' at: ' + localtime+' \n')
		pass
	else:
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write("\nError in extracting annotated polyadenylation sites ...\n")
		print ("\nError in extracting annotated polyadenylation sites ...\n")
		exit()
	if args.expNovel == 1:
		if ExtractPolyAsites.ExtNovelAPA(args.o, args.outPrefix, args.bed, args.fasta, args.md, args.anchor, args.ip, args.a, args.p, args.mode,controls+treated,logfile):
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Completed extracting de-novo polyadenylation sites : '+localdate+' at: ' + localtime+' \n')
			pass
		else:
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write("\nError in extracting de-novo polyadenylation sites ...\n")
			print ("\nError in extracting de-novo polyadenylation sites ...\n")
			exit()
	else:
		ExtractPolyAsites.makeSAF(args.o, args.outPrefix)
	
	###################################
	# Module 3: Make APA count matrix #
	###################################
	if MakeAPAMatrix.MakeMatrix(args.o, args.p, args.outPrefix, args.pa_p, args.pa_a, args.pa_m, controls, treated, args.apa_min, args.gene_min, args.mode,controls+treated,logfile) == 1:
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Completed adstracting APA proportions : '+localdate+' at: ' + localtime+' \n')
		pass
	else:
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write("\nError in adstracting APA proportions ...\n")
		print ("\nError in adstracting APA proportions ...\n")
		exit()

	###################################
	# Module 4: Gene level PolyA Index 
	###################################
	if GenePolyAIndex.VectorPro(args.o.rstrip("/")+"/", args.outPrefix, nc, nt, args.p, args.o.rstrip("/")+"/"+args.outPrefix+ '_APA.CountMatrix.GFil.PA.PR.txt', args.o.rstrip("/")+"/LibSize.txt"):
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Completed computing PolyA Index : '+localdate+' at: ' + localtime+' \n')
		pass
	else:
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write("\nError in computing PolyA Index ...\n")
		print ("\nError in computing PolyA Index ...\n")
		exit()

	###################################
	# Module 5: Stat BetaBinomial iNMF 
	###################################
	if args.t == "BB":
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Using BB model: ' + localdate + ' at: ' + localtime + ' \n')
		if STest.runBBtest(args.o.rstrip("/")+"/"+args.outPrefix+ '_APA.CountMatrix.GFil.PA.PR.txt', nc, nt, args.o, args.outPrefix, args.p,logfile):
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Completed beta-binomail testing: ' + localdate + ' at: ' + localtime + ' \n')
		else:
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write("\nError in beta-binomail testing ...\n")
			print ("\nError in beta-binomail testing ...\n")
			exit()

	if args.t =='iNMF':
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Using iNMF model: ' + localdate + ' at: ' + localtime + ' \n')
		if STest.iNMFtest(args.o,args.outPrefix, nc, nt, args.i,args.p, args.o.rstrip("/")+"/"+args.outPrefix+ '_APA.CountMatrix.GFil.PA.PR.txt', logfile):
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Completed iNMF testing: ' + localdate + ' at: ' + localtime + ' \n')
		else:
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write("\nError in iNMF testing ...\n")
			print ("\nError in iNMF testing ...\n")
			exit()

	# Merge 4 and 5 add gene symbols #
	pdata=pd.read_csv(args.o.rstrip("/")+"/"+args.outPrefix+'_PolyA-miner.Results.txt',sep="\t",header=0,index_col=None)
	stat_data=pd.read_csv(args.o.rstrip("/")+"/"+args.outPrefix+"_Gene_Stats.txt",sep="\t",header=0,index_col=None)
	pdata=pd.merge(pdata,stat_data,left_on=["Gene"],right_on=["Gene"],how="outer")
	
	genes=pd.read_csv(args.bed,sep="\t",header=None,index_col=None)
	genes.columns=["Chr","Start","End","Gene","Symbol","Strand"]
	if "." in pdata['Gene'].iloc[0]:
		pdata=pd.merge(pdata,genes,on=["Gene"],how="left")
		pdata=pdata.drop(columns=["Chr","Start","End","Strand"])
	else:
		genes[['Gene','Version']] = genes['Gene'].str.split('.',expand=True)
		pdata=pd.merge(pdata,genes,on=["Gene"],how="left")
		pdata=pdata.drop(columns=["Chr","Start","End","Strand","Version"])
	pdata.to_csv(args.o.rstrip("/")+"/"+args.outPrefix+'_PolyA-miner.Results.txt',sep="\t",header=True,index=False)
	

	# Summary #
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# PolyA-miner results summary: '+localdate+' at: ' + localtime+' \n')
	nsig_c=pdata[pdata['AdjG-Pval']<=0.05].shape[0]
	nsig_s=pdata[(pdata['AdjG-Pval']<=0.05) & (pdata['PolyAIndex'] <0)].shape[0]
	nsig_l=pdata[(pdata['AdjG-Pval']<=0.05) & (pdata['PolyAIndex'] >0)].shape[0]
	logfile.write('# Significant PolyA changes:\t'+str(nsig_c)+"\n"+'# Significant shortening changes:\t'+str(nsig_s)+"\n"+'# Significant lengthening changes:\t'+str(nsig_l)+"\n")
	logfile.write('# Finished PolyA-miner : '+localdate+' at: ' + localtime+' \n')
	logfile.close()

	# Tidy up  #
	files=[args.o.rstrip("/")+"/LibSize.txt",args.o.rstrip("/")+"/"+args.outPrefix+"_Gene_Stats.txt",args.o.rstrip("/")+"/"+args.outPrefix+'.APSitesDB.bed',args.o.rstrip("/")+"/"+args.outPrefix+"_APA.CountMatrix.GFil.PA.PR.txt",args.o.rstrip("/")+"/"+args.outPrefix+"_denovoAPAsites.bed"]
	log=subprocess.run(['rm','-R']+files,stderr=subprocess.DEVNULL,shell=False)
	try:
		subprocess.run(['rm',dummy_refPA]+files,stderr=subprocess.DEVNULL,shell=False)
	except:
		pass
		
if __name__ == "__main__":
	main()

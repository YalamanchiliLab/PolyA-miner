# PolyA-miner.py -help #
# Hari Krishna Y et al., last update 05/10/2020 #

import os, sys, glob, time, uuid, argparse, pip

def main():
	# parse arguments #
	parser = argparse.ArgumentParser(description='''PolyA-miner: accurate assessment of differential alternative poly-adenylation 
	from 3'Seq data using vector projections and non-negative matrix factorization - Yalamanchili H.K. et al. \n''')
	parser.add_argument('-d',help='Base directory of input fastq files ',required='True',type=str)
	parser.add_argument('-o',help='Out put directory',required='True',type=str,default='PolyA-miner_OUT')
	parser.add_argument('-pa',help='PolyA annotations file',required='True',type=str)
	parser.add_argument('-index',help='Reference genome bowtie2 index',required='True',type=str)
	parser.add_argument('-fasta',help='Reference fasta sequence',required='True',type=str)
	parser.add_argument('-bed',help='Reference genes bed file',required='True',type=str)
	parser.add_argument('-c1',help='Comma-separated list of condition1 fastq files', nargs='+',required='True',type=str)
	parser.add_argument('-c2',help='Comma-separated list of condition2 fastq files', nargs='+',required='True',type=str)
	
	# optional defaults #
	parser.add_argument('-i',help='No. of NMF iterations ',type=int,default=100)
	parser.add_argument('-p',help='No. of processors to use',type=int,default=4)
	parser.add_argument('-expNovel',help='Explore novel APA sites',action='store_true', default=False)
	parser.add_argument('-ip',help='Internal priming window',type=int, default=50)
	parser.add_argument('-pa_p',help='pOverA filter: P ',type=float, default=0.6)
	parser.add_argument('-pa_a',help='pOverA filter: A ',type=int, default=5)
	parser.add_argument('-pa_m',help='pOverA filter: M ',type=int, default=2)
	parser.add_argument('-apa_min',help='Min. proportion per APA',type=float, default=0.05)
	parser.add_argument('-gene_min',help='Min counts per Gene',type=int, default=10)
	parser.add_argument('-mdapa',help='Cluster distance for annotated polyA sites',type=int, default=0)
	parser.add_argument('-md',help='Cluster distance for de-novo polyA sites',type=int, default=25)
	parser.add_argument('-apaBlock',help='Block size for annotated polyA sites',type=int, default=50)
	parser.add_argument('-anchor',help='Anchor length for de-novo polyA sites',type=int, default=1)
	parser.add_argument('-key',help='Output file key', default="PolyAminer.Out",type=str)
	args=parser.parse_args()
	
	# check outdir #
	if (os.path.isdir(args.o)):
		pass
	else:
		os.system("mkdir "+args.o)
	args.o=args.o+"/"
	
	# Log  #
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile = open(args.o+args.key+'.PolyA-miner.log.txt','w')
	lf=args.o+args.key+'.PolyA-miner.log.txt'
	logfile.write('# Starting PolyA-miner : '+localdate+' at: ' + localtime+' \n')
	
	# Check dependency #
	dependency=['pandas','Cython','pybedtools','scipy','sklearn','statsmodels']
	installed_packages = pip.get_installed_distributions()
	ip=''
	for p in installed_packages:
		ip=ip+str(p)+" "
	for package in dependency:
		if package in ip:
			pass
		else:
			logfile.write("Package "+package+" not found. Installing ....\n")
			try:		
				os.system("pip3 install "+package)
			except:
				logfile.write("Failed installing "+package+" ....\n")
	
	# import lib #
	sys.path.append("/".join(sys.argv[0].split("/")[0:-1])+"/lib")
	import DataProcessing, ExtAnnotatedAPA, ExtNovelAPA, MakeMatrix, NMFtest
	
	# Check arguments #
	controls="".join(args.c1).replace(" ","").split(",")
	nc=len(controls)
	treated="".join(args.c2).replace(" ","").split(",")
	nt=len(treated)
	
	# check basedir #:
	if os.path.isdir(args.d):
		pass
	else:
		logfile.write("\nError cannot locate "+args.d+"\n")
		exit()
		
	# check ref #:
	rg="/".join(args.index.split('/')[:-1])
	if os.path.isdir(rg):
		pass
	else:
		logfile.write("\nError cannot locate "+rg+"\n")
		exit()

	# check files #
	cck=[]
	for c in controls:
		cck.append(args.d+"/"+c)
	tck=[]
	for t in treated:
		tck.append(args.d+"/"+t)
	
	checkfiles=cck+tck+[args.fasta]+[args.bed]+[args.pa]
	for checkfile in checkfiles:
		if os.path.exists(checkfile):
			pass
		else:
			logfile.write("\nError cannot locate "+checkfile+"\n")
			exit()
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Arguments checked : '+localdate+' at: ' + localtime+' \n')
	logfile.close()

	# Module 1: Data preprocessing only if bams or not found #
	samples=controls+treated
	flag=1
	for s in samples:
		s=s.replace(".fastq.gz","")
		if (os.path.exists(args.o+s+"/"+s+".sorted.bam")):
			pass
		else:
			flag=0

	if flag ==1:
		logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Already found bam files moving on : '+localdate+' at: ' + localtime+' \n')
		logfile.close()
		pass
	else:
		if DataProcessing.DataProcessing(args.d, args.o, args.p, args.index, controls, treated, args.key, lf) == 1:
			logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Completed data processing : '+localdate+' at: ' + localtime+' \n')
			logfile.close()
			pass
		else:
			logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Error in data processing module ... \n')
			logfile.close()
			print ("\nError in data processing module ...\n")
			exit()


	# Module2: De novo or DB based PolyA sites #
	if args.expNovel:
		if ExtNovelAPA.ExtNovelAPA(args.o, args.key, args.pa, args.bed, args.fasta, args.apaBlock, args.mdapa, args.md, args.anchor, args.ip, args.p, lf) == 1:
			logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Completed de novo identification of polyadenylation sites : '+localdate+' at: ' + localtime+' \n')
			logfile.close()
			pass
		else:
			logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write("\nError in de novo identification of polyadenylation sites ...\n")
			logfile.close()
			print ("\nError in de novo identification of polyadenylation sites ...\n")
			exit()
	else:
		if ExtAnnotatedAPA.ExtAnnotatedAPA(args.o, args.key, args.pa, args.apaBlock, args.mdapa, lf) == 1:
			logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write('# Completed extracting annotated polyadenylation sites : '+localdate+' at: ' + localtime+' \n')
			logfile.close()
			pass
		else:
			logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
			localdate = time.strftime('%a %m/%d/%Y')
			localtime = time.strftime('%H:%M:%S')
			logfile.write("\nError in extracting annotated polyadenylation sites ...\n")
			logfile.close()
			print ("\nError in extracting annotated polyadenylation sites ...\n")
			exit()


	# Module 3: Make matrix #
	if MakeMatrix.MakeMatrix(args.o, args.p, args.key, args.pa_p, args.pa_a, args.pa_m, controls, treated, args.apa_min, args.gene_min, lf) == 1:
		logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Completed adstracting APA proportions : '+localdate+' at: ' + localtime+' \n')
		logfile.close()
		pass
	else:
		logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write("\nError in adstracting APA proportions ...\n")
		logfile.close()
		print ("\nError in adstracting APA proportions ...\n")
		exit()
	

	# Module 4: itterative NMF - test  #
	matrix=args.o + args.key + '_APA.CountMatrix.GFil.PA.PR.txt'
	status=NMFtest.NMFtest(args.o, args.key, nc, nt, args.i, args.p, matrix, args.o+"LibSize.txt", False, lf)
	if status == 1:
		logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write('# Completed NMF module : '+localdate+' at: ' + localtime+' \n')
		logfile.close()
		pass
	else:
		logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
		localdate = time.strftime('%a %m/%d/%Y')
		localtime = time.strftime('%H:%M:%S')
		logfile.write("\nError in NMF module ...\n")
		logfile.close()
		print ("\nError in NMF module ...\n")
		exit()


	# tidy up  #
	os.system("rm "+args.o+"*APSitesDB.* "+args.o+"LibSize.txt")
	os.system("mv "+glob.glob(args.o+"*.CountMatrix.GFil.PA.PR.txt")[0]+" "+glob.glob(args.o+"*.CountMatrix.GFil.PA.PR.txt")[0].replace(".GFil.PA.PR",""))
	if args.expNovel:
		os.system("rm "+args.o+"*_denovoAPAsites.bed")
	logfile = open(args.o+args.key+'.PolyA-miner.log.txt','a')
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	logfile.write('# Finishing PolyA-miner : '+localdate+' at: ' + localtime+' \n')
	logfile.close()

if __name__ == "__main__":
	main()

# Check PolyA-miner dependency packages #
import os, sys, glob
import pkg_resources,subprocess


def checkDep(logfile):
	dependency=['pandas','numpy','Cython','pybedtools','scipy','sklearn','statsmodels','uuid','rpy2']
	installed_packages = " ".join([str(d) for d in pkg_resources.working_set])
	for package in dependency:
		if package in installed_packages:
			pass
		else:
			logfile.write("Package "+package+" not found. Installing ....\n")
			try:		
				os.system("pip3 install "+package)
			except:
				logfile.write("Err: Failed installing "+package+" ....\n")
				print("Failed installing "+package+" ....\n")
				return(0)

	if subprocess.call(['which','featureCounts'],stderr=subprocess.PIPE,stdout=subprocess.PIPE,shell=False) ==1:
		logfile.write("Err: Missing featureCounts.")
		print("Err: Missing featureCounts ...")
		return(0)
	
	if subprocess.call(['which','bedtools'],stderr=subprocess.PIPE,stdout=subprocess.PIPE,shell=False) == 1:
		logfile.write("Err: Missing bedtools.")
		print("Err: Missing bedtools. Install bedtools.")
		return(0)
		
	import rpy2
	from rpy2.robjects.packages import importr
	rpy2.robjects.r['options'](warn=-1)

	try:
		countdata=importr('countdata')
	except:
		logfile.write("Err: Missing R package 'countdata'.")
		print("Err: Missing R package 'countdata'.")
		return(0)

	return(1)


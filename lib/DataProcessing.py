# Processing raw fastq files #
import os, sys, glob
import subprocess


def process_rawfastq(args_d, args_o, s, umi, mc):
	fastq4fastp = args_d + '/' + s
	umi_pat=""
	if umi > 0:
		for i in range(0, umi):
			umi_pat = umi_pat + 'N'
		fastq4fastp = args_o + s.replace('.fastq.gz', '') + '.UMIex.fastq.gz'
		log = subprocess.run(['umi_tools','extract','--extract-method=string','--bc-pattern='+umi_pat,'-L',args_o + s.replace('.fastq.gz', '') + '_UMIextract.log.txt','-I',args_d + '/' + s,'-S',fastq4fastp],stderr=subprocess.DEVNULL,shell=False)
	if int(mc) > 16:
		mc = '16'
	log = subprocess.run(['fastp','-i',fastq4fastp,'-a',"AGATCGGAAGAGC",'-f','6','-g','-l','40','-j',args_o + s.split('/')[(-1)] + '_trim-report.json','-h',args_o + s.split('/')[(-1)] + '_trim-report.html','-w',mc,'-Q','-o',args_o + s.split('/')[(-1)] + '_trim.fastq'],stderr=subprocess.DEVNULL,shell=False)
	log = subprocess.run(['pigz','-p',mc,args_o + s.split('/')[(-1)] + '_trim.fastq'],stderr=subprocess.DEVNULL,shell=False)
	return (1)

def mapping_bowtie2(args_o, s, mc, ref_genome, umi):
	log = subprocess.run(["rm"]+glob.glob(args_o+"*trim-report.json"),stderr=subprocess.DEVNULL,shell=False)
	#d = args_o + s.split('/')[(-1)].replace('.fastq.gz', '')
	cmd = 'bowtie2 -p ' + str(mc) + ' --very-sensitive-local -x ' + ref_genome + ' -U ' + args_o + s.split('/')[(-1)] + '_trim.fastq.gz -S ' + args_o + s.split('/')[(-1)] + '.sam'
	log = subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
	cmd = 'samtools view -@ ' + str(mc) + ' -bS ' + args_o + s.split('/')[(-1)] + '.sam' + ' -o ' + args_o + s.split('/')[(-1)] + '.bam'
	log = subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
	cmd = 'samtools sort -@ ' + str(mc) + ' ' + args_o + s.split('/')[(-1)] + '.bam' + ' -o ' + args_o + s.split('/')[(-1)] + '.sorted.bam'
	log = subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
	cmd = 'samtools index -@ ' + str(mc) + ' ' + args_o + s.split('/')[(-1)] + '.sorted.bam'
	log = subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
	if umi > 0:
		cmd = 'umi_tools dedup -I ' + args_o + s.split('/')[(-1)] + '.sorted.bam' + ' -L ' + args_o + s.split('/')[(-1)] + '.DeDuplog.txt' + ' -S ' + args_o + s.split('/')[(-1)] + '.dedup.sorted.bam' + ' --umi-separator="_"'
		os.system(cmd)
		#log = subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
		subprocess.run(["rm",args_o + s.split('/')[(-1)] + '.sorted.bam.bai'],stderr=subprocess.DEVNULL,shell=False)
		cmd = 'mv ' + args_o + s.split('/')[(-1)] + '.dedup.sorted.bam ' +args_o + s.split('/')[(-1)] + '.sorted.bam'
		subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
		cmd = 'samtools index -@ ' + str(mc) + ' ' + args_o + s.split('/')[(-1)] + '.sorted.bam'
		log = subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
		subprocess.run(["rm",args_o + s.replace('.fastq.gz', '') + '.UMIex.fastq.gz'],stderr=subprocess.DEVNULL,shell=False)
	cmd='rm ' + args_o + s.split('/')[(-1)] + '_trim.fastq.gz ' + args_o + s.split('/')[(-1)] + '.sam ' + args_o + s.split('/')[(-1)] + '.bam'
	log = subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
	cmd = 'mv ' + args_o + s.split('/')[(-1)] + '.sorted.bam ' +args_o + s.split('/')[(-1)] + '.bam'
	log = subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
	cmd = 'mv ' + args_o + s.split('/')[(-1)] + '.sorted.bam.bai ' +args_o + s.split('/')[(-1)] + '.bam.bai'
	log = subprocess.run(cmd.split(" "),stderr=subprocess.DEVNULL,shell=False)
	return (1)

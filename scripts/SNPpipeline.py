#!/usr/bin/env python
import sys, argparse
import os
import json

#import sys,os
#import pandas as pd
#from Bio import SeqIO
#from pysam import VariantFile
#from Bio.SeqRecord import SeqRecord
#from Bio.Seq import Seq

#from pyfaidx import Fasta
#from pyfaidx import FastaVariant
#import vcf

parser=argparse.ArgumentParser(prog='SNPpipeline.py',  
							   usage="%(prog)s [options] -i jsoninputfilename -o outputfolder",
							   description='''Script that compares the given genomes (fastq files) with a reference genome based on SNP number differences.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the json file containing all information about the genomes.')
parser.add_argument('-o','--out', help='The folder for the results.') 
parser.add_argument('-prefix','--prefix', default="allsnps", help='The folder for the results.') 

args=parser.parse_args()
inputfilename= args.input

prefix= args.prefix
resultfoldername= args.out


def SNP_calling_longshot(genome,reference):
	refasm = reference["asm"]
	refasm = reference["asm"]
	if not refasm.startswith("/"):
		refasm=dataprefix + refasm
	reffai = reference["fai"]
	if not reffai.startswith("/"):
		reffai=dataprefix + reffai	
	#index the reference with samtools faidx 
	command = samtools + " faidx " + refasm
	if not os.path.exists(reffai):
		print(command)
		os.system(command)
	fastqfile=genome['fastq']
	if not fastqfile.startswith("/"):
		fastqfile=dataprefix + fastqfile
	result = resultfoldername + "/" + genome['id']
	samfile=result + ".sam"
	bamfile=result + ".bam"
	baifile=result + ".bam.bai"
	vcffile=result + ".vcf"
	snpfile=result + ".snp.vcf"
	filteredsnpfile=result + ".filtered.snp.vcf"
	
	#align each genome to the reference
	command = minimap2 + " -a " +  refasm + " " + fastqfile + " > "  + samfile
	#command = "bwa mem -M -R \"@RG\\tID:group1\\tSM:sample" + samplename + "\\tPL:illumina\\tLB:lib1\\tPU:unit1\" " + referenceidx + " " + result + "_R1_P.fastq.gz " + result + "_R2_P.fastq.gz > " + result + ".sam"
	if not os.path.exists(samfile):
		print("Map the reads to the reference: " + command)
		os.system(command)	
	#sort the sam file	
	command = samtools + " view -bS " + samfile + " | " + samtools + " sort -O BAM -o " + bamfile
	if os.path.exists(samfile) and not os.path.exists(bamfile):
		print("Sort the sam file: " + command)
		os.system(command)
	#index bam file
	command=samtools + " index " + bamfile
	if not os.path.exists(baifile):
		print("Index bam file: " + command)
		os.system(command)
	#search for variants	
	command=longshot + " --bam " + bamfile + " --ref " + refasm + " --out " + vcffile
	if not os.path.exists(vcffile):
		print("Look for variants: " + command)
		os.system(command)	
		os.system("sed -i \'s/SAMPLE/" + genome['id'] + "/g\' " + vcffile)
	#remove indels
	command = 	vcftools + " --vcf " + vcffile + " --remove-indels --recode --recode-INFO-all"
	if not os.path.exists(snpfile):
		print("Remove indels: " + command)
		os.system(command)
		os.system("mv out.recode.vcf " + snpfile)
		
	#Filtering out the SNP with QUAL <20
	command =  vcftools + " --vcf " + snpfile + " --minQ 20 --recode --recode-INFO-all"
	if not (os.path.exists(filteredsnpfile) or os.path.exists(filteredsnpfile + ".gz") ):
		print("Filter the snp: " + command)
		os.system(command)
		os.system("mv out.recode.vcf " + filteredsnpfile)
		
def SNP_calling_gatk(genome,reference):	
	refidx = reference["idx"]
	if not refidx.startswith("/"):
		refidx=dataprefix + refidx
	refasm = reference["asm"]
	if not refasm.startswith("/"):
		refasm=dataprefix + refasm
	reffai =reference["fai"]
	if not reffai.startswith("/"):
		reffai=dataprefix + reffai	
	refdict = reference["dict"]
	if not refdict.startswith("/"):
		refdict=dataprefix + refdict
	#index the reference
	command=bwa + " index -p " + refidx + " " + refasm
	if not os.path.exists(refidx + ".bwt"):
		print("Index the reference: " + command)
		os.system(command)
	#index the reference with samtools faidx 
	command = samtools + " faidx " + refasm
	if not os.path.exists(reffai):
		print("Index the reference: " + command)
		os.system(command)
	#create a dictionary for reference
	command= gatk + " CreateSequenceDictionary -R " + refasm + " -O " + refdict
	if not os.path.exists(refdict):
		print("Make a dictionary out of the reference: " + command)
		os.system(command)
		
	fastq1=genome["fastq1"]
	fastq2=genome["fastq2"]
	if not fastq1.startswith("/"):
		fastq1=dataprefix + fastq1
	if not fastq2.startswith("/"):
		fastq2=dataprefix + fastq2	
	genomeid = genome["id"]
	result = resultfoldername + "/" + genomeid
	samfile= result + ".sam"
	bamfile=result + ".bam"
	markedfile=result + ".marked.bam"
	markedtxtfile= result + ".marked.txt"
	vcffile=result + ".vcf"
#	reportfile=result + "_report.table"
#	bqsrbamfile=result + "_bqsr.bam"
#	bqsrvcffile=result + "_bqsr.vcf"
	snpfile=result + ".snp.vcf"
	filteredsnpfile=result + ".filtered.snp.vcf"
#	indelfile=result + ".indel.vcf"
#	filteredindelfile=result + ".filtered.indel.vcf"
	#mapping
	command = bwa + " mem -M -R \"@RG\\tID:group1\\tSM:" + genomeid + "\\tPL:illumina\\tLB:lib1\\tPU:unit1\" " + refidx + " " + fastq1 + " " + fastq2 + " > " + samfile
	if not os.path.exists(samfile):
		print("Map to the reference: " + command)
		os.system(command)
	#2, sort the SAM file
	command = samtools + " view -bS " + samfile + " | samtools sort -O BAM -o " + bamfile
	if os.path.exists(samfile) and not os.path.exists(bamfile):
		print("Sort the samfile: " + command)
		os.system(command)
	#3, Mark duplicates
	command = picard + " MarkDuplicates I=" + bamfile + " O=" + markedfile + " M=" +  markedtxtfile + " VALIDATION_STRINGENCY=LENIENT"
	if os.path.exists(bamfile) and not os.path.exists(markedfile):
		print("Mark duplicates: " + command)		
		os.system(command)
	#4, index the marked bam file
	command = samtools + " index " + markedfile
	if not os.path.exists(markedfile + ".bai"):
		print("Index duplicate marked bam file: " + command)
		os.system(command)
	#5, Predict the vcf file
	command = gatk + " HaplotypeCaller -R " + refasm + " -I " + markedfile + " -O " + vcffile
	if os.path.exists(markedfile) and not os.path.exists(vcffile):
		print("Predict vcf :" + command)	
		os.system(command)
#	#6, BaseRecalibrator
#	command = gatk + " BaseRecalibrator -R " + refasm + " -I " + markedfile + "  -O " + reportfile + " -known-sites " + vcffile + " -use-original-qualities" 
#	if os.path.exists(vcffile) and not os.path.exists(reportfile):
#		print("BaseRecalibrator:" + command)
#		os.system(command)
#	#7, ApplyBQSR
#	command = gatk + " ApplyBQSR -R " + refasm + " -I " + markedfile + " -bqsr-recal-file " + reportfile + " -O " + bqsrbamfile
#	if os.path.exists(reportfile) and not os.path.exists(bqsrbamfile):
#		print("ApplyBQSR:" + command)
#		os.system(command)
#	#8, Predict vcf file based on BQSR
#	command = gatk + " HaplotypeCaller -R " + refasm + " -I " + bqsrbamfile + " -O " + bqsrvcffile
#	if os.path.exists(bqsrbamfile) and not os.path.exists(bqsrvcffile):
#		print("Predict vcf based on BQSR:" + command)
#		os.system(command)
		
	#9, Select the snp
	command = gatk + " SelectVariants -select-type SNP -V " + vcffile + " -O " + snpfile
	if os.path.exists(vcffile) and not os.path.exists(snpfile):
		print("Select the snp: " + command)
		os.system(command)
		
	#10, Filter the snp
	command = gatk + " VariantFiltration -V " + snpfile + " --filter-expression \"QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 ||ReadPosRankSum < -8.0\" --filter-name \"Filter\" -O " + filteredsnpfile
	if os.path.exists(snpfile) and not (os.path.exists(filteredsnpfile) or os.path.exists(filteredsnpfile +".gz")):
		print("Filter the snp: " + command)
		os.system(command)
#	#11, Select the indel
#	command = gatk + " SelectVariants -select-type INDEL -V " + bqsrvcffile + " -O "+ indelfile
#	if os.path.exists(bqsrvcffile) and not  os.path.exists(indelfile): 
#		print("Select the indel :"  + command)
#		os.system(command)
#	#12, Filter the indel
#	command = gatk + " VariantFiltration -V " + indelfile + " --filter-expression \"Q< 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5|| ReadPosRankSum < -8.0\" --filter-name \"Filter\" -O " + filteredindelfile
#	if os.path.exists(indelfile) and not os.path.exists(filteredindelfile):
#		print("Filter the indel: " + command)
#		os.system(command)

def Download(genome):
	fastq=""
	fastq1=""
	fastq2=""
	if "fastq" in genome.keys():
		fastq= genome["fastq"]
		if not fastq.startswith("/"):
			fastq=dataprefix + fastq
	if "fastq1" in genome.keys():	
		fastq1=genome["fastq1"]
		if not fastq1.startswith("/"):
			fastq1=dataprefix + fastq1
	if "fastq2" in genome.keys():	
		fastq2= genome["fastq2"]
		if not fastq2.startswith("/"):
			fastq2=dataprefix + fastq2	
	if os.path.exists(fastq):
		return
	elif os.path.exists(fastq1) and os.path.exists(fastq2):
		return
	genomefolder=dataprefix + genome["id"]
	if os.path.exists(genomefolder):
		command="mkdir " + genomefolder
		print(command)
		os.system(command)
	sra=""
	if "sra" in genome.keys():
		sra=genome["sra"]
	if sra!= "":
		print("Download sra for " + genome["id"])
		#download sra
		srafilename= sra + ".sra"
		command=prefetch + " " + sra
		print(command)
		os.system(command)
		if os.path.exists("~/ncbi/public/sra/" + srafilename):
			#move srafile to the working folder
			os.system("mv ~/ncbi/public/sra/" + srafilename + " " + genomefolder + "/")
		else:
			command = "curl \"https://sra-pub-run-odp.s3.amazonaws.com/sra/" + sra + "/" + sra +"\" > " + genomefolder + "/" + srafilename
			print(command)
			os.system(command)				
		command = fasterq_dump + " --split-files " + genomefolder + "/" + srafilename + " -O " + genomefolder
		print(command)
		os.system(command)
		#remove sra file
		os.system("rm " + genomefolder + "/" + srafilename)	
	if not (os.path.exists(fastq) or (os.path.exists(fastq1) and os.path.exists(fastq2))):	
		print("Cannot find sra file for " + genome['id'])
		logfile=open("SNPpipeline.log","w")
		logfile.write("Cannot download sra for " + genome["id"] + "\n")
		logfile.close()
	
def TabixVcf(genome):
	genomeid = genome["id"]
	result = resultfoldername + "/" + genomeid
	snpfile=result + ".filtered.snp.vcf"
	gsnpfile=snpfile + ".gz" 
	tabixsnpfile=snpfile + ".gz.tbi" 
	cmd="bgzip " + snpfile
	if os.path.exists(snpfile) and (not os.path.exists(gsnpfile)):
		print(cmd)
		os.system(cmd)
	cmd=tabix + " -p vcf " + gsnpfile
	if os.path.exists(gsnpfile) and (not os.path.exists(tabixsnpfile)):
		print(cmd)
		os.system(cmd)
	return gsnpfile

####################MAIN####################
with open(inputfilename) as data_file:
	data = json.load(data_file)	
dataprefix=data["dataprefix"]	
bwa=data['bwa']
#trim_file=data['trim_file']
minimap2=data['minimap2']	
longshot=data['longshot']	
samtools=data['samtools']	
trim_file=data['trim_file']
gatk=data['gatk']
vcftools=data['vcftools']
vcfmerge=data['vcfmerge']
tabix=data['tabix']

picard=data['picard']	
prefetch=""
if "prefetch" in data.keys():
	prefetch=data["prefetch"]
fasterq_dump =""
if "fasterq-dump" in data.keys():
	fasterq_dump = data["fasterq-dump"]
genomes=data['genomes']	
reference=data['reference']
dataprefix=data["dataprefix"]

if not os.path.exists(resultfoldername):
	os.system("mkdir " + resultfoldername)
#clades=clade.split(",")
vcffilenames=" "
#numbers={}
samples=[]
longreads=[]
shortreads=[]
for key in genomes.keys():
	genome=genomes[key]
	#download fastq file if not existed
	Download(genome)
	fastq=""
	fastq1=""
	fastq2=""
	if "fastq" in genome.keys():
		fastq= genome["fastq"]
		if not fastq.startswith("/"):
			fastq=dataprefix + fastq
	if "fastq1" in genome.keys():	
		fastq1=genome["fastq1"]
		if not fastq1.startswith("/"):
			fastq1=dataprefix + fastq1
	if "fastq2" in genome.keys():	
		fastq2= genome["fastq2"]
		if not fastq2.startswith("/"):
			fastq2=dataprefix + fastq2	
	print("Looking for the SNPs of " + key)	
	#SNP calling
	if os.path.exists(fastq):
		SNP_calling_longshot(genome,reference)
		longreads.append(key)
	elif os.path.exists(fastq1) and os.path.exists(fastq2):
		SNP_calling_gatk(genome,reference)
		shortreads.append(key)
	samples.append(key)	
	#Tabix vcf file
	gvcffile=TabixVcf(genome)	
	vcffilenames=vcffilenames + gvcffile + " "
	
print("The number of long read sequencing samples: " + str(len(shortreads)))	
print(shortreads)
	
print("The number of long read sequencing samples: " + str(len(longreads)))
print(longreads)

#Merge the snps
allsnps=resultfoldername + "/" + prefix + ".vcf.gz"	
if not os.path.exists(allsnps):
	cmd = vcfmerge + vcffilenames + " | bgzip -c > " + allsnps
	print(cmd)
	os.system(cmd)
	os.system("gunzip " + allsnps)
	print("The vcf are saved in file " + allsnps +".")
else:
	print("File " +  allsnps + " exists. Please delete the file " + allsnps + " if you want to create a new vcf file.")

	
	

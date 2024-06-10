#!/usr/bin/env python
# FILE: convertTab2Json.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 Nay 2024
import sys, argparse
import os
import json


parser=argparse.ArgumentParser(prog='convertTab2Json.py',  
							   usage="%(prog)s [options] -i tabdelimitedfilename -o jsonfilenane",
							   description='''Script that converts tab delimited file name containing information about genomes into json format file.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the tab delimited file name containing information about genomes.')
parser.add_argument('-o','--out', help='the json file containing information about the genomes which can be used as the input of the SNP pipeline.') 

args=parser.parse_args()

inputfilename= args.input
outputfilename= args.out


def Tab2Json(inputfilename):
	inputdict={}
	
	inputdict.setdefault("dataprefix","")
	inputdict.setdefault("longshot","longshot")
	inputdict.setdefault("bwa","bwa")
	inputdict.setdefault("minimap2","minimap2")
	inputdict.setdefault("samtools","samtools")
	inputdict.setdefault("vcftools","vcftools")
	inputdict.setdefault("vcfmerge","vcf-merge")
	inputdict.setdefault("tabix","tabix")
	inputdict.setdefault("gatk","gatk")
	inputdict.setdefault("prefetch","prefetch")
	inputdict.setdefault("fasterq-dump","fasterq-dump")
	inputdict.setdefault("picard","picard.jar")
	reference={}
	reference.setdefault("id","")
	reference.setdefault("asm","")
	reference.setdefault("fai","")
	reference.setdefault("idx","")
	reference.setdefault("dict","")
	inputdict.setdefault("reference",reference)
	genomes={}
	inputfile=open(inputfilename)
	headers=next(inputfile).split("\t")
	for line in inputfile:
		genome={}
		texts=line.split("\t")
		genomeid=""
		for i in range(len(texts)):
			if (i >= len(headers) or i >= len(texts)):
				continue
			genome[headers[i].rstrip()]=texts[i].rstrip()
			if genomeid=="":
				genomeid=texts[i].rstrip()
		fastq=""
		fastq1=""
		fastq2=""
		if "fastq" in genome.keys():
			fastq=genome["fastq"]
		if "fastq1" in genome.keys():
			fastq1=genome["fastq1"]
		if "fastq2" in genome.keys():
			fastq2=genome["fastq2"]	
		if "sra" in genome.keys():
			sra=genome["sra"]
			if sra!="":
				if fastq1=="":
					fastq1=genomeid + "/" + sra + "_1.fastq" 
				if fastq2=="":
					fastq2=genomeid + "/" + sra + "_2.fastq"
				if fastq=="":
					fastq=genomeid + "/" + sra + ".fastq"		
		genome["fastq"]=fastq	
		genome["fastq1"]=fastq1	
		genome["fastq2"]=fastq2
		genomes.setdefault(genomeid,genome)
	inputdict.setdefault("genomes",genomes)	
	return inputdict
	
####################MAIN####################
inputdict=Tab2Json(inputfilename)
#save the input dictionary
with open(outputfilename,"w") as json_file:
	json.dump(inputdict,json_file,indent=2)

  

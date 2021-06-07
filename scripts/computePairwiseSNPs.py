#!/usr/bin/env python
# FILE: computePairwiseSNPs.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 May 2021
import sys, argparse
import os
from Bio import SeqIO
import json

parser=argparse.ArgumentParser(prog='computePairwiseSNPs.py',  
							   usage="%(prog)s [options] -i jsoninputfilename -f fastafile -o outputfile",
							   description='''Script that compares the SNP number differences.''',
							   epilog="""Written by Duong Vu d.vu@wi.knaw.nl""",
   )

parser.add_argument('-i','--fasta', help='The fasta alignment file.') 
parser.add_argument('-o','--out', help='The folder for the results.') 
args=parser.parse_args()
outputfilename= args.out
fastafilename=args.fasta

record_dict = SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))

outputfile = open(outputfilename, "w")
firstline = ""
for name1 in record_dict.keys():
	firstline = firstline + "\t" + name1
outputfile.write(firstline + "\n")	
for name1 in record_dict.keys():
	seq1 = str(record_dict[name1].seq)
	print(name1 + "\t" + str(len(seq1)) + "\t" + str(len(seq1.replace("N",""))))
	line = name1
	for name2 in record_dict.keys():
		seq2=str(record_dict[name2].seq)
		if len(seq1) != len(seq2):
			print("Not a good alignment file!")
			break
		count = 0
		for i in range(0,len(seq1)-1):
			if seq1[i] != seq2[i]:
				count = count + 1
		line = line + "\t" +  str(count)
	outputfile.write(line + "\n")

outputfile.close()

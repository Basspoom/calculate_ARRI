#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 20230208
@author: hanmz
@version: V1.0
@Anhui Medical University
@Calculate the ARRI for evalating the quality of water based on the paired-end metagenomic reads.
'''
from __future__ import division
import glob
import os
import datetime
import types
import csv
import sys
import argparse
import time
csv.field_size_limit(500 * 1024 * 1024)
def main():
	print ("Usage: python2 calculate_ARRI.py -i1 Sample_R1.fastq.gz -i2 Sample_R2.fastq.gz -sa samplename -t number -s similarity -c coverage -o1 ARRI.value")
	parser = argparse.ArgumentParser()
	parser.add_argument("-i1","--inputfile1")
	parser.add_argument("-i2","--inputfile2")
	parser.add_argument("-sa","--samplename")
	parser.add_argument("-t","--thread",type=int,default=48)
	parser.add_argument("-s","--similarity",type=float,default=80)
	parser.add_argument("-c","--coverage",type=float,default=70)
	parser.add_argument("-o1","--output1",default="ARRI.value")
	ARRI = parser.parse_args()
	#At the begain of analysis, several bioinformatic tools should be installed, such as fastp, MEGAHIT, Prodigal, BLAST, DeePaC(version=0.10.1), seqkit, metabat2, bowtie2.
	#Please download the databases from ourlab or contract to us to obtain the essential databases, including SARG, IS, aclame, and integron databases.
	#In this pepline, we assumed that you have properly installed the relevant tools and databases.
	#If you have any questions, please do not hesitate to contract with us, thanks.
	#1. Conduct the QC of sequencing reads with fastp.
	time1=time.strftime('%Y%m%d',time.localtime(time.time()))
	outfile1=open(r"1fastp%s.sh"%(time1),'w')
	outfile1.writelines("#!/bin/bash\n")
	outfile1.writelines("# -*- coding: utf-8 -*-\n")
	outfile1.writelines("##Author Your_Name\n")
	outfile1.writelines("##Time is %s for QC!\n"%time1)
	outfile1.writelines("fastp -i %s -o %s_R1.fastq.gz -I %s -O %s_R2.fastq.gz -w %s -h %s.fastp.html -j %s.fastp.json\n"%(ARRI.inputfile1, ARRI.samplename, ARRI.inputfile2, ARRI.samplename, ARRI.thread, ARRI.samplename, ARRI.samplename))
	outfile1.close()
	os.system("sh 1fastp%s.sh"%time1)
	#2. Assembly the paired-end reads to obtain the contig fasta
	os.system("mkdir assembly results Abundance")
	time1=time.strftime('%Y%m%d',time.localtime(time.time()))
	outfile1=open(r"2megahit%s.sh"%(time1),'w')
	outfile1.writelines("#!/bin/bash\n")
	outfile1.writelines("# -*- coding: utf-8 -*-\n")
	outfile1.writelines("##Author Your_Name\n")
	outfile1.writelines("##Time is %s for assembly!\n"%time1)
	outfile1.writelines("megahit -1 %s_R1.fastq.gz -2 %s_R2.fastq.gz -o ./assembly/%s/ -m 0.9 --min-contig-len 500 -t %s --continue --presets meta-large --out-prefix %s\n"%(ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.thread, ARRI.samplename))
	outfile1.close()
	os.system("sh 2megahit%s.sh"%time1)
	#3. Predict the genes and proteins
	time1=time.strftime('%Y%m%d',time.localtime(time.time()))
	outfile1=open(r"3prodigal%s.sh"%(time1),'w')
	outfile1.writelines("#!/bin/bash\n")
	outfile1.writelines("# -*- coding: utf-8 -*-\n")
	outfile1.writelines("##Author Your_Name\n")
	outfile1.writelines("##Time is %s for predicting genes and proteins!\n"%time1)
	outfile1.writelines("prodigal -i ./assembly/%s/%s.contigs.fa -o ./assembly/%s/%s.genes -d ./assembly/%s/%s.nucl.faa -a ./assembly/%s/%s.proteins.faa -p meta\n"%(ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.close()
	os.system("sh 3prodigal%s.sh"%time1)
	#4. Predicting the potential ARGs againest SARG.
	time1=time.strftime('%Y%m%d',time.localtime(time.time()))
	outfile1=open(r"4predict_ARGs%s.sh"%(time1),'w')
	outfile1.writelines("#!/bin/bash\n")
	outfile1.writelines("# -*- coding: utf-8 -*-\n")
	outfile1.writelines("##Author Your_Name\n")
	outfile1.writelines("##Time is %s for predicting ARGs!\n"%time1)
	outfile1.writelines("blastp -query ./assembly/%s/%s.proteins.faa "%(ARRI.samplename, ARRI.samplename))
	outfile1.writelines("-db ~/database/sarg ")
	outfile1.writelines("-out ./results/%s_ARG_blastp.txt "%ARRI.samplename)
	outfile1.writelines("-evalue 1e-5 -max_target_seqs 1 -num_threads %s -outfmt '6 qseqid qseq qlen qstart qend sseqid sgi ppos qcovs pident length sseq evalue bitscore stitle'\n"%ARRI.thread)
	outfile1.close()
	os.system("sh 4predict_ARGs%s.sh"%time1)
	#5 Predicting the MGE components.
	time1=time.strftime('%Y%m%d',time.localtime(time.time()))
	outfile1=open(r"5predict_MGEs%s.sh"%(time1),'w')
	outfile1.writelines("#!/bin/bash\n")
	outfile1.writelines("# -*- coding: utf-8 -*-\n")
	outfile1.writelines("##Author Your_Name\n")
	outfile1.writelines("##Time is %s for predicting MGEs!\n"%time1)
	##5.1 MGE_IS.out
	outfile1.writelines("blastn -db ~/database/ISfinder_database/IS.nucl.database ")
	outfile1.writelines("-query ./assembly/%s/%s.nucl.faa -out ./results/%s_IS.blast.txt "%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.writelines("-evalue 1e-10 -max_target_seqs 1 -num_threads %s "%ARRI.thread)
	outfile1.writelines("-outfmt '6 qseqid qseq qlen qstart qend sseqid sgi ppos  qcovs pident length sseq evalue bitscore stitle'\n")
	##5.2 MGE_aclame.out
	outfile1.writelines("blastn -db ~/database/aclame_database/aclame.nucl.database ")
	outfile1.writelines("-query ./assembly/%s/%s.nucl.faa -out ./results/%s_aclame.blast.txt "%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.writelines("-evalue 1e-10 -max_target_seqs 1 -num_threads %s "%ARRI.thread)
	outfile1.writelines("-outfmt '6 qseqid qseq qlen qstart qend sseqid sgi ppos  qcovs pident length sseq evalue bitscore stitle'\n")
	##5.3 MGE_integron.out
	outfile1.writelines("blastn -db ~/database/Integron_database/integron.nucl.database ")
	outfile1.writelines("-query ./assembly/%s/%s.nucl.faa -out ./results/%s_Integron.blast.txt "%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.writelines("-evalue 1e-10 -max_target_seqs 1 -num_threads %s "%ARRI.thread)
	outfile1.writelines("-outfmt '6 qseqid qseq qlen qstart qend sseqid sgi ppos qcovs pident length sseq evalue bitscore stitle'\n")
	outfile1.close()
	os.system("sh 5predict_MGEs%s.sh"%time1)
	#6 Predicting HBPs
	time1=time.strftime('%Y%m%d',time.localtime(time.time()))
	outfile1=open(r"6predict_HBPs%s.sh"%(time1),'w')
	outfile1.writelines("#!/bin/bash\n")
	outfile1.writelines("# -*- coding: utf-8 -*-\n")
	outfile1.writelines("##Author Your_Name\n")
	outfile1.writelines("##Time is %s for predicting HBPs!\n"%time1)
	outfile1.writelines("deepac predict -s ./assembly/%s/%s.contigs.fa \n"%(ARRI.samplename, ARRI.samplename))
	outfile1.writelines("deepac filter ./assembly/%s/%s.contigs.fa ./assembly/%s/%s.contigs_predictions.npy -t 0.90 -p -o ./results/%s_deepac_90.fasta\n"%(ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.writelines("grep '>' ./results/%s_deepac_90.fasta | awk -F \" \" '{printf(\"%s\\n\",$1)}' | awk -F '>' '{printf(\"%s\\n\",$2)}' > ./results/%s_deepac_90_name.txt\\n"%(ARRI.samplename))
	outfile1.close()
	os.system("sh 6predict_HBPs%s.sh"%time1)
	#7 deal with the BLAST results
	def deal_blast_results(infile_name,outfile_name):
		file2= list(csv.reader(file('./results/%s'%infile_name),delimiter="\t"))
		nrow2 = len(open('./results/%s'%infile_name).readlines())
		outfile3=open(r"./results/%s_contigs_name.txt"%outfile_name,'w')
		outfile4=open(r"./results/%s_genes_name.txt"%outfile_name,'w')
		contigs_name=[]
		for i in range (0,nrow2):
			if (float(file2[i][7])>=80 and float(file2[i][8])>=70):
				gene_name=file2[i][0]
				outfile4.writelines("%s\n"%gene_name)
				str2=file2[i][0].split("_")
				str3=str2[0]+"_"+str2[1]
				if str3 not in contigs_name:
					contigs_name.append(str3)
				else:
					continue
		outfile4.close()
		for item1 in contigs_name:
			outfile3.writelines("%s\n"%item1)
		outfile3.close()
	deal_blast_results("%s_ARG_blastp.txt"%ARRI.samplename,"%s_ARG"%ARRI.samplename)
	deal_blast_results("%s_IS.blast.txt"%ARRI.samplename,"%s_IS"%ARRI.samplename)
	deal_blast_results("%s_aclame.blast.txt"%ARRI.samplename,"%s_aclame"%ARRI.samplename)
	deal_blast_results("%s_Integron.blast.txt"%ARRI.samplename,"%s_Integron"%ARRI.samplename)
	os.system("awk '{print $1}' ./results/%s_ARG_blastp.txt > ./results/%s_ARG_genes_name.txt"%(ARRI.samplename, ARRI.samplename))
	os.system("awk -v FS='/_' '{print $1,$2}' %s_ARG_genes_name.txt > %s_ARG_contigs_name_notuniq.txt"%(ARRI.samplename, ARRI.samplename))
	os.system("cat %s_ARG_contigs_name_notuniq.txt | awk -v OFS='_' '{print $1,$2}' | uniq > %s_ARG_contigs_name.txt"%(ARRI.samplename, ARRI.samplename))
	###多个文件取并集
	os.system("cat ./results/%s_IS_contigs_name.txt ./results/%s_aclame_contigs_name.txt ./results/%s_Integron_contigs_name.txt | sort | uniq -u > ./results/%s_MGE_contigs_name.txt"%(ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.samplename))
	os.system("cat ./results/%s_IS_genes_name.txt ./results/%s_aclame_genes_name.txt ./results/%s_Integron_genes_name.txt | sort | uniq -u > ./results/%s_MGE_genes_name.txt"%(ARRI.samplename, ARRI.samplename, ARRI.samplename, ARRI.samplename))
	def calculate_contigs_abundance(number, samplename):
		outfile1=open(r"calculate_contigs_abundance_%s.sh"%samplename,'w')
		outfile1.writelines("#!/bin/bash\n")
		outfile1.writelines("# -*- coding: utf-8 -*-\n")
		outfile1.writelines("##Author Your_Name\n")
		outfile1.writelines("##Time is for calculate_contigs_abundance!\n")
		outfile1.writelines("bowtie2-build --threads %s ./assembly/%s/%s.contigs.fa ./assembly/%s/%s_contigs_fa\n"%(number, samplename, samplename, samplename, samplename))
		outfile1.writelines("bowtie2 -p %s --local --very-sensitive-local -x ./assembly/%s/%s_contigs_fa -1 %s_R1.fastq.gz -2 %s_R2.fastq.gz -S ./assembly/%s/%s_contigs.sam\n"%(number, samplename, samplename, samplename, samplename, samplename, samplename))
		outfile1.writelines("samtools sort --threads %s ./assembly/%s/%s_contigs.sam -o ./assembly/%s/%s_contigs.sort.bam\n"%(number, samplename, samplename, samplename, samplename))
		outfile1.writelines("jgi_summarize_bam_contig_depths --outputDepth ./results/%s_contigs_coverage.txt ./assembly/%s/%s_contigs.sort.bam\n"%(samplename, samplename, samplename))
		outfile1.close()
		os.system("sh calculate_contigs_abundance_%s.sh"%samplename)
	def calculate_genes_abundance(number, samplename):
		outfile1=open(r"calculate_genes_abundance_%s.sh"%samplename,'w')
		outfile1.writelines("# -*- coding: utf-8 -*-\n")
		outfile1.writelines("##Author Your_Name\n")
		outfile1.writelines("##Time is for calculate_genes_abundance!\n")
		outfile1.writelines("bowtie2-build --threads %s ./assembly/%s/%s.nucl.faa ./assembly/%s/%s_nucl_fa\n"%(number, samplename, samplename, samplename, samplename))
		outfile1.writelines("bowtie2 -p %s --local --very-sensitive-local -x ./assembly/%s/%s_nucl_fa -1 %s_R1.fastq.gz -2 %s_R2.fastq.gz -S ./assembly/%s/%s_nucl.sam\n"%(number, samplename, samplename, samplename, samplename, samplename, samplename))
		outfile1.writelines("samtools sort --threads %s ./assembly/%s/%s_nucl.sam -o ./assembly/%s/%s_nucl.sort.bam\n"%(number, samplename, samplename, samplename, samplename))
		outfile1.writelines("jgi_summarize_bam_contig_depths --outputDepth ./results/%s_nucl_coverage.txt ./assembly/%s/%s_nucl.sort.bam\n"%(samplename, samplename, samplename))
		outfile1.close()
		os.system("sh calculate_genes_abundance_%s.sh"%samplename)
	calculate_contigs_abundance("%s"%ARRI.thread, "%s"%ARRI.samplename)
	calculate_genes_abundance("%s"%ARRI.thread, "%s"%ARRI.samplename)
	#两个文件取交集
	#awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1],$2}' 1.txt 2.txt >3.txt
	#######ARG-MGE
	os.system("cat ./results/%s_ARG_contigs_name.txt ./results/%s_MGE_contigs_name.txt | sort | uniq -d > ./results/%s_ARG_MGE_contigs_name.txt"%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	#######ARG-HBP
	os.system("cat ./results/%s_ARG_contigs_name.txt ./results/%s_deepac_90_name.txt | sort | uniq -d > ./results/%s_ARG_HBP_contigs_name.txt"%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	#######ARG-MGE-HBP
	os.system("cat ./results/%s_ARG_HBP_contigs_name.txt ./results/%s_MGE_contigs_name.txt | sort | uniq -d > ./results/%s_ARG_MGE_HBP_contigs_name.txt"%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	######obtain the ARG, MGE, HBP, genes and contigs abundance
	outfile1=open(r"obtain_various_genes_contigs_coverage_%s.sh"%ARRI.samplename,'w')
	outfile1.writelines("# -*- coding: utf-8 -*-\n")
	outfile1.writelines("##Author Your_Name\n")
	outfile1.writelines("##Time is for obtain ARG, MGE,genes_abundance!\n")
	outfile1.writelines("grep -w -F -f ./results/%s_ARG_genes_name.txt ./results/%s_nucl_coverage.txt > ./results/%s_ARG_gene_coverage.txt\n"%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.writelines("grep -w -F -f ./results/%s_MGE_genes_name.txt ./results/%s_nucl_coverage.txt > ./results/%s_MGE_gene_coverage.txt\n"%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.writelines("grep -w -F -f ./results/%s_deepac_90_name.txt ./results/%s_contigs_coverage.txt > ./results/%s_deepac_contigs_coverage.txt\n"%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	###
	outfile1.writelines("grep -w -F -f ./results/%s_ARG_MGE_contigs_name.txt ./results/%s_contigs_coverage.txt > ./results/%s_ARG_MGE_contigs_coverage.txt\n"%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.writelines("grep -w -F -f ./results/%s_ARG_HBP_contigs_name.txt ./results/%s_contigs_coverage.txt > ./results/%s_ARG_HBP_contigs_coverage.txt\n"%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.writelines("grep -w -F -f ./results/%s_ARG_MGE_HBP_contigs_name.txt ./results/%s_contigs_coverage.txt > ./results/%s_ARG_MGE_HBP_contigs_coverage.txt\n"%(ARRI.samplename, ARRI.samplename, ARRI.samplename))
	outfile1.close()
	os.system("sh obtain_various_genes_contigs_coverage_%s.sh"%ARRI.samplename)
	def obtain_abundance(samplename, genetype1, genetype2):
		a=os.popen("awk -F '\t' '{sum += $2};END {print sum}' ./results/%s_%s_%s_coverage.txt"%(samplename,genetype1, genetype2))
		#awk -F '\t' '{sum += $3};END {print sum}' S01_coverage.txt   求一列的和
		x=a.read()
		b=x.split(':')
		c=b[0].split('.')[0]
		e=c.split(' ')[0]
		a.close()
		return e
	Abundance_ARG=float(obtain_abundance("%s"%ARRI.samplename, "ARG", "gene"))
	Abundance_MGE=float(obtain_abundance("%s"%ARRI.samplename, "MGE", "gene"))
	Abundance_HBP=float(obtain_abundance("%s"%ARRI.samplename, "deepac", "contigs"))
	Abundance_ARG_MGE=float(obtain_abundance("%s"%ARRI.samplename, "ARG_MGE", "contigs"))
	Abundance_ARG_HBP=float(obtain_abundance("%s"%ARRI.samplename, "ARG_HBP", "contigs"))
	Abundance_ARG_MGE_HBP=float(obtain_abundance("%s"%ARRI.samplename, "ARG_MGE_HBP", "contigs"))
	ARRI_value=(Abundance_ARG_MGE+Abundance_ARG_HBP+Abundance_ARG_MGE_HBP)/(Abundance_ARG+Abundance_MGE+Abundance_HBP)*10000
	print("ARRI value for %s is %s\n"%(ARRI.samplename,ARRI_value))
	outfile1=open(r"%s.txt"%ARRI.output1,'w')
	outfile1.writelines("%s\t%s\n"%(ARRI.samplename, ARRI_value))
	outfile1.close()

if __name__ == "__main__":
	main()

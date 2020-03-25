#!/bin/python
#coding:utf-8
import os

index1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACG'
index2='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGT'
index3='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGC'
index4='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCA'
index5='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTG'
index6='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAAT'
index7='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATC'
index8='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGA'
index9='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAG'
index10='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTT'
index11='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTAC'
index12='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTA'
index13='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACA'
index14='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGT'
index15='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGA'
index16='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCC'
index18='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCAC'
index19='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACG'
index20='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTT'
index21='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGA'
index22='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTA'
index23='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGAT'
index25='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATAT'
index27='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTT'
primer5='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT'
SR_primer='GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT'

#os.system('ls *.R1.fastq.gz | while read id; do(echo "$id" >> dict.csv); done')

sample_dict = {}
with open("dict.ultra.csv","r") as f:
	for i in f:	
		j = i.strip().split("\t")
		key = j[0]
		sample_dict[key] = eval("index" + j[1])
print(sample_dict)

for i in sample_dict.keys():
	index = sample_dict[i]
	print(index)
	print(index1)
	print(i)
	os.system("cutadapt -a {} -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 10 -e 0.1 -O 5 -m 13 -o {}.trimmed.R1.fq.gz -p {}.trimmed.R2.fq.gz {}.R1.fastq.gz {}.R2.fastq.gz".format(index,i,i,i,i))
f.close()

os.system("ls *.trimmed.R1.fq.gz | while read id; do(hisat2 --pen-noncansplice 1000000  -x /home/l/backup1/refgenome/homo_sapiens/rRNA_tRNA/rRNA_tRNA -1 $id -2 $(basename $id '.R1.fq.gz').R2.fq.gz -S $(basename $id '.trimmed.R1.fq.gz').sam); done")
for i in sample_dict.keys():
        os.system("samtools view -f 12 {}.sam > {}.unaligned.sam".format(i,i))
os.system("sh unalign.sh")

#map to known_gene from UCSC
for i in sample_dict.keys():
        os.system("hisat2 -p 30 -x /home/l/backup1/refgenome/homo_sapiens/CDS/CDS_3UTR/18bp/CDS_3UTR_18 -f -U {}.unaligned.fasta -S {}.UCSC.hisat2.sam".format(i,i))
#os.system("bowtie2 -N 1 -x /home/l/backup1/refgenome/homo_sapiens/CDS/CDS_3UTR/18bp/CDS_3UTR_18 -f -U $_.unaligned.fasta -S $_.UCSC.bowtie2.sam")
        os.system("samtools view -F 12 {}.UCSC.hisat2.sam > {}.mapped.UCSC.sam".format(i,i))
        os.system("cat {}.mapped.UCSC.sam | grep NH:i:1 > {}.mapped_1site.UCSC.sam".format(i,i))




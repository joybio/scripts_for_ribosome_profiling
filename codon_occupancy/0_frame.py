#!/bin/miniconda3/bin/python
#coding=utf-8

"""extract unmapped reads from sam"""

__date__ = "2019-10-2"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import re
import os
import optparse
import bisect
#sort  package
from optparse import OptionParser

parser = OptionParser('Usage: %prog ')
parser.add_option('-s','--sam',
                dest='sam',
                help='mapped 1site of UCSC.format.fa in sam file')

parser.add_option('-o','--out',
                dest='out',
                help='generate bed file of 0_frame')

parser.add_option('-f','--fa',
                dest='fa',
                help='fa file')

parser.add_option('-t','--temp',
                dest='temp',
                help='Asite seq in temp.csv')

parser.add_option('-c','--csv',
                dest='csv',
                help='out csv')
(options,args) = parser.parse_args()
data = open("/home/l/backup1/refgenome/homo_sapiens/CDS/CDS_3UTR/18bp/CDS_3UTR_18.human.fa","r")
data_dict = {}
for i in data:
	if i.startswith(">"):
		key = i.lstrip(">").strip()
	else:
		data_dict[key] = len(i)

sam = open(options.sam,'r')
sam_dict = {}
fastq_dict = {}
out = open(options.out,'w')

for i in sam:
	i = i.strip().split('\t')
	n = len(i[9])
	m = int(i[3])
	gene_id = i[2]
	if m+n < data_dict[gene_id]:
		l = m+n
	else:
		l = data_dict[gene_id]-1	
	if m > 15 and n >= 28 and n <= 31:
		if m%3 == 1:
#		bedtools getfasta 提取bed文件start的后一位，因此要减一
#		例如：sam文件里2414是NTGTTGTGAGTTGCAGAGACAGGAAAG。
#		但是bedtools getfasta提取出来则是TGTTGTGAGTTGCAGAGACAGGAAAGATAGAAGACGTTCCAT
#		相当于从2015开始提取
#		再例如：bedtools提取CDS的1-15位，提取的是TGTCGTCCGTTTCA，而不是ATGTCGTCCGTTTCA
#		再例如：bedtools提取CDS的1-28位，提取的是TGTCGTCCGTTTCAGGAAATGCCTCTT，总长为27
#		所以如果余数为1则应该向左补1，即m-1；如果余数为2则应该向右补1，即m+1；如果余数为0，则不需增减
			out.write(i[2] + '\t' + str(m-1) + '\t' + str(l) + '\n')# + \
#				i[9] + "\n")
		elif m%3 == 2:
			out.write(i[2] + '\t' + str(m+1) + '\t' + str(l) + '\n')# + \
#				"N" + i[9] + "\n")
		elif m%3 == 0:
			out.write(i[2] + '\t' + str(m) + '\t' + str(l) + '\n')# + \
#				i[9][1:] + "\n")
sam.close()
out.close()

os.system("bedtools getfasta -fi /home/l/backup1/refgenome/homo_sapiens/CDS/CDS_3UTR/18bp/CDS_3UTR_18.human.fa -bed {} -fo {}".format(options.out,options.fa))

Asitecsv = open(options.csv,'w')
temp = open(options.temp,'w')
with open(options.fa,'r') as f:
	for i in f:
		i = i.strip()
		if i.startswith('>'):
			pass
		else:
			i = i.upper()
			i = list(i)
			if len(i) >= 26:
				temp.write(i[15] + i[16] + i[17] + '\t' \
					+ i[18] + i[19] + i[20] + '\t' \
					+ i[21] + i[22] + i[23] + '\t' \
					+ i[24] + i[25] + i[26] + '\t' \
					+ i[12] + i[13] + i[14] + '\t' \
					+ i[9] + i[10] + i[11] + '\t' \
					+ i[6] + i[7] + i[8] + '\t' \
					+ i[3] + i[4] + i[5] + '\n')
				
f.close()
temp.close()

Asitecsv.write("Asite_codon" + "\t" + "Asite_freq" + "\t" + "Asite1_freq" + "\t" + "Asite2_freq" + "\t" \
	+ "Asite3_freq" + "\t" + "Psite_freq" + "\t" + "Esite_freq" + "\t" + "Esite_1_freq" + "\t" \
	+ "Esite_2_freq" + "\n")

Asite = []
Asite1 = []
Asite2 = []
Asite3 = []
Psite = []
Esite = []
Esite_1 = []
Esite_2 = []

def all_list(arr):
	result = {}
	for i in set(arr):
        	result[i] = str(arr.count(i))
	return result

with open(options.temp,'r') as f:
	for i in f:
		i = i.strip()
		if re.search("N",i):
			pass
		else:
			i = i.split("\t")
			Asite.append(i[0])
			Asite1.append(i[1])
			Asite2.append(i[2])
			Asite3.append(i[3])
			Psite.append(i[4])
			Esite.append(i[5])
			Esite_1.append(i[6])
			Esite_2.append(i[7])

f.close()

A = all_list(Asite)
A1 = all_list(Asite1)
A2 = all_list(Asite2)
A3 = all_list(Asite3)
P = all_list(Psite)		
E = all_list(Esite)
E1 = all_list(Esite_1)
E2 = all_list(Esite_2)

for i in A.keys():
	if i not in A1.keys():
		A1[i] = str("0")
	if i not in A2.keys():
                A2[i] = str("0")
	if i not in A3.keys():
                A3[i] = str("0")
	if i not in P.keys():
                P[i] = str("0")
	if i not in E.keys():
                E[i] = str("0")
	if i not in E1.keys():
                E1[i] = str("0")
	if i not in E2.keys():
                E2[i] = str("0")
	Asitecsv.write(i + '\t' + A[i] + '\t' + A1[i] + '\t' + A2[i] + '\t' + A3[i] + '\t' + P[i] + '\t' + E[i] + '\t' + E1[i] + '\t' + E2[i] + '\n') 
Asitecsv.close()



		



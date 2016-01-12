import re
import os
import sys

# fastqc_parse_v1.0.py
# Version: 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: ##------ Fri Nov  6 14:27:25 2015 ------ #
# Last Modified: ##------ Wed Nov 11 16:48:00 2015 ------##
# Usage: python fastqc_parse_v1.0.py <path to project directory>
# Function: This script goes in the specified project folder and summarizes the fastqc output from all the samples
# It creates an output directory at the project level to store the summary

rootdir = sys.argv[1]

# make output directory
outdir = rootdir+'/output'
if not os.path.exists(outdir):
	os.makedirs(outdir)

# usage python fastqcParse.py <directory with fastqc output files/directories>
with open(os.path.join(outdir,'fastqc_summary.csv'),'w') as fout:
	a = (','.join(('Filename','Encoding','Total Sequences','Sequence length','Per base sequence quality','Per sequence quality','Per base sequence content','Per base GC content','Per sequence GC content','Per base N content','Sequence Duplication levels','Total Duplicate Percentage','Overrepresented','Kmer content')))
	fout.write(''.join((a,'\n')))
  	for root, subFolders, files in os.walk(rootdir):
		for file in files:
			if (file == 'fastqc_data.txt'):
				#with open('fastqc_data.txt','r') as fin:
				with open(os.path.join(root, file), 'r') as fin:
					for lines in fin:
						line = lines.rstrip('\n')
						if (line.startswith('Filename')):
							aa,aavalue = line.split('\t',1)
							avalue = re.sub("_rnaseqFastqInput*","",aavalue)
							# avalue,blah,blaaaa = aavalue.split('.',3)
							# avalue = re.sub(".txt.bz2","",aavalue) for txt.bz2 files
						elif (line.startswith('Encoding')):
							bb,bvalue = line.split('\t',1)
						elif (line.startswith('Total Sequences')):
							cc,cvalue = line.split('\t',1)
						elif (line.startswith('Sequence length')):
							nn,nvalue = line.split('\t',1)
						elif (line.startswith('>>Per base sequence quality')):
							dd,dvalue = line.split('\t',1)	
						elif (line.startswith('>>Per sequence quality')):
							ee,evalue = line.split('\t',1)
						elif (line.startswith('>>Per base sequence content')):
							ff,fvalue = line.split('\t',1)
						elif (line.startswith('>>Per base GC content')):
							gg,gvalue = line.split('\t',1)
						elif (line.startswith('>>Per sequence GC content')):
							hh,hvalue = line.split('\t',1)	
						elif (line.startswith('>>Per base N content')):
							ii,ivalue = line.split('\t',1)
						elif (line.startswith('>>Sequence Duplication Levels')):
							jj,jvalue = line.split('\t',1)
						elif (line.startswith('#Total Duplicate Percentage')):
							kk,kvalue = line.split('\t',1)
						elif (line.startswith('>>Overrepresented')):
							ll,lvalue = line.split('\t',1)
						elif (line.startswith('>>Kmer Content')):
							mm,mvalue = line.split('\t',1)
					z = (','.join((avalue,bvalue,cvalue,nvalue,dvalue,evalue,fvalue,gvalue,hvalue,ivalue,jvalue,kvalue,lvalue,mvalue)))
					y = (''.join((z.split())))
					fout.write(''.join((y,'\n')))

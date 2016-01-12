# getgenecoord.R
# Version: 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: ##------ Fri Nov  6 14:27:25 2015 ------ #
# Last Modified: ##------ Wed Nov 11 16:48:00 2015 ------##
# Usage: Rscript getgenecoord.R convertGTF2txt_v1.0.pl-output.txt
# Function: This script is called by convertGTF2txt_v1.0.pl to convert GTF file format to an easy to use tab-delimited format.
# The output of this script can be used to annotate genes.

# get gene level coordinates
library(plyr)
args <- commandArgs(trailingOnly = TRUE)
filename = args[1];
gtf = read.delim(filename,header=F)
gtf = ddply(gtf, .(V1,V4,V5,V6,V7), summarise, start= min(V2), end= max(V3))
if(length(grep('gencode',filename))){
gtf = gtf[,c(1,6,7,3,2,5,4)]
gtf$V5 = sub('[.][0-9]{1,2}','',gtf$V5)
} else
gtf = gtf[,c(1,6,7,3,2,4,5)]
colnames(gtf) <- c('Chr','Start','End','ID','Strand','Symbol','Biotype')
write.table(gtf,filename,quote=F,row.names=F,col.names=T,sep='\t')


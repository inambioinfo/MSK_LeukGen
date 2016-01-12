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


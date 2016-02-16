# ERCC92toGTF.R
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: 2015-11-13_11:52:36
# Last Modified: 2015-11-13_11:52:36
# Usage: Rscript ERCC92toGTF.R /path/to/ercc92.gtf /path/to/output.gtf
# Function: This script converts the publicly available ERCC92 gtf to ensembl and gencode format

arg <- commandArgs(trailingOnly=T)
input <- arg[1] # input ERCC92 control gtf
output <- arg[2] # output file name

# read default ERCC92 gtf file
dat <- read.table(input,sep=" ")

# edit columns to follow ensembl format
# edit the gene id column
dat$V2 <- sub(";","",dat$V2)
dat$V2 <- paste0("\"",dat$V2,"\";")

# edit the transcript id column
dat$V4 <- sub(";","",dat$V4)
dat$V4 <- paste0("\"",dat$V4,"\";")

# create a gene name and gene biotype column
# change biotype to type for gencode
dat$V5 <- 'gene_name'
dat$V6 <- dat$V2
dat$V7 <- 'gene_biotype'
dat$V8 <- paste0("\"","spikein","\";")

# write the gtf out
write.table(x = dat, file = output, row.names = F, col.names = F, quote = F)

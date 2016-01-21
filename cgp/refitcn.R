# refitcn.R
# Author: Komal Rathi, Gunes Gundem
# Institute: Memorial Sloan Kettering Cancer Center
# Created ##------ Thu Jan 14 10:13:54 2016 ------##
# Last Modified ##------ Thu Jan 14 10:13:54 2016 ------##
# Function: 
# This script takes in the chromosome, position, expected Major and Minor allele copy numbers
# It re-calculates rho and psi and re-fits them to get a new copy number profile

# get environment variables for required paths
TMPBBG <- Sys.getenv('TMPBBG') # tmpBattenberg
TNAME <- Sys.getenv('TNAME') # tumor name
OUTDIR <- Sys.getenv('OUTDIR') # output dir
PROGRESS <- Sys.getenv('PROGRESS') # progress dir

# get environment variables for chr, pos, major and minor allele
chr <- as.numeric(Sys.getenv('CHR')) # chr
pos <- as.numeric(Sys.getenv('POS')) # pos
ref_seg_nMajor <- as.numeric(Sys.getenv('MJA')) # major allele
ref_seg_nMinor <- as.numeric(Sys.getenv('MNA')) # minor allele

# get subclones.txt
dat <- read.delim(paste(TMPBBG,'/',TNAME,'_subclones.txt',sep=''))

# get BAF and LogR of the chr and position
ref_seg_BAF <- dat[dat$chr==chr & dat$startpos==pos,'BAF']
ref_seg_LogR <- dat[dat$chr==chr & dat$startpos==pos,'LogR']

# run script to adjust rho and psi
gamma_param = 1 # 0.55 for array data

# get corrected rho and psi
rho = (2*ref_seg_BAF - 1) / (2*ref_seg_BAF - ref_seg_BAF*(ref_seg_nMajor + ref_seg_nMinor) - 1 + ref_seg_nMajor)
psi = (rho*(ref_seg_nMajor + ref_seg_nMinor) + 2 - 2*rho)/(2^(ref_seg_LogR/gamma_param))
psi = (psi-2*(1-rho))/rho

# read err file to get the command
library(stringr)
f <- readLines(paste(OUTDIR,'/err/bbg_fitcn.1.err',sep=''))
cmd <- grep("Errors from command:",f,value=TRUE)
cmd <- sub('Errors from command: ','',cmd)

# add rho and psi to cmd
cmd <- paste(cmd,rho,psi) 

# call the command
exitcode <- system(cmd)

# if the command runs succesfully, make a progress file
if(exitcode==0){
  file <- paste(PROGRESS,'/Sanger::CGP::Battenberg::Implement::battenberg_refitcn.0',sep='')
  file.create(file=file)
}

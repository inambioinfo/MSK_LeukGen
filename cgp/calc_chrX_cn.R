# calc_chrX_cn.R
# Author: Komal Rathi, Gunes Gundem
# Institute: Memorial Sloan Kettering Cancer Center
# Created ##------ Mon Jan 18 10:07:27 2016 ------##
# Last Modified ##------ Mon Jan 18 10:07:27 2016 ------##
# Function:
# This script calculates copy number for chrX since Battenberg does not call CN for chrX
# Usage: Rscript /ifs/work/leukgen/local/bin/r/segment_chrX.R

# get environment variables for required paths
TNAME <- Sys.getenv('TNAME') # tumor sample name
OUTDIR <- Sys.getenv('OUTDIR') # output dir
is_male = Sys.getenv('SEX') # -is-male or empty
BBGLIB <- Sys.getenv('BBGLIB') # Battenberg lib directory

# determine input and output files
input_file = paste(OUTDIR,'/',TNAME,'_mutantLogR.tab',sep='') # input file
samplename = as.character(TNAME) # sample name
output_file = paste(OUTDIR,'/',TNAME,'_caveman_cn_file_chrX.txt',sep='') # output file
print(output_file)

# get psi and rho from OUTDIR/TNAME_rho_and_psi.txt
# use FRAC_GENOME row to extract the values
dat <- read.delim(paste(OUTDIR,'/',TNAME,'_rho_and_psi.txt',sep=''))
ploidy = dat[2,2]
cellularity = dat[2,1]

# source fastPCF.R
source(paste(BBGLIB,'/fastPCF.R',sep=''))
input_dir = dirname(input_file)

if(is_male=="-is-male") {
  is_male = TRUE
} else {
  is_male = FALSE
}

logr.chrx.file = paste(getwd(),'/', samplename, "_mutantLogR_chrX.tab", sep="")
system(paste("grep Chromosome ", input_file, " > ", logr.chrx.file, "; zless ", input_file, " | awk '$2==\"X\"' >> ", logr.chrx.file, sep=" "))
gamma = 25
data = read.table(logr.chrx.file,sep="\t",header=T)
data = data[!is.na(data[,3]) & !is.nan(data[,3]) & !is.infinite(data[,3]),]
sdev = getMad(data[,3],k=25)
res= selectFastPcf(data[,3],6,gamma*sdev,T)
segLogR = res$yhat
out.data = cbind(data,segLogR)

logRadjustment = log2(ploidy/2)

out.data$CN = (2^(out.data$segLogR+logRadjustment)-1 + cellularity)/cellularity
unique.CN = unique(out.data$CN)
min.max.pos = array(NA,c(length(unique.CN),2))
for(i in 1:length(unique.CN)){
  min.max.pos[i,1] = min(out.data$Position[out.data$CN==unique.CN[i]])
  min.max.pos[i,2] = max(out.data$Position[out.data$CN==unique.CN[i]])
}
ord = order(min.max.pos[,1])

#avoid negative copy numbers! - now moved down
out.data$CN[out.data$CN<0] = 0
unique.CN[unique.CN<0] = 0

out.data2 = cbind(min.max.pos[ord,],round(unique.CN[ord]))
out.data3 = array(NA,c(0,3))
start.seg = -1
end.seg = -1
CN = -10000000
for(i in 1:nrow(out.data2)){
  if(out.data2[i,3]!=CN){
    if(start.seg != -1){
      out.data3 = rbind(out.data3,c(start.seg,end.seg,CN))
    }
    start.seg = out.data2[i,1]              
    CN = out.data2[i,3]
  }
  end.seg = out.data2[i,2]
}
if(is_male) {
  write.table(cbind("X", out.data3[,1:2],1,0,out.data3[,3],0), output_file, sep="\t",quote=F,row.names=F,col.names=F)
} else {
  write.table(cbind("X", out.data3[,1:2],2,0,out.data3[,3],0), output_file, sep="\t",quote=F,row.names=F,col.names=F)
}
system(paste('rm ', logr.chrx.file, sep=''))

# simulate_copynumber.R
# Author: Komal Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created ##------ Tue Mar 15 18:02:11 2016 ------##
# Last Modified ##------ Tue Mar 15 18:02:11 2016 ------##
# Function:
# This script simulates copy number for caveman
# Usage: Rscript /ifs/work/leukgen/local/bin/r/simulate_copynumber.R

TCP <- Sys.getenv('TCP') # tumor copy number
NCP <- Sys.getenv('NCP') # normal copy number
CAVOUT <- Sys.getenv('CAVOUT') # caveman output dir
TNAME <- Sys.getenv('TNAME') # tumor name

tumour_out_file = paste(CAVOUT,'/',TNAME,'_caveman_tumour_cn_input.txt',sep='') # input file
normal_out_file = paste(CAVOUT,'/',TNAME,'_caveman_normal_cn_input.txt',sep='') # input file

dat <- structure(list(V1 = structure(1:23, .Label = c("1", "10", "11", "12", "13", "14", "15", "16", 
                                                      "17", "18", "19", "2", "20", "21", "22", "3", 
                                                      "4", "5", "6", "7", "8", "9", "X"), class = "factor"), 
                      V2 = c(762601L, 94426L, 218141L, 188285L, 19455957L, 20433516L, 
                      20001861L, 84170L, 3104L, 125371L, 282753L, 10797L, 20000786L, 
                      15273232L, 16857983L, 60363L, 11870L, 849605L, 147750L, 80500L, 
                      159286L, 203761L, 60112L), 
                      V3 = c(249201268L, 135226330L, 134946396L, 133839356L, 114938134L, 107260733L, 102377717L, 
                      89998957L, 79998834L, 78017073L, 59097308L, 242985695L, 62917667L, 
                      48092076L, 49992906L, 197846280L, 190790246L, 180715810L, 
                      170919682L, 159122659L, 146300622L, 140969861L, 155177250L)), 
                 .Names = c("V1", "V2", "V3"), row.names = c(NA, -23L), class = "data.frame")

dat$V4 <- TCP
write.table(dat, tumour_out_file, sep="\t", quote=F, row.names=F, col.names=F)

dat$V4 <- NCP
write.table(dat, normal_out_file, sep="\t", quote=F, row.names=F, col.names=F)
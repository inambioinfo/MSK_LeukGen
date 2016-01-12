library(DESeq)

# function to do differential expression
source('DECompare.R') 

# commandline arguments
# argument 1 = directory where counts are stored
# argument 2 = sample info file
args <- commandArgs(TRUE)

# get command line arguments
dir <- args[1]
sample.info <- args[2]
sample.info <- read.csv(sample.info, stringAsFactors = F)

# get counts
dat.list <- list.files(path=dir,pattern=".counts",full.names=T)
nms <- sub('_htseq_.*','',sub('.*/','',dat.list))

dat <- do.call('cbind',lapply(dat.list,read.delim,header=F))
row.names(dat) <- dat$V1
dat <- dat[,grep("V2", colnames(dat))]
dat <- dat[grep('^ENS',rownames(dat)),] # remove unwanted features

# based on sample info get the treatment info
condition <- sample.info[match(nms, sample.info$sample_name),'sample_source']
names(dat) <- condition

# call DECompare
res <- DECompare(dat, 'tumor', 'normal', condition)
res.sig <- res[which(res$padj < 0.05),] # only get genes with pvalues < 0.05

# create boxplots

# create PCA plot

# create heatmaps
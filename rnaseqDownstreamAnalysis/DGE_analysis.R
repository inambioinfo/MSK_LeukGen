# DGE_analysis.R
# Version: 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: ##------ Wed Nov 18 17:30:07 2015 ------##
# Last Modified: ###------ Tue Nov 24 12:51:14 2015 ------##
# Function: This script takes count data and produces the following files in the output directory:
# 1. csv file of normalized counts
# 2. csv file of limma results
# 3. boxplot of normalized counts
# 4. pca of top 1% genes 
# 5. heatmap of top 100 differentially expressed genes
# Usage: Rscript DGE_analysis.R <input_directory> <sample information> <path to gtf file>

library(limma)
library(grid)
library(gridExtra)

# read arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- toString(args[1]) # project directory
sample.info <- toString(args[2]) # sample information
annotation <- toString(args[3]) # gtf annotation
sample.info <- read.csv(sample.info)

# get required scripts
source('/Users/rathik/scripts/voomApply.R')
source('/Users/rathik/scripts/limmaApply.R')
source('/Users/rathik/scripts/makePCA.R')
source('/Users/rathik/scripts/makeBoxplot.R')
source('/Users/rathik/scripts/makeCorrPlot.R')
source('/Users/rathik/scripts/makeHeatmap.R')
source('/Users/rathik/scripts/DECompare.R')

# get counts
dat.list <- list.files(path=input_dir,pattern="htseq_ensembl.counts",full.names=T)  
nms <- sub('_rnaseqFastqInput_.*','',sub('.*/','',dat.list))

# subset dat.list for samples that you want to test
dat.list <- dat.list[nms %in% sample.info$Source_name]

dat <- do.call('cbind',lapply(dat.list,read.delim,header=F))
row.names(dat) <- dat$V1
dat <- dat[,grep("V2", colnames(dat))]
dat <- dat[grep('^ENS',rownames(dat)),] # remove unwanted features

# pre normalization filtering
# remove all genes have low counts before applying voom normalization
# dat <- dat[apply(dat, MARGIN = 1,FUN = function(x) all(x>0)),] 
# remove all genes that have total read count across all samples < 10
dat <- dat[rowSums(dat) > 10,]

# set conditions
condition <- sample.info[na.omit(match(nms, sample.info$Source_name)),'Source_type']

# set names
names(dat) <- sample.info[na.omit(match(nms, sample.info$Source_name)),'Source_name']

# perform normalization
normCounts <- voomApply(counts = dat, groupA = 'Tumor', groupB = 'Normal', conditions = condition, sample.info = sample.info)
normCounts <- normCounts - min(normCounts)

# add pseudo counts to normalized counts i.e. smallest non-zero count in each library
# log2 transform

# generate correlation plot of the normalized samples
corrplot <- makeCorrPlot(normCounts = normCounts)
corrplot

# generate PCA of samples
pca <- makePCA(normCounts = normCounts, sample.info = sample.info)
pca

# generate boxplot of samples
boxplot <- makeBoxplot(normCounts = normCounts, sample.info = sample.info)
boxplot

# get gene annotation
annotation <- read.delim(annotation)

# compute differential gene expression
res.limma <- limmaApply(normCounts = normCounts, groupA = 'Tumor', groupB = 'Normal', conditions = condition, sample.info = sample.info, annotation = annotation)
write.csv(res.limma,'data/VoomLimma_results.csv',quote=F,row.names = F)

res.deseq <- DECompare(counts = dat, groupA = 'Tumor', groupB = 'Normal', conditions = condition, sample.info = sample.info, annotation = annotation)
write.csv(res.deseq,'data/DESeq_results.csv',quote=F,row.names = F)

# p-value distribution P-value vs. Frequency

# generate heatmap
heatmap <- makeHeatmap(results = res.limma, normCounts = normCounts, sample.info = sample.info, title = 'Top Differentially Expressed Genes\n', n = 50)

# output results
output <- paste(input_dir,'DGE_summary_report.pdf',sep='/')
pdf(file = output, width = 25, height = 20)
grid.arrange(arrangeGrob(main=textGrob('Plots: MPL vs MigR1 CD11b',gp = gpar(fontsize=20,fontface='bold')),
                         pca,boxplot,heights = c(1,4,4)),heatmap$gtable,ncol=2)
invisible(dev.off())

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

# load libraries
library(limma)
library(grid)
library(gridExtra)
library(plyr)

# get required scripts
source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/voomApply.R')
source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/limmaApply.R')
source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/makePCA.R')
source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/makeBoxplot.R')
source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/makeCorrPlot.R')
source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/makeHeatmap.R')
source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/DECompare.R')

# read arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- toString(args[1])   # project directory
sample.info <- toString(args[2]) # sample information
annotation <- toString(args[3])  # gtf annotation
conf <- toString(args[4])        # config file

# read sample info
sample.info <- read.csv(sample.info)

# read config file
conf <- read.csv(conf, stringsAsFactors = F)

# get gene annotation
annotation <- read.delim(annotation)

# get counts
dat.list <- list.files(path=input_dir,pattern=".counts",full.names=T)
nms <- sub('_rnaseqFastqInput_.*','',sub('.*/','',dat.list))

# subset dat.list for samples that you want to test
dat.list <- dat.list[nms %in% sample.info$source_name]

dat <- do.call('cbind',lapply(dat.list,read.delim,header=F))
row.names(dat) <- dat$V1
dat <- dat[,grep("V2", colnames(dat))]
dat <- dat[grep('^ENS',rownames(dat)),] # remove unwanted features

# pre normalization filtering
# remove all genes have low counts before applying voom normalization
# dat <- dat[apply(dat, MARGIN = 1,FUN = function(x) all(x>0)),]
# remove all genes that have total read count across all samples < 10
dat <- dat[rowSums(dat) > 10,]
rownames(dat) <- sub('[.].*','',rownames(dat))

# set conditions
condition <- sample.info[na.omit(match(nms, sample.info$source_name)),'source_type']

# set names
names(dat) <- sample.info[na.omit(match(nms, sample.info$source_name)),'source_name']

for(i in 1:nrow(conf)){
  groupA <- conf[i,1]
  groupB <- conf[i,2]
  print(paste(groupA,"vs",groupB))

  # perform normalization
  normCounts <- voomApply(counts = dat, groupA = groupA, groupB = groupB, conditions = condition, sample.info = sample.info)
  normCounts <- normCounts - min(normCounts)

  # add pseudo counts to normalized counts i.e. smallest non-zero count in each library
  # log2 transform

  # generate correlation plot of the normalized samples
  corrplot <- makeCorrPlot(normCounts = normCounts)

  # generate PCA of samples
  pca <- makePCA(normCounts = normCounts, sample.info = sample.info, m=10)

  # generate boxplot of samples
  boxplot <- makeBoxplot(normCounts = normCounts, sample.info = sample.info)

  # compute differential gene expression
  res.limma <- limmaApply(normCounts = normCounts, groupA = groupA, groupB = groupB, conditions = condition, sample.info = sample.info, annotation = annotation)
  output <- paste(input_dir,paste(groupA,"vs",groupB,'VoomLimma_results.csv',sep = '_'),sep='/')
  write.csv(res.limma, output, quote=F, row.names = F)

  res.deseq <- DECompare(counts = dat, groupA = groupA, groupB = groupB, conditions = condition, sample.info = sample.info, annotation = annotation)
  output <- paste(input_dir,paste(groupA,"vs",groupB,'DESeq_results.csv',sep = '_'),sep='/')
  write.csv(res.deseq, output, quote=F, row.names = F)

  # p-value distribution P-value vs. Frequency

  # generate heatmap
  heatmap <- makeHeatmap(results = res.limma, normCounts = normCounts, sample.info = sample.info, title = 'Top Differentially Expressed Genes\n', n = 50)

  # output results
  output <- paste(input_dir,paste(groupA,"vs",groupB,'DGE_summary_report.pdf',sep = '_'),sep='/')
  text <- paste('Plots:',groupA,'vs',groupB,sep=' ')
  pdf(file = output, width = 25, height = 20)
  grid.arrange(arrangeGrob(main=textGrob(text,gp = gpar(fontsize=20,fontface='bold')),
                           pca,boxplot,heights = c(1,4,4)),heatmap$gtable,ncol=2)
  invisible(dev.off())
}

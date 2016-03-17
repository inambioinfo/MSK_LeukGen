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
suppressMessages(library(limma, quietly = TRUE))
suppressMessages(library(grid, quietly = TRUE))
suppressMessages(library(gridExtra, quietly = TRUE))
suppressMessages(library(plyr, quietly = TRUE))

# get required scripts
source('/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/r/voomApply.R')
source('/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/r/limmaApply.R')
source('/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/r/makePCA.R')
source('/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/r/makeBoxplot.R')
source('/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/r/makeCorrPlot.R')
source('/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/r/makeHeatmap.R')
source('/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/r/DECompare.R')
source('/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/r/rgsea.R')

# read arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- toString(args[1])               # project directory
sample.info <- toString(args[2])             # sample information
annotation <- toString(args[3])              # gtf annotation
conf <- toString(args[4])                    # config file
npca <- as.integer(args[5])                  # top percent genes in pca
nheatmap <- as.integer(args[6])              # number of genes in heatmap
sample.ann <- toString(args[7])              # sample annotation containing species name

# read sample info
sample.info <- read.csv(sample.info)

# read config file
conf <- read.csv(conf, stringsAsFactors = F)

# read sample annotation file
sample.ann <- read.csv(sample.ann, stringsAsFactors = F)

# get species
species <- unique(sample.ann$species)

# get gene annotation
if(tolower(species)=="human"){
    if(tolower(annotation)=="ensembl")
    {
        annotation <- '/ifs/work/leukgen/ref/homo_sapiens/37/ensembl/75/Homo_sapiens.GRCh37.75.txt'
    }
    if(tolower(annotation)=="gencode"){
        annotation <- '/ifs/work/leukgen/ref/homo_sapiens/37/gencode/18/gencode.v18.annotation.txt'
    }
}
if(tolower(species)=="mouse"){
    if(tolower(annotation)=="ensembl")
    {
        annotation <- '/ifs/work/leukgen/ref/mus_musculus/38/ensembl/80/Mus_musculus.GRCm38.80_canonical_chromosomes_nochr.txt'
    }
    if(tolower(annotation)=="gencode"){
        annotation <- '/ifs/work/leukgen/ref/mus_musculus/38/gencode/M7/gencode.vM7.annotation.txt'
    }
}

# read annotation
annotation <- read.delim(annotation)

# get counts
# dat.list <- list.files(path=input_dir,pattern=".counts",full.names=T)
dat.list <- list.files(path=input_dir,pattern="_ReadsPerGene.out.tab",full.names=T)
nms <- sub('_.*','',sub('.*/','',dat.list))

# subset dat.list for samples that you want to test
# dat.list <- dat.list[nms %in% sample.info$source_name]**
dat.list <- dat.list[nms %in% sample.ann$source_name]

dat <- do.call('cbind',lapply(dat.list,read.delim,header=F))
row.names(dat) <- dat$V1
dat <- dat[,grep("V2", colnames(dat))]
dat <- dat[grep('^ENS',rownames(dat)),] # remove unwanted features

# pre normalization filtering
# remove all genes have low counts before applying voom normalization
# atleast 1 sample should have >20 reads
dat <- dat[apply(dat, MARGIN = 1, function(x) any(x > 20)),]
rownames(dat) <- sub('[.].*','',rownames(dat))

# set conditions
# condition <- sample.info[na.omit(match(nms, sample.info$source_name)),'source_condition']

# set names
# names(dat) <- sample.info[na.omit(match(nms, sample.info$source_name)),'source_name']**
names(dat) <- sample.ann[na.omit(match(nms, sample.ann$source_name)),'ext_id']

for(i in 1:nrow(conf)){
  groupA <- conf[i,1]
  groupB <- conf[i,2]
  print(paste(groupA,"vs",groupB))

  # perform normalization
  # add pseudo counts to normalized counts i.e. smallest non-zero count in each library
  normCounts <- voomApply(counts = dat, sample.info = sample.info)
  normCounts <- normCounts - min(normCounts)

  # generate correlation plot of the normalized samples
  corrplot <- makeCorrPlot(normCounts = normCounts)
  output <- paste(input_dir,paste(groupA,"vs",groupB,'Corrplot.pdf',sep = '_'),sep='/')
  pdf(file = output, width = 15, height = 10)
  print(corrplot)
  invisible(dev.off())

  # generate PCA of samples
  pca <- makePCA(normCounts = normCounts, sample.info = sample.info, m = npca)
  output <- paste(input_dir,paste(groupA,"vs",groupB,'PCA.pdf',sep = '_'),sep='/')
  pdf(file = output, width = 15, height = 10)
  print(pca)
  invisible(dev.off())

  # generate boxplot of samples
  boxplot <- makeBoxplot(normCounts = normCounts, sample.info = sample.info)
  output <- paste(input_dir,paste(groupA,"vs",groupB,'Boxplot.pdf',sep = '_'),sep='/')
  pdf(file = output, width = 15, height = 10)
  print(boxplot)
  invisible(dev.off())

  # compute differential gene expression
  res.limma <- limmaApply(normCounts = normCounts, groupA = groupA, groupB = groupB, sample.info = sample.info, annotation = annotation)
  output <- paste(input_dir,paste(groupA,"vs",groupB,'VoomLimma_results.csv',sep = '_'),sep='/')
  write.csv(res.limma, output, quote=F, row.names = F)

  if(length(groupA)==1 && length(groupB)==1){
    res.deseq <- DECompare(counts = dat, groupA = groupA, groupB = groupB, sample.info = sample.info, annotation = annotation)
    output <- paste(input_dir,paste(groupA,"vs",groupB,'DESeq_results.csv',sep = '_'),sep='/')
    write.csv(res.deseq, output, quote=F, row.names = F)
  }

  # gene set enrichment analysis
  rgsea(input_dir, normCounts, annotation, sample.info, res.limma, groupA, groupB)

  # generate heatmap
  heatmap <- makeHeatmap(results = res.limma, normCounts = normCounts, sample.info = sample.info, title = 'Top Differentially Expressed Genes\n', n = nheatmap)
  output <- paste(input_dir,paste(groupA,"vs",groupB,'Heatmap.pdf',sep = '_'),sep='/')
  pdf(file = output, width = 15, height = 10)
  grid.arrange(heatmap$gtable,ncol=1)
  invisible(dev.off())

  # output results
  output <- paste(input_dir,paste(groupA,"vs",groupB,'DGE_summary_report.pdf',sep = '_'),sep='/')
  text <- paste('Plots:',groupA,'vs',groupB,sep=' ')
  pdf(file = output, width = 25, height = 20)
  grid.arrange(arrangeGrob(main=textGrob(text,gp = gpar(fontsize=20,fontface='bold')),
                           pca,boxplot,heights = c(1,4,4)),heatmap$gtable,ncol=2)
  invisible(dev.off())
}

# variant_correction.R
# Version: 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: ##------ Wed Apr  6 12:06:32 2016 ------##
# Last Modified: ##------ Wed Apr  6 12:06:32 2016 ------##
# Function: This script takes flagged vcf data from pindel output and generates input for and calls variant correction script
# Usage: Rscript variant_correction.R <path to samplelist> <path to output file>

suppressMessages(library(VariantAnnotation, quietly = TRUE))
suppressMessages(library(plyr, quietly = TRUE))
suppressMessages(library(stringi, quietly = TRUE))

# get annotation from INFO
get.info <- function(x){
  tmp <- x$Info

  x$Gene <- gsub('VD=([^|]+)|.', '\\1', tmp)
  x$Transcript <- unlist(strsplit(tmp,"[|]"))[2]
  x$RNA <- unlist(strsplit(tmp,"[|]"))[3]
  x$CDS <- unlist(strsplit(tmp,"[|]"))[4]
  x$Protein <- unlist(strsplit(tmp,"[|]"))[5]
  x$Type <- gsub('VT=([^;]+)|.', '\\1', tmp)
  x$Effect <- gsub('VC=([^;]+)|.', '\\1', tmp)
  x$PC <- gsub('PC=([^;]+)|.', '\\1', tmp)

  return(x)
}

# get arguments
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]  # sample list
output <- args[2] # output file

# read input file
input <- read.delim(input, header=F, stringsAsFactors = F)

tmp <- data.frame()
for(i in 1:nrow(input)){

  samp <- input[i, 1]
  if(length(grep('vcf.gz$',samp))==1){

    # read vcf file
    vcf <- readVcf(samp, "hg19")

    # get tumorname and normalname from vcf header
    vcf.header <- readLines(samp)
    normalname <- vcf.header[grep('ID=NORMAL',vcf.header)]
    normalname <- sub('>','',sub('.*SampleName=','',normalname))
    tumorname <- vcf.header[grep('ID=TUMOUR',vcf.header)]
    tumorname <- sub('>','',sub('.*SampleName=','',tumorname))

    # column names
    col.names <- vcf.header[grep('#CHROM',vcf.header)]
    col.names <- sub('#','',col.names)
    col.names <- unlist(strsplit(col.names, '\t'))

    # read vcf contents
    vcf <- read.table(gzfile(samp), header=F, comment.char = "#", stringsAsFactors = F)
    colnames(vcf) <- stri_trans_totitle(col.names)

    # get only variants that PASS
    vcf <- vcf[vcf$Filter=="PASS",]
    vcf[] <- lapply(vcf, as.character)

    # add annotation
    vcf$Sample <- tumorname
    vcf$Normal <- normalname
    vcf$VariantID <- seq(1:nrow(vcf))
    vcf <- ddply(.data = vcf, .variables = 'Id', .fun = function(x) get.info(x))
    vcf[is.na(vcf)] <- "-"
    vcf[vcf==""] <- "-"
    vcf$AnalysisProc <- NA
    vcf <- vcf[,c("AnalysisProc","Sample", "Normal", "VariantID", "Chrom", "Pos", "Ref", "Alt", "Qual", "Filter", "Gene", "Transcript", "RNA", "CDS", "Protein", "Type", "Effect", "PC")]
    vcf$BAM <- input[i, 2]
    colnames(vcf)[1] <- "#AnalysisProc"
    tmp <- rbind(tmp, vcf)
  }
}

# write merged vcf output
tmp <- tmp[-which(tmp$Chrom=="MT"),]
write.table(x = tmp, file = output, quote=F, row.names = F, col.names = T, sep = "\t")



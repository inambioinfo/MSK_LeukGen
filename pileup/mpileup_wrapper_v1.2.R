# R wrapper script for mpileup -> nucleotide frequencies
# mpileup_wrapper_v1.2.R
# version: 1.2
# Komal S Rathi
# Memorial Sloan Kettering Cancer Center
# Created ##------ Fri Nov  6 14:27:25 2015 ------ #
# Last Modified ##------ Mon Feb  8 11:29:43 2016 ------##
# Usage:
# Rscript mpileup_wrapper.R <input directory or list file> <positions_file> <reference_fasta> <max_depth> <min_base_quality> <min_mapq> <include_insertions> <include_deletions>
# only the first three arguments are mandatory
# Changes in this version:
# takes position range as well as list of positions
# take input directory or file with list of bam files
# writes pileup output for samples with 0 depth
# writes pileup output for positions with 0 coverage
# takes commandline parameters

timestamp()
suppressMessages(library(Rsamtools, quietly = TRUE))
suppressMessages(library(reshape2, quietly = TRUE))
library(optparse)
library(getopt)


option_list = list(make_option(c("-i", "--input"), type="character",
                               help="input directory or sample file list",
                               metavar="character",action="store"),
                   make_option(c("-p", "--pos"), type="character",
                               help="positions file", metavar="character"),
                   make_option(c("-f", "--fasta"), type="character",
                               help="referece fasta path", metavar="path"),
                   make_option(c("-d", "--depth"), type="integer",
                               help="max depth [default %default]",
                               default=1000, metavar="integer"),
                   make_option(c("-bq", "--minbq"), type="integer",
                               help="min base quality [default %default]",
                               default=0, metavar="integer"),
                   make_option(c("-mq", "--minmq"), type="integer",
                               help="min mapq [default %default]",
                               default=0, metavar="integer"),
                   make_option(c("-ins", "--insertions"), type="logical",
                               default=FALSE,
                               help="include insertions? [default %default]",
                               metavar="logical"),
                   make_option(c("-del", "--deletions"), type="logical",
                               default=FALSE,
                               help="include deletions? [default %default]",
                               metavar="logical"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# get argument values
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Provide input file (list of bams) or input directory", call.=FALSE)
}

if (is.null(opt$pos)){
  print_help(opt_parser)
  stop("Provide positions file", call.=FALSE)
}

if (is.null(opt$fasta)){
  print_help(opt_parser)
  stop("Provide reference fasta", call.=FALSE)
}

positions_file <- opt$pos
input <- opt$input
fasta_file <- opt$fasta
max_depth <- opt$depth
min_base_quality <- opt$minbq
min_mapq <- opt$minmq
include_insertions <- opt$insertions
include_deletions <- opt$deletions

# read positions file & bam files
# read positions file & bam files
positions_file <- read.delim(positions_file, header = F)

# if input is directory, process all bams in that directory
# else if input is file, process all bams listed in the text file
if(dir.exists(input)){
    files <- list.files(path = input,
                        pattern = "*.sorted.bam$",
                        full.names = T,
                        recursive = T)
} else if(file.exists(input)){
    files <- read.table(input, comment.char = "")
    files <- as.vector(files$V1)
}

# read fasta reference
fasta_file <- FaFile(file = fasta_file)

# format positions file and get reference base
if(ncol(positions_file) >= 3 && nrow(positions_file) == 1){
    # range file only one range
    tags <- positions_file$V4
    positions_file <- data.frame(V1 = positions_file$V1,
                                 V2 = seq(from = positions_file$V2,
                                          to = positions_file$V3))
    refbase <- getSeq(fasta_file,
                      GRanges(positions_file$V1,
                              IRanges(start = as.numeric(positions_file$V2),
                                      end = as.numeric(positions_file$V2))))

    refbase <- as.data.frame(refbase)$x
    positions_file$REF <- refbase

} else if(ncol(positions_file) >= 3 && nrow(positions_file) > 1){
    # when you have multi range file
    # tags <- positions_file$V4
    mylist <- apply(positions_file,1,function(x) data.frame(V1 = as.numeric(x[1]),V2 = seq(from = as.numeric(x[2]), to = as.numeric(x[3])), V3 = as.character(x[4])))
    positions_file <- do.call(rbind, lapply(mylist, data.frame, stringsAsFactors=FALSE))
    tags <- positions_file$V3
    refbase <- getSeq(fasta_file,
                      GRanges(positions_file$V1,
                              IRanges(start = as.numeric(positions_file$V2),
                                      end = as.numeric(positions_file$V2))))
    refbase <- as.data.frame(refbase)$x
    positions_file$REF <- refbase
} else if(ncol(positions_file) < 3 && nrow(positions_file) > 1){
    # when you have chr pos list
    refbase <- getSeq(fasta_file,
                      GRanges(positions_file$V1,
                              IRanges(start = as.numeric(positions_file$V2),
                                      end = as.numeric(positions_file$V2))))
    refbase <- as.data.frame(refbase)$x
    positions_file$REF <- refbase
}

# get pileup for each file
for(i in files){
    print(paste("Processing...", basename(i), sep = ''))
    bamfile <- i
    bf <- BamFile(bamfile)
    param <- ScanBamParam(which = GRanges(positions_file$V1,
                          IRanges(start = as.numeric(positions_file$V2),
                                  end = as.numeric(positions_file$V2))))

    # change max depth, strand specification, various cut-offs
    p_param <- PileupParam(distinguish_strand = TRUE,
                           distinguish_nucleotides = TRUE,
                           max_depth = max_depth,
                           include_deletions = include_deletions,
                           include_insertions = include_insertions,
                           min_base_quality = min_base_quality,
                           min_mapq = min_mapq,
                           min_nucleotide_depth = 0)

    # call pileup function
    res <- pileup(bf, scanBamParam = param, pileupParam = p_param)

    # if there is no pileup for this file, go to next
    # if(nrow(res) == 0) next
    # if there is no pileup, generate 0 values
    if(nrow(res) == 0){
        res <- data.frame(seqnames=positions_file$V1,pos=positions_file$V2,which_label=paste(positions_file$V1,":",positions_file$V2,"-",positions_file$V2,sep=''))
        res$strand <- rep(c('-','+'),
                          length.out=nrow(positions_file))
        res$nucleotide <- rep(c('A','T','G','C'),
                              length.out=nrow(positions_file))
        res$count <- 0
        res <- res[,c(1,2,4,5,6,3)]
    }

    # get reference base
    res <- merge(res, positions_file, by.x = c('seqnames', 'pos'), by.y = c('V1', 'V2'),all.y=TRUE)
    res$strand <- factor(res$strand,levels=c('-','+'))
    res$nucleotide <- factor(res$nucleotide,levels=c('A','T','G','C'))
    res$seqnames <- factor(res$seqnames, levels=unique(res$seqnames))

    # process and write the output
    results <- dcast(res, seqnames+pos~nucleotide+strand,
                     value.var = "count", fill = 0, drop=FALSE)
    results$REF <- positions_file$REF
    results <- results[,c(1,2,ncol(results),3:(ncol(results)-1))]

    # remove columns that have NA in the name
    if(length(grep("NA$|^NA",colnames(results)))>0){
        results <- results[, -grep("NA$|^NA", colnames(results))]
    }

    # calculate depth
    results$D <- apply(results[, 4:ncol(results)], 1, sum)
    results$D_forward <- apply(results[, grep('[+]', colnames(results))],
                               1, sum)
    results$D_reverse <- apply(results[, grep('[-]', colnames(results))],
                               1, sum)
    results$TAG <- tags
    results <- results[, c(1, 2, 3, ncol(results),
                           grep('D', colnames(results)),
                           grep('[+]', colnames(results)),
                           grep('[-]', colnames(results)))]
    colnames(results) <- sub('[+]', 'forward', colnames(results))
    colnames(results) <- sub('[-]', 'reverse', colnames(results))
    colnames(results)[1:2] <- c('CHR', 'POS')

    # output generation
    outfile <- sub('[.]bam', '.out', sub('.*/', '', i))
    write.table(x = results, file = outfile, quote = F,
                row.names = F, sep = '\t')
}

print('Total time taken...')
time <- proc.time()
print(paste(time[[1]], 'secs', sep = ' '))
timestamp()

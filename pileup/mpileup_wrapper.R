# R wrapper script for mpileup -> nucleotide frequencies
# mpileup_wrapper.R
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: ##------ Fri Nov  6 14:27:25 2015 ------ #
# Last Modified: ##------ Wed Nov 11 16:48:00 2015 ------##
# Function: This script takes in a number of bam files and generates pileups for a given list of positions. 
# The pileups are then converted to easy to read format of nucleotide frequencies.
# Usage:
# Rscript mpileup_wrapper.R <input_dir> <positions_file> <reference_fasta> <max_depth> <min_base_quality> <min_mapq> <include_insertions> <include_deletions>
# only the first three arguments are mandatory

timestamp()
suppressMessages(library(Rsamtools,quietly=TRUE))
suppressMessages(library(reshape2,quietly=TRUE))

# command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- toString(args[1])
positions_file <- toString(args[2])
fasta_file <- toString(args[3])
max_depth <- 1000
min_base_quality <- 0
min_mapq <- 0
include_insertions <- FALSE
include_deletions <- FALSE

if(length(args)>=4){
	max_depth <- as.integer(args[4])
	if(length(args)>=5){
		min_base_quality <- as.integer(args[5])
		if(length(args)>=6){
			min_mapq <- as.integer(args[6])
			if(length(args)>=7){
				include_deletions <- as.logical(args[7])
				if(length(args)>=8){
					include_insertions <- as.logical(args[8])
				}
			}
		}
	}
}

# read positions file & bam files
positions_file <- read.delim(positions_file,header=F)
files <- list.files(path=input_dir, pattern="*.bam$", full.names=T, recursive=T)

# read fasta reference 
# get reference base
fasta_file <- FaFile(file=fasta_file)
refbase <- getSeq(fasta_file,GRanges(positions_file$V1,IRanges(start=as.numeric(positions_file$V2),end=as.numeric(positions_file$V2))))
refbase <- as.data.frame(refbase)$x
positions_file$REF <- refbase

# get pileup for each file
for(i in files){
	print(paste("Processing...",basename(i),sep=''))
	bamfile <- i
	bf <- BamFile(bamfile)
	param <- ScanBamParam(which=GRanges(positions_file$V1,IRanges(start=as.numeric(positions_file$V2),end=as.numeric(positions_file$V2))))

	# change max depth, strand specification, various cut-offs
	p_param <- PileupParam(distinguish_strand=TRUE,distinguish_nucleotides=TRUE,
				max_depth=max_depth,include_deletions=include_deletions,
				include_insertions=include_insertions,min_base_quality=min_base_quality,min_mapq=min_mapq)
	
	# call pileup function
	res <- pileup(bf, scanBamParam=param, pileupParam=p_param)
	
	# get reference base
	res <- merge(res,positions_file,by.x=c('seqnames','pos'),by.y=c('V1','V2'))
	
	# process and write the output
	results <- dcast(res,seqnames+pos+REF~nucleotide+strand,value.var="count",fill=0)
	results$D <- apply(results[,4:ncol(results)],1,sum)
	results$D_forward <- apply(results[,grep('[+]',colnames(results))],1,sum)
	results$D_reverse <- apply(results[,grep('[-]',colnames(results))],1,sum)

	results <- results[,c(1,2,3,grep('D',colnames(results)),grep('[+]',colnames(results)),grep('[-]',colnames(results)))]
	colnames(results) <- sub('[+]','forward',colnames(results))
	colnames(results) <- sub('[-]','reverse',colnames(results))
	colnames(results)[1:2] <- c('CHR','POS')
	
	# temporary file generation
	# outfile <- sub('[.]bam','.out',i)
	outfile <- sub('[.]bam','.out',sub('.*/','',i))	
	write.table(x = results, file = outfile, quote = F, row.names = F, sep = '\t')
}

print('Total time taken...')
time <- proc.time() 
print(paste(time[[1]],'secs',sep=' '))
timestamp()

# R wrapper script for mpileup -> nucleotide frequencies
# mpileup_wrapper_v1.0.R
# version: 1.0
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

timestamp()
suppressMessages(library(Rsamtools, quietly = TRUE))
suppressMessages(library(reshape2, quietly = TRUE))

# command line arguments
args <- commandArgs(trailingOnly = TRUE)
input <- toString(args[1])
positions_file <- toString(args[2])
fasta_file <- toString(args[3])
max_depth <- 1000
min_base_quality <- 0
min_mapq <- 0
include_insertions <- FALSE
include_deletions <- FALSE

if(length(args) >= 4){
	max_depth <- as.integer(args[4])
	if(length(args) >= 5){
		min_base_quality <- as.integer(args[5])
		if(length(args) >= 6){
			min_mapq <- as.integer(args[6])
			if(length(args) >= 7){
				include_deletions <- as.logical(args[7])
				if(length(args) >= 8){
					include_insertions <- as.logical(args[8])
				}
			}
		}
	}
}

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

# get reference base
if(nrow(positions_file) == 1){
    positions_file <- data.frame(V1 = positions_file$V1,
                                 V2 = seq(from = positions_file$V2,
                                          to = positions_file$V3))
    refbase <- getSeq(fasta_file,
                      GRanges(positions_file$V1,
                              IRanges(start = as.numeric(positions_file$V2),
                                      end = as.numeric(positions_file$V2))))
    refbase <- as.data.frame(refbase)$x
    positions_file$REF <- refbase
} else if(nrow(positions_file) > 1){
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
                           min_mapq = min_mapq)

	# call pileup function
	res <- pileup(bf, scanBamParam = param, pileupParam = p_param)

    # if there is no pileup for this file, go to next
    if(nrow(res) =  = 0) next

	# get reference base
	res <- merge(res, positions_file, by.x = c('seqnames', 'pos'), by.y = c('V1', 'V2'))

	# process and write the output
	results <- dcast(res, seqnames+pos+REF~nucleotide+strand,
                     value.var = "count", fill = 0)
	results$D <- apply(results[, 4:ncol(results)], 1, sum)
	results$D_forward <- apply(results[, grep('[+]', colnames(results))],
                               1, sum)
	results$D_reverse <- apply(results[, grep('[-]', colnames(results))],
                               1, sum)

	results <- results[, c(1, 2, 3,
                           grep('D', colnames(results)),
                           grep('[+]', colnames(results)),
                           grep('[-]', colnames(results)))]
	colnames(results) <- sub('[+]', 'forward', colnames(results))
	colnames(results) <- sub('[-]', 'reverse', colnames(results))
	colnames(results)[1:2] <- c('CHR', 'POS')

	# temporary file generation
	# outfile <- sub('[.]bam', '.out', i)
	outfile <- sub('[.]bam', '.out', sub('.*/', '', i))
	write.table(x = results, file = outfile, quote = F,
                row.names = F, sep = '\t')
}

print('Total time taken...')
time <- proc.time()
print(paste(time[[1]], 'secs', sep = ' '))
timestamp()

# gene counts output using STAR
STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which
correspond to different strandedness options:
column 1: gene ID
column 2: counts for unstranded RNA-seq
column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s
yes)
column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s
reverse)

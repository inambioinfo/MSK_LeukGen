# For Spiked-in controls
# DESeq/Limma+Voom Normalizing Rna-Seq Data Using Ercc Spike-In

################################### Using DESeq ###################################
library(DESeq)
# get the size factors from the spiked-ins
# dat.spike.in is count data of only spike ins
cds <- newCountDataSet(countData = dat.spike.in, conditions = as.character(condition))
cds <- estimateSizeFactors(cds)
nf <- sizeFactors(cds) 

# add the size factors to data without spike-ins
# dat.without.spike.in is count data where spike-ins are removed
cds <- newCountDataSet(countData = dat.without.spike.in, conditions = as.character(condition))
sizeFactors(cds) <- nf 

# differential expression test
cds <- estimateDispersions( cds, method= 'per-condition', sharingMode='maximum', fitType='local') 
deseq.res <- nbinomTest( cds, "A", "B")

################################### Using Limma+Voom ###################################
library(limma)

groups = as.factor(condition)
design <- model.matrix(~0+groups)
colnames(design) = levels(groups)
# nf <- calcNormFactors(dat.spike.in) # calculation norm factors using spike ins only
# voom.norm <- voom(dat.without.spike.in, design, lib.size = colSums(countsBuffer) * nf) # add that to voom
# alternative solution
N <- colSums(dat.without.spike.in)
nf <- calcNormFactors(dat.spike.in, lib.size = N)
voom.norm <- voom(dat.without.spike.in, design, lib.size = N * nf)


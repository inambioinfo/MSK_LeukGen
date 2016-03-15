# Run R-GSEA
source('/Users/rathik/Downloads/GSEA-P-R/GSEA.1.0.R', verbose = T)
source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/make_cls.R')
source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/make_gct.R')

rgsea <- function(input_dir, normCounts, annotation, sample.info, res.limma, groupA, groupB) {
  # output directory for R-GSEA analysis
  output <- paste(input_dir,'gsea-out/',sep='')
  
  # make output directory
  system(paste('mkdir',output))
  
  # files
  gctFile <- paste(output,'sample.gct',sep='')
  clsFile <- paste(output,'sample.cls',sep='')
  
  # make gct and cls files
  gct <- make.gct(normCounts = normCounts, annotation = annotation, sample.info = sample.info, test = 'both', limmaout = res.limma)
  cls <- make.cls(sample.conditions, groupA, groupB)
  write.table(x = gct, file = gctFile, quote = F, row.names = F, col.names = F, sep = '\t')
  write.table(x = cls, file = clsFile, quote = F, row.names = F, col.names = F)
  
  # Run GSEA
  GSEA(input.ds = gctFile,
       input.cls = clsFile,
       gs.db = '~/Downloads/GSEA-P-R/GeneSetDatabases/c5.all.v5.0.symbols.gmt',
       output.directory = output)
}

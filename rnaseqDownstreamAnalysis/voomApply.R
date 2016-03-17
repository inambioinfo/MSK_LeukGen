suppressMessages(library(edgeR, quietly = TRUE))

voomApply <- function(counts, sample.info){
  countsBuffer = counts[,match(colnames(counts),sample.info$source_name)]
  condsBuffer = sample.info$source_condition
  groups = as.factor(condsBuffer)
  designCombined <- model.matrix(~0+groups)
  colnames(designCombined) = levels(groups)
  nf <- calcNormFactors(countsBuffer)
  datCombined <- voom(countsBuffer, designCombined, lib.size=colSums(countsBuffer) * nf, plot=FALSE)
  # datCombined <- voom(countsBuffer, designCombined, plot=FALSE)
  return(datCombined$E)
}

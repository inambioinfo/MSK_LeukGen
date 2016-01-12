library(edgeR)
voomApply <- function(counts,groupA,groupB,conditions,sample.info){
  countsBuffer = counts[,match(colnames(counts),sample.info$Source_name)]
  condsBuffer = sample.info$Source_type
  groups = as.factor(condsBuffer)
  designCombined <- model.matrix(~0+groups)
  colnames(designCombined) = levels(groups)
  nf <- calcNormFactors(countsBuffer)
  print(nf)
  datCombined <- voom(countsBuffer, designCombined, plot=T, lib.size=colSums(countsBuffer) * nf)
  # datCombined <- voom(countsBuffer, designCombined, plot=FALSE)
  return(datCombined$E)
}

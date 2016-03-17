limmaApply <- function(normCounts,groupA,groupB,sample.info,annotation){
  countsBuffer = normCounts[,match(colnames(normCounts),sample.info$source_name)]
  condsBuffer = sample.info$source_condition
  groups = as.factor(condsBuffer)
  designCombined <- model.matrix(~0+groups)
  colnames(designCombined) = levels(groups)

  # differential expression using limma
  fit <- lmFit(countsBuffer, design =  designCombined)

  # change groups to fit the model
  groupA <- paste('(',groupA,')',sep = '')
  groupB <- paste('(',groupB,')',sep = '')

  cont.matrix <- makeContrasts(grp = noquote(paste(groupA,"-",groupB,sep='')), levels = designCombined)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, coef = 1, number = Inf, p.value = 1, adjust.method = "BH")
  res$FC <- 2^res$logFC
  res$FC <- ifelse(res$FC<1,-1*(1/res$FC),res$FC)
  res <- res[,c('AveExpr','FC','logFC','P.Value','adj.P.Val')]
  res <- merge(annotation,res,by.x='ID',by.y='row.names')
  return(res)
}

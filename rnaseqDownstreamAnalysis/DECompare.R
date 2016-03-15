require(DESeq)
require(data.table)
DECompare <- function(counts,groupA,groupB,sample.info,annotation){
#   countsBuffer = counts[,match(colnames(counts),sample.info$source_name)]
#   condsBuffer = sample.info$source_type
#   cds <- newCountDataSet( countsBuffer, condsBuffer )
#   cds <- estimateSizeFactors( cds )
#   cds <- estimateDispersions( cds, method= 'pooled', sharingMode='maximum', fitType='local')
#   # normalizedCounts <- t(t(counts(cds)) / sizeFactors(cds))
#   deseq.res <- nbinomTest( cds, groupA, groupB)
#   setnames(deseq.res,old = c("id","foldChange","log2FoldChange","pval","padj"),new = c("ID","FC","logFC","P.Value","adj.P.Val"))
#   deseq.res <- merge(annotation,deseq.res,by='ID')
  countsBuffer = counts[,match(colnames(counts),sample.info$source_name)]
  condsBuffer = sample.info$source_condition
  
  ga <- unlist(strsplit(groupA,""))
  ga <- ga[!ga %in% "+"]
  gb <- unlist(strsplit(groupB,""))
  gb <- gb[!gb %in% "+"]
  
  condsBuffer[which(condsBuffer %in% ga)] <- "groupA"
  condsBuffer[which(condsBuffer %in% gb)] <- "groupB"
  cds <- newCountDataSet(countsBuffer, condsBuffer)
  cds <- estimateSizeFactors( cds )
  cds <- estimateDispersions( cds, method= 'pooled', sharingMode='maximum', fitType='local')
  
  # normalizedCounts <- t(t(counts(cds)) / sizeFactors(cds))
  deseq.res <- nbinomTest( cds, "groupA", "groupB")
  setnames(deseq.res,old = c("id","foldChange","log2FoldChange","pval","padj"),new = c("ID","FC","logFC","P.Value","adj.P.Val"))
  deseq.res <- merge(annotation,deseq.res,by='ID')
  return(deseq.res)
}
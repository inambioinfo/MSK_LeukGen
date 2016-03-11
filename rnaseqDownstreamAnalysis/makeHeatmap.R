library(RColorBrewer)
library(gplots)
library(ggplot2)
library(pheatmap)
#library(GMD)
source('/Users/rathik/scripts/ggplot2_colors.R')

makeHeatmap <- function(results,normCounts,sample.info,title,n){
  
  # sort
  results <- results[order(results$adj.P.Val),]
  
  # only get specified number of top significant degs
  results <- results[1:n,]
  
  # merge limma results with expression matrix
  results <- merge(results,normCounts,by='ID',by.y='row.names')
  
  # annotation row
  annotation_row <- data.frame(Chr=results$Chr,Biotype=results$Biotype) 
  rownames(annotation_row) <- results$Symbol
  
  results <- results[,colnames(results) %in% sample.info$source_name | colnames(results)=='Symbol']

  # make gene symbols as row names
  rownames(results) <- results$Symbol
  results <- results[,-1]
  
  # annotation column
  annotation_col <- data.frame(source_type=sample.info$source_type)
  rownames(annotation_col) <- sample.info$source_name
  
  title <- paste(title,'n = ',n,sep='')
  p <- pheatmap(results, display_numbers = F, scale = 'row',clustering_distance_rows = 'correlation',clustering_distance_cols='correlation',
              main = paste(title), annotation_col = annotation_col, annotation_row = annotation_row,
              cellwidth = 70 , cellheight = 15, fontsize = 14,fontsize_number = 10)
  
  return(p)
  
  }


library(NMF)

res <- limmaApply(normCounts,'Tumor','Normal',condition,sample.info)

n = 50

res <- merge(res,annotation,by.x='row.names',by.y='ID')
res <- merge(res,normCounts,by='Row.names',by.y='row.names')
res <- res[,colnames(res) %in% sample.info$Source_name | colnames(res)=='Symbol']
rownames(res) <- res$Symbol
res <- res[,-1]
res <- res[1:n,]

nr = nrow(res)
nc = ncol(res)
aheatmap(x = as.matrix(res), scale = 'row', distfun = dist2, 
         color = "YlGnBu", 
         main = 'Heatmap of ... Genes')

# distance function
annotation_col <- data.frame(Source_type=sample.info$Source_type)
rownames(annotation_col) <- sample.info$Source_name
p<-pheatmap(res, display_numbers = T, scale = 'row',clustering_distance_rows = 'correlation',
         main = "Differentially Expressed Genes\n", annotation_col = annotation_col,
         cellwidth = 70 , cellheight = 15,srtCol=30)


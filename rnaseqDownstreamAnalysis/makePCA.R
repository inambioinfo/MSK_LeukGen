# make PCA plot before and after normalization
library(ggplot2)
source('/Users/rathik/scripts/theme.R')

makePCA <- function(normCounts, sample.info){

  # top 10% of most variable genes
  n <- nrow(normCounts)*10/100
  
  var <- apply(normCounts,1,var)
  var <- sort(var,decreasing = T)
  var <- data.frame(var = var[1:n])
  var <- dat[rownames(normCounts) %in% rownames(var),]
  pca <- prcomp(t(var))
  scores <- data.frame(Source_type = sample.info$Source_type, Source_name = row.names(pca$x), pca$x)
  tmp <- summary(pca)
  
  pc1 <- round(tmp$importance[2,1]*100,2)
  pc2 <- round(tmp$importance[2,2]*100,2)

	# jpeg(file = paste(title,'.jpg',sep = ''),800,800)
# 	p <- ggplot(data = scores, aes(x = PC1, y = PC2)) + theme.obj + ggtitle(paste("\nTop", '10% most variable', "genes\n", sep = ' ')) +  
# 	  geom_text(aes(label = Source_name, color = Source_type), cex = 7,fontface='bold') 
  p <- ggplot(data = scores, aes(x = PC1, y = PC2)) + theme.obj + ggtitle(paste("\nPCA: Top", '10% most variable', "genes\n", sep = ' ')) +  
    geom_point(aes(shape = Source_type, color = Source_name),cex=7)
	return(p)
	# dev.off()
}


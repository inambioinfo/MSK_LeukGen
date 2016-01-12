# correlation matrix for any data

library(ggplot2)
library(reshape2)

makeCorrPlot <- function(normCounts){
  
  # make correlation matrix
  cormat <- round(cor(normCounts),2)
  
  # reorder
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  # cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  melted_cormat <- melt(upper_tri)
  melted_cormat <- na.omit(melted_cormat)
  
  corrplot <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(melted_cormat$value), 
                         limit = c(min(melted_cormat$value),max(melted_cormat$value)), name="Pearson\nCorrelation") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1),
                            axis.text.y = element_text(size=14),
                            legend.title = element_text(size=14),
                            legend.text = element_text(size=14),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank()) + coord_fixed() + 
    xlab('\nSamples') + ylab('Samples\n') + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) 

  return(corrplot)
}

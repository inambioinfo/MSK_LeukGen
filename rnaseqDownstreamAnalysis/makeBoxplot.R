# make boxplot before and after normalization
library(ggplot2)
library(reshape2)

source('/Users/rathik/scripts/theme.R')
		
makeBoxplot <- function(normCounts, sample.info){
  normCounts.m <- melt(normCounts,varnames = c('ID','Source_name'))
  normCounts.m <- merge(normCounts.m,sample.info,by='Source_name')
	grps <- as.character(sample.info[match(colnames(normCounts), sample.info$Source_name),'Source_type'])
	p <- ggplot(data = normCounts.m, aes(x = Source_name, y = value))+ 
	  geom_boxplot(aes(fill = Source_type)) + theme.obj + ggtitle("\nBoxplot: Voom Normalized Intensities\n") + 
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# 	p <- ggplot(data = normCounts.m, aes(x = Source_name, y = value, color=Source_type,fill=Source_type)) + geom_point(position=position_jitterdodge(dodge.width=0.9)) +
# 	  geom_boxplot(fill="white",outlier.colour = NA, 
# 	               position = position_dodge(width=0.9))
	return(p)
}




# make boxplot before and after normalization
library(ggplot2)
library(reshape2)

source('/Users/rathik/projects/MSK_LeukGen/rnaseqDownstreamAnalysis/theme.R')
		
makeBoxplot <- function(normCounts, sample.info){
  normCounts.m <- melt(normCounts,varnames = c('ID','source_name'))
  normCounts.m <- merge(normCounts.m,sample.info,by='source_name')
	grps <- as.character(sample.info[match(colnames(normCounts), sample.info$source_name),'source_condition'])
	p <- ggplot(data = normCounts.m, aes(x = source_name, y = value))+ 
	  geom_boxplot(aes(fill = source_condition)) + theme.obj + ggtitle("\nBoxplot: Voom Normalized Intensities\n") + 
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# 	p <- ggplot(data = normCounts.m, aes(x = Source_name, y = value, color=Source_type,fill=Source_type)) + geom_point(position=position_jitterdodge(dodge.width=0.9)) +
# 	  geom_boxplot(fill="white",outlier.colour = NA, 
# 	               position = position_dodge(width=0.9))
	return(p)
}




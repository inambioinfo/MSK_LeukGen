theme.set <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, colour = "black"),
                   axis.text.y = element_text(size = 14, colour = "black"), 
                   axis.title.x = element_blank(), axis.title.y = element_text(size=14),
                   title=element_text(size=14,face="bold"))


theme.obj <- theme(axis.text = element_text(vjust = 0.5, colour = "black", size = 14),
                   axis.title = element_text(face = "bold", colour = "black", size = 14),
                   legend.text = element_text(size = 14),
                   legend.title = element_text(size = 14),
                   plot.title = element_text(face = 'bold', size = 16),
                   strip.text = element_text(face = 'bold', size = 12))
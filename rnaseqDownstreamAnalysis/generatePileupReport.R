######## This script reads in multiple mpileup outputs and creates a summary table and plots #######

# process_pileup.R
# Komal S Rathi
# Memorial Sloan Kettering Cancer Center
# Created ##------ Fri Nov  6 14:27:25 2015 ------ #
# Last Modified ##------ Wed Nov 11 16:48:00 2015 ------##
# Note: I have used tidyr extensively in this script. If you are unfamiliar with the package, please refer the package manual.

########################## Usage: Rscript process_pileup.R <project_folder> <Type: Custom or VCF> ##########################

suppressMessages(require(ggplot2,quietly = TRUE))
suppressMessages(require(tidyr,quietly = TRUE))
suppressMessages(require(dplyr,quietly = TRUE))
suppressMessages(require(scales,quietly = TRUE))
suppressMessages(require(tools,quietly = TRUE))
suppressMessages(require(gridExtra,quietly = TRUE))
suppressMessages(require(limma,quietly = TRUE))
suppressMessages(require(reshape2,quietly = TRUE))
suppressMessages(require(stringr,quietly = TRUE))
suppressMessages(require(VennDiagram,quietly = TRUE))

timestamp()

# import dependencies
# imports ggplot2 themes & default color scheme
source('/Users/rathik/scripts/theme.R') 
source('/Users/rathik/scripts/ggplot2_colors.R') 

# transform from long format to wide format
transform.long2wide <- function(tt)
{
  n <- names(tt)
  rr <- reshape(tt, direction = 'wide', sep = '~', idvar = n[c(1:2,4:5)], v.names = n[6:9], timevar = n[3])
  idx <- grepl('~', names(rr))
  names(rr)[idx] <- sapply(strsplit(names(rr)[idx], '~'), function(x) paste0(rev(x), collapse = '_'))
  return(rr)
}

# read commandline arguments i.e. project directory
args <- commandArgs(trailingOnly = TRUE)
input_dir <- toString(args[1])
datatype <- toString(args[2])

# get all the files with .out extension
files <- list.files(path = input_dir, pattern = "*_pileup.out", full.names = T, recursive = TRUE)
pileup <- do.call("rbind", lapply(files, function(fn) data.frame(SAMPLE=fn, read.table(fn,stringsAsFactors=FALSE,header=T))))

# process pileup to long format
pileup.long <- pileup %>%
  gather(NRA, value, -SAMPLE, -CHR, -POS, -REF, -D, -D_forward, -D_reverse) %>% # melt except SAMPLE...D_reverse
  mutate(NRA = sub('_.*','',NRA), SAMPLE = gsub('.*/|_.*','',SAMPLE)) %>% # add column for alternate allele
  subset(value != 0) %>% # remove 0 counts
  group_by(SAMPLE, CHR, POS, NRA) %>%
  mutate(DIR = n(), DNRA = sum(value)) %>% # get Depth of NRA, and directionality
  mutate(prop = value/DNRA) %>% # proportion of directionality
  group_by(SAMPLE, CHR, POS, NRA) %>%
  mutate(PDIR = min(prop)) %>% # take minimum of the two proportions, in case of bi-directional
  mutate(PDIR = ifelse(DIR == 1, 0, PDIR)) %>% # unidirectional proportions are 0
  select(CHR, POS, SAMPLE, REF, NRA, D, DNRA, DIR, PDIR) %>% # select only these columns
  filter(NRA!=REF) %>% # remove ref==nra
  unique %>% # get unique rows
  as.data.frame() 

# tranform long format to wide format
pileup.wide <- transform.long2wide(pileup.long)
output <- paste(input_dir,'mpileup_summary_report.txt',sep='/')
write.table(x = pileup.wide,file = output, quote=F, row.names = F, sep='\t')

# modify data -> pileup.plot -> use for making plots
pileup.plot <- pileup.long %>% 
  group_by(CHR,POS,REF,NRA) %>%
  mutate(SAMPLES = paste(SAMPLE,collapse=','),SAMPLE.COUNT=n(),GROUPS=paste(DIR,collapse=',')) %>%
  mutate(SAMPLE.COUNT = ifelse(SAMPLE.COUNT == 1, SAMPLES, paste('Shared.By.',SAMPLE.COUNT,'T',sep=''))) %>%
  mutate(SAMPLE.NORMAL = ifelse(grepl('N1',SAMPLES),'Shared By Normal','Not Shared')) %>%
  mutate(SAMPLE.COUNT = ifelse(SAMPLE.NORMAL=='Shared By Normal' & GROUPS!=1,'Shared in Normal',SAMPLE.COUNT)) %>%
  mutate(GROUPS = lapply(GROUPS, function(x) paste(sort(unlist(strsplit(x, ","))), collapse = ""))) %>% # make strings of 1s and 2s
  mutate(Unidirectional = as.numeric(lapply(GROUPS, function(x) str_count(x,'1'))), # count the number of 1s
         Bidirectional = as.numeric(lapply(GROUPS, function(x) str_count(x,'2')))) %>% # count the number of 2s
  mutate(PERC = DNRA/D*100) %>% # variant allele frequency
  as.data.frame()

# define groups based on how many samples have uni/bidirectional variants
pileup.plot$Direction <- paste(paste(pileup.plot$Unidirectional,'Unidirectional'), paste(pileup.plot$Bidirectional,'Bidirectional'), sep='; ')
pileup.plot$Direction[pileup.plot$Unidirectional == 0] <- 'All Bidirectional'
pileup.plot$Direction[pileup.plot$Bidirectional == 0] <- 'All Unidirectional'

# plot1: summary barplot
print('Creating Barplot...')
bar.plot <- ggplot(pileup.plot, aes(factor(SAMPLE.COUNT),fill=Direction)) + geom_bar(color='black',size=0.5) + theme_bw() +
  theme.obj + xlab('') + ylab('Count\n') + ggtitle('VARIANT SUMMARY PLOT\n') + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot2: VAF density plot
print('Creating Density plot...')
density.plot <- ggplot(pileup.plot, aes(x=PERC,fill=factor(DIR))) + geom_density(alpha=0.7,size=0.5) + 
  ylab('Density\n') + xlab('\n%VAF') + theme_bw() + theme.obj + ggtitle('VARIANT ALLELE FREQUENCY\n') + 
  scale_fill_discrete(name="Direction", breaks=c("1", "2"), labels=c("Unidirectional", "Bidirectional")) + 
  facet_wrap(~SAMPLE.COUNT,nrow = 2)

# draw venn diagram
print('Creating Venn Diagram...')
dd <- dcast(pileup.plot,CHR+POS+REF+NRA~SAMPLE,length,value.var='PERC')
dd <- dd[,grep('T',colnames(dd))]
col <- gg_color_hue(3)
p3 <- draw.triple.venn(area1 = nrow(dd[which(dd[[1]]==1),]),
                       area2 = nrow(dd[which(dd[[2]]==1),]),
                       area3 = nrow(dd[which(dd[[3]]==1),]),
                       n12 = nrow(dd[which(dd[[1]]==1 & dd[[2]]==1),]),
                       n23 = nrow(dd[which(dd[[2]]==1 & dd[[3]]==1),]),
                       n13 = nrow(dd[which(dd[[1]]==1 & dd[[3]]==1),]),
                       n123 = nrow(dd[which(dd[[1]]==1 & dd[[2]]==1 & dd[[3]]==1),]),
                       category = colnames(dd), fill = col, cat.cex = 1, cex = 1, cat.dist = c(0.01,0.01,0.1),main='Variants')

# text to be printed on the report
txt <- paste(paste('Number of Samples:',length(files)), 
             paste('Normal1:',unique(pileup.long$SAMPLE)[1],
                   '\nTumor1:',unique(pileup.long$SAMPLE)[2],
                   '\nTumor2:',unique(pileup.long$SAMPLE)[3],
                   '\nTumor3:',unique(pileup.long$SAMPLE)[4]),
             paste('Number of Variants:',nrow(unique(pileup[,c('CHR','POS')]))), 
             paste('Type:',datatype),
             sep='\n')

# print plots to PDF
print('Writing output files...')
output <- paste(input_dir,'mpileup_summary_report.pdf',sep='/')
pdf(output, width = 20, height = 20)
grid.arrange(bar.plot,textGrob(txt, just = 'center', x = 0.5, y=0.5, gp = gpar(font=2)),gTree(children = p3),density.plot,
             layout_matrix = rbind(c(1,1,2), 
                                   c(1,1,3), 
                                   c(4,4,4),
                                   c(4,4,4)))
invisible(dev.off())

print('Total time taken...')
time <- proc.time() 
print(paste(time[[1]],'secs',sep=' '))
timestamp()
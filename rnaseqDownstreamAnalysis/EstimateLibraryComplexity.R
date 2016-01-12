
EstimateLibraryComplexity <- function(input_dir){

    # read files in and add Sample names
    getFile <- function(x){
    d <- read.delim(x, skip=10)
    colnames(d)[2] <- 'Reads'
    d$Sample <- sub('_rnaseqFastqInput_*','',sub('.*/','',x))
    return(d)
  }
  
  # get all files with libcomp.txt extension
  d.list <- list.files(path = input_dir, pattern = "*libcomp.txt", full.names = T)  
  d1 <- do.call('rbind',lapply(d.list,getFile))
  
  # density plot
  p <- ggplot(data = d1, aes(x = log2(Reads),fill = Sample)) + geom_density(alpha = 0.4) + ylab('Duplication Density') + theme.obj
  
  return(p)  
}

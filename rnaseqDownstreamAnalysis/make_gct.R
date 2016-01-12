make.gct <- function(normCounts,annotation,sample.info,limmaout,test){

  # subset the dataframe based on what set you want to test
  if(test=='both'){
    normCounts <- normCounts[rownames(normCounts) %in% limmaout[limmaout$adj.P.Val<0.05,'ID'],]  
  } else if(test=='up'){
    normCounts <- normCounts[rownames(normCounts) %in% limmaout[limmaout$adj.P.Val<0.05 & limmaout$FC > 0,'ID'],]
  } else if(test=='down'){
    normCounts <- normCounts[rownames(normCounts) %in% limmaout[limmaout$adj.P.Val<0.05 & limmaout$FC < 0,'ID'],]
  }
  
  # get the expression values & symbols
  nr <- nrow(normCounts)
  nsample <- ncol(normCounts)
  annotation$DESCRIPTION <- 'na'
  gct <- merge(annotation, normCounts, by.x='ID', by.y='row.names')
  gct <- gct[,colnames(gct) %in% sample.info$Source_name | colnames(gct)=='DESCRIPTION' | colnames(gct)=='Symbol']
  gct$Symbol <- toupper(gct$Symbol)
  colnames(gct)[1] <- c('NAME')
  
  # make the initial part before expression values 
  dt <- matrix('',ncol=2+nsample,nrow = 3)
  dt[1,1] <-  '#1.2'
  dt[2,1:2] <- c(nr,nsample)
  dt[3,] <- colnames(gct)
  dt <- data.frame(dt,stringsAsFactors=F)
  colnames(dt) <- colnames(gct)
  
  # rbind the two parts
  gct <- rbind(dt,gct)
  return(gct)
}
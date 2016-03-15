make.cls <- function(sample.info, groupA, groupB){
  
  condsBuffer = sample.info$source_condition
  ga <- unlist(strsplit(groupA,""))
  ga <- ga[!ga %in% "+"]
  gb <- unlist(strsplit(groupB,""))
  gb <- gb[!gb %in% "+"]
  
  condsBuffer[which(condsBuffer %in% ga)] <- "groupA"
  condsBuffer[which(condsBuffer %in% gb)] <- "groupB"
  
  first <- paste(nrow(sample.info),nlevels(factor(condsBuffer)),1)
  second <- paste('#',paste(as.character(levels(factor(condsBuffer))),collapse = ' '))
  third <- paste(as.character(condsBuffer),collapse = ' ')
  cls <- data.frame(V1=c(first,second,third))

  return(cls)
}
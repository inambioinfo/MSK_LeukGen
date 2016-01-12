make.cls <- function(sample.info){
  first <- paste(nrow(sample.info),nlevels(sample.info$Source_type),1)
  second <- paste('#',paste(as.character(levels(sample.info$Source_type)),collapse = ' '))
  third <- paste(as.character(sample.info$Source_type),collapse = ' ')
  cls <- data.frame(V1=c(first,second,third))
  return(cls)
}
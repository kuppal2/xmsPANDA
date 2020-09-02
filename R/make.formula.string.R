make.formula.string <-
function(factors, do.interactions=FALSE){
  fs = "1"
  if(length(factors)){
    fs = paste(factors, collapse=" + ")
    if(do.interactions && length(factors) > 1)
      fs = paste(unlist(lapply(as.data.frame(t(combinations(length(factors), 2, factors)), stringsAsFactors=F), paste, collapse="*")), collapse = " + ")
  }
  return(fs)
}

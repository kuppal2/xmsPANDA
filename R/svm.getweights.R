svm.getweights <-
function(model){
  w=0
  if(model$nclasses==2){
    w=t(model$coefs)%*%model$SV
  }else{
    
    #when we deal with OVO svm classification
    ## compute start-index
    start <- c(1, cumsum(model$nSV)+1)
    start <- start[-length(start)]
    calcw <- function (i,j) {
      ## ranges for class i and j:
      ri <- start[i] : (start[i] + model$nSV[i] -1)
      rj <- start[j] : (start[j] + model$nSV[j] -1)
      ## coefs for (i,j):
      coef1 <- model$coefs[ri, j-1]
      coef2 <- model$coefs[rj, i]
      ## return w values:
      w=t(coef1)%*%model$SV[ri,]+t(coef2)%*%model$SV[rj,]
      return(w)
    }
    W=NULL
    for (i in 1 : (model$nclasses - 1)){
      for (j in (i + 1) : model$nclasses){
        wi=calcw(i,j) 
        W=rbind(W,wi) 
      } 
    } 
    w=W 
  } 
  return(w) 
}

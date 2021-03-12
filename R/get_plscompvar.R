get_plscompvar <-
function(p2,nvar,h) {
  
  
  b<-c(p2$loadings$Y)[1:h]
  T<-p2$variates$X[,1:h]
  SS<-b^2 * colSums(T^2)
  
  
  return(SS)
}

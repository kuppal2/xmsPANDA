get_plscompvar <-
function(p2,nvar,h) {
<<<<<<< HEAD
  
  
  b<-c(p2$loadings$Y)[1:h]
  T<-p2$variates$X[,1:h]
  SS<-b^2 * colSums(T^2)
  
  
  return(SS)
=======
    
    
    b<-c(p2$loadings$Y)[1:h]
    T<-p2$variates$X[,1:h]
    SS<-b^2 * colSums(T^2)
    
    
    return(SS)
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

docortest <-
function(n,r){
<<<<<<< HEAD
  #n<-length(x)
  t<-r*sqrt((n-2)/(1-r^2))
  pvalue<-2*pt(-abs(t),df=n-2)
  
  return(pvalue)
  
=======
	#n<-length(x)
    t<-r*sqrt((n-2)/(1-r^2))
    pvalue<-2*pt(-abs(t),df=n-2)
	
	return(pvalue)
	
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

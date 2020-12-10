do_rsd <-
function(x){
<<<<<<< HEAD
  
  sd_val<-sd(x,na.rm=TRUE)
  mean_val<-mean(x,na.rm=TRUE)
  cv_val<-100*(sd_val/(mean_val))
  
  return(cv_val)
=======
	
	sd_val<-sd(x,na.rm=TRUE)
	mean_val<-mean(x,na.rm=TRUE)
	cv_val<-100*(sd_val/(mean_val))
	
	return(cv_val)
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

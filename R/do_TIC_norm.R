do_TIC_norm <-
function(data_m){
<<<<<<< HEAD
  
  
  sum_int<-apply(data_m,2,function(x){sum(x,na.rm=TRUE)})
  
  median_int<-apply(data_m,2,function(x){median(x,na.rm=TRUE)})
  
  data_m<-sweep(data_m,2,(sum_int/median_int),'/')
  return(data_m)
=======
      
     
      sum_int<-apply(data_m,2,function(x){sum(x,na.rm=TRUE)})
      
      median_int<-apply(data_m,2,function(x){median(x,na.rm=TRUE)})
      
      data_m<-sweep(data_m,2,(sum_int/median_int),'/')
      return(data_m)
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

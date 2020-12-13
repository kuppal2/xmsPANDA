do_TIC_norm <-
function(data_m){
  
  
  sum_int<-apply(data_m,2,function(x){sum(x,na.rm=TRUE)})
  
  median_int<-apply(data_m,2,function(x){median(x,na.rm=TRUE)})
  
  data_m<-sweep(data_m,2,(sum_int/median_int),'/')
  return(data_m)
}

do_log2transform_norm <-
function(data_m,log2.transform.constant=1){
  
  data_m<-log2(data_m+log2.transform.constant)
  return(data_m)
  
}

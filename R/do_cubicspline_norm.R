do_cubicspline_norm <-
function(data_m){
  
  data_m<-normalize.qspline(data_m)
  return(data_m)
}

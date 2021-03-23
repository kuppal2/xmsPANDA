do_rangescaling_norm <-
function(data_m)
{
  #range scaling
  data_m <- apply(data_m,1,function(x){(x-min(x))/(max(x)-min(x))})
  
  data_m<-t(data_m)
  return(data_m)
}

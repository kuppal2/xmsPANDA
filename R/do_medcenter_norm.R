do_medcenter_norm <-
function(data_m)
{
  colmedians=apply(data_m,1,function(x){median(x,na.rm=TRUE)})
  data_m=sweep(data_m,1,colmedians)
  return(data_m)
  
}

do_medcenter_norm <-
function(data_m)
<<<<<<< HEAD
{
  colmedians=apply(data_m,1,function(x){median(x,na.rm=TRUE)})
  data_m=sweep(data_m,1,colmedians)
  return(data_m)
  
}
=======
  {
      colmedians=apply(data_m,1,function(x){median(x,na.rm=TRUE)})
      data_m=sweep(data_m,1,colmedians)
      return(data_m)
      
  }
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3

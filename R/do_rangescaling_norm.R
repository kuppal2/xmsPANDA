do_rangescaling_norm <-
function(data_m)
<<<<<<< HEAD
{
  #range scaling
  data_m <- apply(data_m,1,function(x){(x-min(x))/(max(x)-min(x))})
  
  data_m<-t(data_m)
  return(data_m)
}
=======
  {
      #range scaling
      data_m <- apply(data_m,1,function(x){(x-min(x))/(max(x)-min(x))})
      
      data_m<-t(data_m)
      return(data_m)
  }
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3

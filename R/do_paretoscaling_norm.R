do_paretoscaling_norm <-
function(data_m)
{
<<<<<<< HEAD
  #pareto scaling
  data_m <- apply(data_m, 1, function(x){(x - mean(x))/sqrt(sd(x))})
  
  data_m<-data_m/mean(abs(data_m))
  data_m<-data_m+1000
  data_m<-t(data_m)
  
  return(data_m)
  
=======
      #pareto scaling
      data_m <- apply(data_m, 1, function(x){(x - mean(x))/sqrt(sd(x))})
      data_m<-t(data_m)

       return(data_m)
      
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

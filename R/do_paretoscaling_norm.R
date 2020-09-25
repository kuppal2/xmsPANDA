do_paretoscaling_norm <-
function(data_m)
{
      #pareto scaling
      data_m <- apply(data_m, 1, function(x){(x - mean(x))/sqrt(sd(x))})
      data_m<-t(data_m)

       return(data_m)
      
}

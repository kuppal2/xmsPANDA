do_znormtransform_norm <-
function(data_m)
{
  #autoscaling; z-scaling
  data_m<-scale(t(data_m))
  data_m<-t(data_m)
  return(data_m)
}

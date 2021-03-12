do_loess_norm <-
function(data_m)
{
  data_m<-normalizeCyclicLoess(data_m)
  return(data_m)
}

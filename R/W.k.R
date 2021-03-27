W.k <-
function(x, clus,d.power=1) {
  n <- nrow(x);ii <- seq_len(n)
  res<-suppressWarnings(0.5 * sum(vapply(split(ii, clus), function(I) {
    xs <- x[I, , drop = FALSE]
    sum(dist(xs)^d.power/nrow(xs))
  }, 0)))
  
 print(res)
  logW.k<-log(res+0.001)
  return(logW.k)
}

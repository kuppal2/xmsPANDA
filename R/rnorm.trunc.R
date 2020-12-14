rnorm.trunc <-
function (n, mu, sigma, lo=-Inf, hi=Inf){
  # Calculates truncated noraml
  p.lo = pnorm (lo, mu, sigma)
  p.hi = pnorm (hi, mu, sigma)
  u = runif (n, p.lo, p.hi)
  return (qnorm (u, mu, sigma))
}

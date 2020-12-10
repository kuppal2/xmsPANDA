rnorm.trunc <-
function (n, mu, sigma, lo=-Inf, hi=Inf){
<<<<<<< HEAD
  # Calculates truncated noraml
=======
# Calculates truncated noraml
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
  p.lo = pnorm (lo, mu, sigma)
  p.hi = pnorm (hi, mu, sigma)
  u = runif (n, p.lo, p.hi)
  return (qnorm (u, mu, sigma))
}

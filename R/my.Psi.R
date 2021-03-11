my.Psi <-
function(x, my.pi){
  # calculates Psi
  exp(log(1-my.pi)  + dnorm(x, 0, 1, log=T) - log(my.pi + (1 - my.pi) * pnorm(x, 0, 1) ))
}

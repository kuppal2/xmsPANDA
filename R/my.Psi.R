my.Psi <-
function(x, my.pi){
<<<<<<< HEAD
  # calculates Psi
  exp(log(1-my.pi)  + dnorm(x, 0, 1, log=T) - log(my.pi + (1 - my.pi) * pnorm(x, 0, 1) ))
=======
# calculates Psi
exp(log(1-my.pi)  + dnorm(x, 0, 1, log=T) - log(my.pi + (1 - my.pi) * pnorm(x, 0, 1) ))
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

my.Psi.dash <-
function(x, my.pi){
  # calculates the derivative of Psi
  -my.Psi(x, my.pi) * (x + my.Psi(x, my.pi))
}

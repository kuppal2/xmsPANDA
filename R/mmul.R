mmul <-
function(A, B){
  # multiply square matrix by rectangle with NA's (???)
  X = A
  Y = B
  X[is.na(A)] = Y[is.na(B)] = 0
  R = X %*% Y
  R[is.na(B)] = NA
  R
}

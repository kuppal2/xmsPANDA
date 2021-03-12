check_model <-
function(model.fit,dependent.var){
  
  
  residuals <- resid(model.fit)
  plot(fitted(model.fit), residuals)
  abline(0,0)
  
  plot(fitted(model.fit), dependent.var)
  
  qqnorm(residuals)
  qqline(residuals)
}

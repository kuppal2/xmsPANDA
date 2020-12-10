check_model <-
function(model.fit,dependent.var){
<<<<<<< HEAD
  
  
  residuals <- resid(model.fit)
  plot(fitted(model.fit), residuals)
  abline(0,0)
  
  plot(fitted(model.fit), dependent.var)
  
  qqnorm(residuals)
  qqline(residuals)
=======
    
    
    residuals <- resid(model.fit)
    plot(fitted(model.fit), residuals)
    abline(0,0)
    
    plot(fitted(model.fit), dependent.var)
    
    qqnorm(residuals)
    qqline(residuals)
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

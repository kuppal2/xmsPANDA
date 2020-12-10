ggadjustedcurves.single <-
function(data, fit) {
<<<<<<< HEAD
  time <- surv <- variable <- NULL
  
  pred <- survexp(~1, data = data, ratetable = fit)
  
  curve <- data.frame(time = c(0,pred$time),
                      variable = "total",
                      surv = c(1, pred$surv))
  
=======
    time <- surv <- variable <- NULL
    
    pred <- survexp(~1, data = data, ratetable = fit)
    
    curve <- data.frame(time = c(0,pred$time),
    variable = "total",
    surv = c(1, pred$surv))
    
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

ggadjustedcurves.single <-
function(data, fit) {
  time <- surv <- variable <- NULL
  
  pred <- survexp(~1, data = data, ratetable = fit)
  
  curve <- data.frame(time = c(0,pred$time),
                      variable = "total",
                      surv = c(1, pred$surv))
  
}

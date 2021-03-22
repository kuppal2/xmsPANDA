ggadjustedcurves.average <-
function(data, fit, variable) {
  time <- surv <- NULL
  
  lev <- sort(unique(data[,variable]))
  pred <- survexp(as.formula(paste("~", variable)), data = data,
                  ratetable = fit)
  
  curve <- data.frame(time = rep(c(0,pred$time), length(lev)),
                      variable = factor(rep(lev, each=1+length(pred$time))),
                      surv = c(rbind(1, pred$surv)))
  
}

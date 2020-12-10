ggadjustedcurves.average <-
function(data, fit, variable) {
<<<<<<< HEAD
  time <- surv <- NULL
  
  lev <- sort(unique(data[,variable]))
  pred <- survexp(as.formula(paste("~", variable)), data = data,
                  ratetable = fit)
  
  curve <- data.frame(time = rep(c(0,pred$time), length(lev)),
                      variable = factor(rep(lev, each=1+length(pred$time))),
                      surv = c(rbind(1, pred$surv)))
  
=======
    time <- surv <- NULL
    
    lev <- sort(unique(data[,variable]))
    pred <- survexp(as.formula(paste("~", variable)), data = data,
    ratetable = fit)
    
    curve <- data.frame(time = rep(c(0,pred$time), length(lev)),
    variable = factor(rep(lev, each=1+length(pred$time))),
    surv = c(rbind(1, pred$surv)))
    
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

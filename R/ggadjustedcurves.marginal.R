ggadjustedcurves.marginal <-
function(data, fit, variable) {
  time <- surv <- NULL
  
  lev <- sort(unique(data[,variable]))
  ndata <- data[rep(1:nrow(data), each=length(lev)),
                setdiff(colnames(data), variable)]
  ndata[,variable] = rep(lev, nrow(data))
  
  pred <- survexp(as.formula(paste("~", variable)), data = ndata,
                  ratetable = fit)
  # remove leading zeros
  # while survexp returns non monotonic results
  if (length(dim(pred$surv)) == 2) {
    for (i in 1:ncol(pred$surv))
      for (j in nrow(pred$surv):2)
        if (pred$surv[j,i] > pred$surv[j - 1,i])
          pred$surv[j - 1,i] <- 1
  }
  
  curve <- data.frame(time = rep(c(0,pred$time), length(lev)),
                      variable = factor(rep(lev, each=1+length(pred$time))),
                      surv = c(rbind(1, pred$surv)))
}

ROCinternal.panda <-
function(test, resp, plot, ...)
{
  dotsCall <- substitute(list(...))
  ll <- eval(dotsCall)
  if(!hasArg(xlab)) ll$xlab <- "Threshold for assignment to class 1"
  if(!hasArg(ylab)) ll$ylab <- "specificity for class 0"
  if(!hasArg(main)) ll$main <- "Receiver Operator Characteristic"
  if(!hasArg(lwd)) ll$lwd <- 2
  m <- as.matrix(table(test, resp))
  fv <- as.numeric(row.names(m))
  nr <- dim(m)[1]
  a <- apply(m, 2, sum)
  m <- addmargins(m, 2)
  m <- apply(m[nr:1, ], 2, cumsum)[nr:1, ]
  sn <- c(m[, 2]/a[2], 0)
  sp <- c((a[1] - m[, 1])/a[1], 1)
  pvp <- c(m[, 2]/m[, 3], 1)
  pvn <- (a[1] - m[, 1])/(sum(a) - m[, 3])
  pvn <- c(pvn, rev(pvn)[1])
  res <- data.frame(cbind(sn, sp, pvp, pvn, c(NA, fv)))
  auc <- sum((res[-1, 1] +res[-nr, 1])/2 * diff(res[, 2]))
  #xl <- range(test)
  #ll$x <- xl
  #ll$y <- 0:1
  #ll$xlim <- xl
  #ll$ylim <- 0:1
  #ll$type <- "n"
  ll$x <- 1-res[,2]
  ll$y <- res[,1]
  ll$xlim <- 0:1
  ll$xlab <- "1-specificity"
  ll$ylim <- 0:1
  ll$ylab <- "Sensitivity"
  ll$type <- "n"
  if(plot){
    do.call("plot", args=ll)
    #plot(xl, 0:1, xlim = xl, xlab = paste(deparse(substitute(test)),
    #    "(grid at deciles)"), ylim = 0:1, ylab = " ",
    #    type = "n")
    
    #plot(1 - res[, 2], res[, 1], xlim = 0:1, xlab = "1-Specificity",
    #       ylim = 0:1, ylab = "Sensitivity", type = "n", ...)
    #   if (grid)
    #       abline(h = 0:10/10, v = 0:10/10, col = gray(0.9))
    #   abline(0, 1, col = gray(0.4))
    box()
    ll$type <- "l"
    do.call("lines", args = ll)
    
    
    #box()
    #ll$xlim <- NULL
    #ll$ylim <- NULL
    #ll$type <- "l"
    #ll$x <- fv
    #ll$y <- res[, 2]
    #ll$left <- TRUE
    #ll$order <- TRUE
    #do.call("steplines", args=ll)
    plot(function(x) x, from=0, to=1, lty="dashed", add=TRUE)
    text(0.8, 0.1, cex=2, label=paste("AUC=", round(auc,3), sep=""))
  }
  names(auc) <- "auc"
  return(list("plotcoordinates"=ll,"auc"=auc))
}

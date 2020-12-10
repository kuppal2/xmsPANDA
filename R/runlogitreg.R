runlogitreg <-
function(X,Y,fdrmethod="BH",fdrthresh=0.05,pvalue.thresh=0.05,robust.estimate=FALSE){
<<<<<<< HEAD
  options(warn=-1)
  logistic_reg=TRUE
  res<-runlmreg(X=X,Y=Y,fdrmethod=fdrmethod,fdrthresh=fdrthresh,pvalue.thresh=pvalue.thresh,logistic_reg=TRUE,robust.estimate=robust.estimate)
  options(warn=0)
  return(res)
=======
     options(warn=-1)
    logistic_reg=TRUE
    res<-runlmreg(X=X,Y=Y,fdrmethod=fdrmethod,fdrthresh=fdrthresh,pvalue.thresh=pvalue.thresh,logistic_reg=TRUE,robust.estimate=robust.estimate)
     options(warn=0)
    return(res)
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

do_ridge_lasso_elasticnet <-
function(X,classlabels,alpha.val=1){ #,lambda=10^seq(10, -2, length = 100)){
  
  #temp1<-cbind(classlabels,t(X))
  #temp1<-as.data.frame(temp1)
  
  x<-t(X)
  x<-as.matrix(x)
  classlabels=as.numeric(as.factor(classlabels))
  cvfit<-cv.glmnet(x,classlabels, alpha = alpha.val) #, lambda = lambda)	
  
  plot(cvfit)
  res<-coef(cvfit, s = "lambda.min")	
  
  res<-res[-c(1),]
  
  #res_sum<-summary(res)
  
  return(res) #_sum)
}

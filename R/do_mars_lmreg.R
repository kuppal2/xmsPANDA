do_mars_lmreg <-
function(dataA){
  
  
  dataA<-as.data.frame(dataA)
  
  
  #print(dim(dataA))
  mars_res<-new("list")
  
  
  mars_res[[1]]<-try(earth(dataA$Response~.,data=as.data.frame(dataA),degree=1,ncross=10,nfold=10),silent=TRUE)
  mars_res[[2]]<-try(earth(dataA$Response~.,data=as.data.frame(dataA),degree=2,ncross=10,nfold=10),silent=TRUE)
  mars_res[[3]]<-try(earth(dataA$Response~.,data=as.data.frame(dataA),degree=3,ncross=10,nfold=10),silent=TRUE)
  
  gcv_list<-c(mars_res[[1]]$gcv,mars_res[[2]]$gcv,mars_res[[3]]$gcv)
  
  min_gcv<-which(gcv_list==min(gcv_list,na.rm=TRUE))
  
  marsfitdeg<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=min_gcv[1],ncross=kfold,nfold=10)
  
  print(summary(marsfitdeg))
  #marsfitdeg<-earth(mz~Strain+Treatment+Batch,data=as.data.frame(curdata),degree=3,ncross=10,nfold=10)
  #marsfitdeg<-try(earth(dataA$Response~.,data=as.data.frame(dataA),degree=1,ncross=10,nfold=10),silent=TRUE)
  
  if (is(marsfitdeg, "try-error")){
    lmtest_pval<-NA
  }else{
    #print(summary(marsfitdeg))
    #evdeg1 <- evimp(marsfitdeg)
    #plot(evdeg1)
    qqnorm(marsfitdeg$residuals)
    qqline(marsfitdeg$residuals)
    shapiro.test(marsfitdeg$residuals)
    bx1 <- model.matrix(marsfitdeg)
    #deg.lm <- glm(as.vector(curdata[,1]) ~ bx1[,-1]) # -1 to drop intercept
    deg.lm <- try(glm(as.numeric(dataA$Response) ~ bx1[,-1]),silent=TRUE)
    if (is(deg.lm, "try-error")){
      lmtest_pval<-1
    }else{
      slm<-summary(deg.lm)
      lmtest_pval=slm$coefficients[,4]
      #print(slm$coefficients)
      
    }
    
    #slm<-summary(deg.lm) # yields same coeffs as above summary
    #aicdeg1[i]<-s1$aic
    #plot(effect('Dept:Gender', berk.mod2), multiline=TRUE)
  }
  return(list("marsummary"=summary(marsfitdeg), "pvalues"=lmtest_pval))
  
}

diffexplmreg <-
function(dataA,logistic_reg=FALSE,poisson_reg=FALSE,robust.estimate=FALSE,vcovHC.type="HC3"){
  
  dataA<-as.data.frame(dataA)
  
  ###savedataA,file="lmreg_func.Rda")
  if(logistic_reg==TRUE){
    cnames1<-colnames(dataA)
    cnames1[2]<-"Class"
    colnames(dataA)<-cnames1
    
    labels_1<-levels(as.factor(dataA$Class))
    
    dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[1]),0)
    dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[2]),1)
    
    
    a1 <- glm(dataA$Class ~ .,family=binomial(logit),data=dataA)
  }else{
    
    if(poisson_reg==TRUE){
      
      cnames1<-colnames(dataA)
      cnames1[2]<-"Class"
      colnames(dataA)<-cnames1
      
      ####savedataA,file="temp1.Rda")
      labels_1<-levels(as.factor(dataA$Class))
      dataA$Class<-as.numeric(dataA$Class)
      
      a1 <- glm(dataA$Class ~ .,family=poisson(log),data=dataA)
      
    }else{
      
      cnames1<-colnames(dataA)
      cnames1[2]<-"Class"
      colnames(dataA)<-cnames1
      dataA$Class<-as.numeric(dataA$Class)
      a1 <- lm(dataA$Response ~ .,data=dataA) # aov(dataA$Response ~ .,data=dataA) # + chocolate$Factor1*chocolate$Factor2)
      
    }
  }
  s1<-summary(a1)
  
  
  if(logistic_reg==FALSE){
    
    r2<-s1$adj.r.squared
  }else{
    r2<-NA
  }
  
  if(robust.estimate==FALSE){
    
    
    s1<-s1$coefficients
  }else{
    cov.a1 <- vcovHC(a1,vcovHC.type) #using default: HC3 #, type="HC0")
    std.err <- sqrt(diag(cov.a1))
    s1 <- cbind(Estimate= coef(a1), "Robust SE" = std.err, "z value"=coef(a1)/std.err,
                "Pr(>|z|)" = 2 * pnorm(abs(coef(a1)/std.err), lower.tail=FALSE),
                LL = coef(a1) - 1.96 * std.err,
                UL = coef(a1) + 1.96 * std.err)
    
    
  }
  
  #save(s1,file="s1.Rda")
  
  if(nrow(s1)>1){
    
    s1<-s1[-c(1),]
    
    if(dim(dataA)[2]<3){ # && dim(dataA)[1]<3){
      #s1<-as.data.frame(s1)
      s1<-t(s1)
      
    }
    
    
    confint_lower<-s1[,1]-(1.96*s1[,2])
    confint_upper<-s1[,1]+(1.96*s1[,2])
    
    
    return(list("mainpvalues"=s1[,4],"estimates"=s1[,1],"statistic"=s1[,3],"stderr"=s1[,2],"r2"=r2,"confint"=c(confint_lower,confint_upper)))
  }else{
    
    return(list("mainpvalues"=NA,"estimates"=NA,"statistic"=NA,"stderr"=NA,"r2"=r2,"confint"=c(NA,NA)))
  }
  
}

diffexponewayanova <-
function(dataA,covar.matrix=NA){
  
  dataA<-as.data.frame(dataA)
  
  ###savedataA,file="dataA.Rda")
  
  a1 <- aov(dataA$Response ~ .,data=dataA) # + chocolate$Factor1*chocolate$Factor2)
  
  posthoc <- TukeyHSD(x=a1, conf.level=0.95,test=univariate())
  anova_res<-anova(a1)
  
  num_rows<-dim(anova_res)
  pvalues_factors<-data.frame(t(anova_res["Pr(>F)"][-c(num_rows),]))
  
  return(list("mainpvalues"=pvalues_factors,"posthocfactor1"=posthoc$Factor1[,4]))
  
}

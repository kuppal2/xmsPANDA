diffexplmtwowayanova <-
function(dataA,covar.matrix=NA){
  dataA<-as.data.frame(dataA)
  
  dataA$Factor1<-as.factor(dataA$Factor1)
  dataA$Factor2<-as.factor(dataA$Factor2)
  
  #Fit a 2-way ANOVA model
  res<-aov(as.numeric(Response) ~ (Factor1) + (Factor2) + (Factor1) * (Factor2),data=dataA)
  
  #get ANOVA results with p-values
  anova_res<-anova(res)
  
  #Perform post-hoc comparisons using the Tukey Honestly Significant Differences (HSD) test
  posthoc <- TukeyHSD(x=res, conf.level=0.95,test=univariate())
  
  num_rows<-dim(anova_res)[1]
  
  pvalues_factors<-data.frame(t(anova_res["Pr(>F)"][-c(num_rows),]))
  
  names(pvalues_factors)<-rownames(anova_res)[-c(num_rows)]
  
  interact_res<-t(c(posthoc$Factor1[,4],posthoc$Factor2[,4],posthoc$'Factor1:Factor2'[,4]))
  
  colnames(interact_res)<-c(rownames(posthoc$Factor1),rownames(posthoc$Factor2),rownames(posthoc$'Factor1:Factor2'))
  
  #return resutls
  return(list("mainpvalues"=pvalues_factors,"posthoc"=interact_res))
  
  
}

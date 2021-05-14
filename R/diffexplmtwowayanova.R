diffexplmtwowayanova <-
function(dataA,covar.matrix=NA,plot.tukeyhsd=FALSE,var.name=NA){
  dataA<-as.data.frame(dataA)
  
  dataA$Factor1<-as.factor(dataA$Factor1)
  dataA$Factor2<-as.factor(dataA$Factor2)
  
  ref_group<-dataA$Factor1:dataA$Factor2
  
  ref_group<-ref_group[1]
  
  #Fit a 2-way ANOVA model
  res<-aov(as.numeric(Response) ~ (Factor1) + (Factor2) + (Factor1) * (Factor2),data=dataA)
  
  #get ANOVA results with p-values
  anova_res<-anova(res)
  
  #Perform post-hoc comparisons using the Tukey Honestly Significant Differences (HSD) test
  posthoc <- TukeyHSD(x=res, conf.level=0.95,test=univariate())
  
  p1<-posthoc
  p2=p1$`Factor1`[grep(rownames(p1$`Factor1`),pattern=paste("-",levels(dataA$Factor1)[1],"$",sep="")),]
  p1$`Factor1`<-p2
  p2=p1$`Factor2`[grep(rownames(p1$`Factor2`),pattern=paste("-",levels(dataA$Factor2)[1],"$",sep="")),]
  p1$`Factor2`<-p2
  p2=p1$`Factor1:Factor2`[grep(rownames(p1$`Factor1:Factor2`),pattern=paste("-",ref_group,"$",sep="")),]
  
  #rownames(p2)<-gsub(rownames(p2),pattern=paste("-",ref_group,"$",sep=""),replacement="")
  
  p1$`Factor1:Factor2`<-p2
  
  ref_groups<-paste("Differences in mean levels \ncompared to the reference group: ",c(levels(dataA$Factor1)[1],levels(dataA$Factor2)[1],levels(dataA$Factor1:dataA$Factor2)[1]),sep="")
  
  #plot(p1,cex.axis=0.5,las=2)
  
  plot_res<-try(plotTukeyHSD1(tukey.res=p1,x.axis.label = "",y.axis.label=ref_groups,var.name=var.name),silent=TRUE)
  
 # print(plot_res)
  num_rows<-dim(anova_res)[1]
  
  pvalues_factors<-data.frame(t(anova_res["Pr(>F)"][-c(num_rows),]))
  
  names(pvalues_factors)<-rownames(anova_res)[-c(num_rows)]
  
  interact_res<-t(c(posthoc$Factor1[,4],posthoc$Factor2[,4],posthoc$'Factor1:Factor2'[,4]))
  
  colnames(interact_res)<-c(rownames(posthoc$Factor1),rownames(posthoc$Factor2),rownames(posthoc$'Factor1:Factor2'))
  
  #return resutls
  return(list("mainpvalues"=pvalues_factors,"posthoc"=interact_res,"plot.tukeyhsd"=plot_res))
  
  
}

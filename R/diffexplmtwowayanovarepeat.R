diffexplmtwowayanovarepeat <-
function(dataA,subject_inf,modeltype="lme.RIRS",covar.matrix=NA){
  
  dataA<-as.data.frame(dataA)
  
  dataA$Factor1<-as.factor(dataA$Factor1)
  dataA$Factor2<-as.factor(dataA$Factor2)
  
  subject_inf<-as.vector(subject_inf)
  
  Subject<-subject_inf
  dataA<-cbind(dataA,subject_inf)
  
  ##save(dataA,file="dataA.Rda")
  
  if(modeltype=="lme.RI"){
    
    
    #call the lme function from the nlme package; random intercept only model
    res <- lme(as.numeric(Response) ~ Factor1 + Factor2 + Factor1 * Factor2, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf))
  }else{
    if(modeltype=="lme.RIRS"){
      
      #call the lme function from the nlme package;random intercept and random slope model
      res <- lme(as.numeric(Response) ~ Factor1 + Factor2 + Factor1 * Factor2, random=~ 1 + Factor2 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #corAR1())
      
    }else{
      if(modeltype=="nlme.RIRS"){
        
        #call the lme function from the nlme package;random intercept and random slope model
        res <- nlme(as.numeric(Response) ~ Factor1 + Factor2 + Factor1 * Factor2, random=~ 1 + Factor2 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #corAR1())
        
      }else{
        if(modeltype=="nlme.RI"){
          
          #call the lme function from the nlme package; random intercept only model
          res <- nlme(as.numeric(Response) ~ Factor1 + Factor2 + Factor1 * Factor2, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf))
          
        }
        
      }
      
    }
    
  }
  
  if(is(res,"try-error")){
    
    return(list("mainpvalues"=NA,"posthoc"=NA))
  }else{
    
    anova_res<-anova(res)
    
    
    
    #using lsmeans package for post-hoc comparisons since version v1.0.7.6
    means.factors=lsmeans(res,specs=c("Factor1","Factor2"))
    posthoc_res=pairs(means.factors,adjust="tukey")
    posthoc_res<-data.frame(posthoc_res)
    posthoc_pvalues<-posthoc_res$p.value
    names(posthoc_pvalues)<-as.character(posthoc_res$contrast)
    
    num_rows<-dim(anova_res)
    pvalues_factors<-data.frame(t(anova_res["p-value"][-c(1),]))
    
    names(pvalues_factors)<-rownames(anova_res)[-c(1)]
    
    return(list("mainpvalues"=pvalues_factors,"posthoc"=posthoc_pvalues))
  }
  
}

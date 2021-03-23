diffexplmonewayanovarepeat <-
function(dataA,subject_inf,analysismode="classification",modeltype="RIRS",covar.matrix=NA){
  
  dataA<-as.data.frame(dataA)
  
  if(analysismode=="classification"){
    dataA$Factor1<-as.factor(dataA$Factor1)
  }
  Subject<-subject_inf
  dataA<-cbind(dataA,subject_inf)
  
  ##savedataA,file="dataA.Rda")
if(FALSE){ 
  if(modeltype=="RI"){
    
    res <- lme(as.numeric(Response) ~ Factor1, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #,silent=TRUE)
  }else{
    if(modeltype=="RIRS"){
      res <- try(lme(as.numeric(Response) ~ Factor1, random = ~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)),silent=TRUE)
    }
  }
}  
  if(modeltype=="lme.RI"){
    
    res <- lme(as.numeric(Response) ~ Factor1, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #,silent=TRUE)
  }else{
    if(modeltype=="lme.RIRS"){
      res <- try(lme(as.numeric(Response) ~ Factor1, random = ~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)),silent=TRUE)
    }else{
      if(modeltype=="nlme.RIRS"){
        res <- try(nlme(as.numeric(Response) ~ Factor1, random = ~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)),silent=TRUE)
      }else{
        if(modeltype=="nlme.RI"){
          res <- nlme(as.numeric(Response) ~ Factor1, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #,silent=TRUE)
        }
      }
    }
    
    
  }
  
  if(is(res,"try-error")){
    return(list("mainpvalues"=NA,"posthoc"=NA))
  }else{
    anova_res<-anova(res)
    num_rows<-dim(anova_res)
    pvalues_factors<-data.frame(t(anova_res["p-value"][-c(1),]))
    
    if(analysismode=="classification"){
      
      
      
      #using lsmeans package for post-hoc comparisons since version v1.0.7.6
      means.factors=lsmeans(res,specs=c("Factor1"))
      posthoc_res=pairs(means.factors,adjust="tukey")
      posthoc_res<-data.frame(posthoc_res)
      posthoc_pvalues<-posthoc_res$p.value
      names(posthoc_pvalues)<-as.character(posthoc_res$contrast)
      
      return(list("mainpvalues"=pvalues_factors,"posthoc"=posthoc_pvalues))
    }else{
      return(list("mainpvalues"=pvalues_factors))
    }
  }
  
}

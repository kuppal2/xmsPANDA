diffexplmregrepeat <-
function(dataA,subject_inf,modeltype="lme.RIRS",covar.matrix=NA){
  
  dataA<-as.data.frame(dataA)
  
  Subject<-subject_inf
  dataA<-cbind(dataA,subject_inf)
  
  if(is.na(covar.matrix)==FALSE){
          if(modeltype=="lme.RI"){
            
            res <- lme(as.numeric(Response) ~ as.numeric(Factor1) + covar.matrix, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #,silent=TRUE)
          }else{
            if(modeltype=="lme.RIRS"){
              res <- try(lme(as.numeric(Response) ~ as.numeric(Factor1) + covar.matrix, random = ~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)),silent=TRUE)
            }else{
              if(modeltype=="nlme.RIRS"){
                res <- try(nlme(as.numeric(Response) ~ as.numeric(Factor1) + covar.matrix, random = ~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)),silent=TRUE)
              }else{
                if(modeltype=="nlme.RI"){
                  res <- nlme(as.numeric(Response) ~ as.numeric(Factor1) + covar.matrix, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #,silent=TRUE)
                }
              }
            }
            
            
          }
  }else{
    
    if(modeltype=="lme.RI"){
      
      res <- lme(as.numeric(Response) ~ as.numeric(Factor1), random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #,silent=TRUE)
    }else{
      if(modeltype=="lme.RIRS"){
        res <- try(lme(as.numeric(Response) ~ as.numeric(Factor1), random = ~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)),silent=TRUE)
      }else{
        if(modeltype=="nlme.RIRS"){
          res <- try(nlme(as.numeric(Response) ~ as.numeric(Factor1), random = ~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)),silent=TRUE)
        }else{
          if(modeltype=="nlme.RI"){
            res <- nlme(as.numeric(Response) ~ as.numeric(Factor1), random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #,silent=TRUE)
          }
        }
      }
      
      
    }
    
    
  }
  if(is(res,"try-error")){
    
    return(NA)
  }else{
    s1<-summary(res)
    
    s2=s1$tTable[2,c("p-value","Value","Std.Error","t-value")]
    return(s2)
  }
  
}

diffexplmregrepeat <-
function(dataA,subject_inf,modeltype="RIRS",covar.matrix=NA){
  
  dataA<-as.data.frame(dataA)
  
  Subject<-subject_inf
  dataA<-cbind(dataA,subject_inf)
  
  
  
  if(modeltype=="RI"){
    
    res <- lme(as.numeric(Response) ~ Factor1, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)) #,silent=TRUE)
  }else{
    if(modeltype=="RIRS"){
      res <- try(lme(as.numeric(Response) ~ Factor1, random = ~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim"),correlation=corCompSymm(form=~1|subject_inf)),silent=TRUE)
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

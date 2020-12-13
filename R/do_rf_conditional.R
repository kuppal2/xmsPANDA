do_rf_conditional <-
function(X,classlabels, ntrees=1000, analysismode="classification"){
  
  #classlabels<-c(rep("A",6),rep("B",6))
  
  dataA<-t(X)
  dataA<-as.data.frame(dataA)
  dataA<-cbind(classlabels,dataA)
  dataA<-as.data.frame(dataA)
  
  colnames(dataA)<-c("Targetvar",paste("mz",seq(1,dim(X)[1]),sep=""))
  
  dataA<-as.data.frame(dataA)
  
  attach(dataA,warn.conflicts=FALSE)
  
  mtry_num<-if (!is.null(classlabels) && !is.factor(classlabels))
    max(floor(ncol(X)/3), 1) else floor(sqrt(ncol(X)))
  
  #return(dataA)
  set.seed(290875)
  
  nvar=dim(X)[1] 
  if(analysismode=="classfication"){
    set.seed(290875)
    rf1<-cforest(factor(Targetvar)~.,data=dataA,control=cforest_control(teststat = "quad", testtype = "Univ",mincriterion = 0,savesplitstats = FALSE,ntree = ntrees, mtry = ceiling(sqrt(nvar)), replace = TRUE,fraction = 0.632, trace = FALSE))
    
  }else{
    set.seed(290875)
    #rf1<-cforest(Targetvar~.,data=dataA,control=cforest_control(teststat = "max",testtype = "Teststatistic",mincriterion = qnorm(0.9),savesplitstats = FALSE,ntree = ntrees, mtry = NULL, replace = TRUE,fraction = 0.632, trace = FALSE))
    
    #rf1<-cforest(Targetvar~.,data=dataA,control=cforest_unbiased(mtry = ceiling(sqrt(nvar)), ntree =ntrees))
    #  ntree = ntrees, mtry = NULL, replace = TRUE,
    #fraction = 0.632, trace = FALSE))
    
    rf1<-cforest(factor(Targetvar)~.,data=dataA,control=cforest_control(teststat = "quad", testtype = "Univ",mincriterion = 0,savesplitstats = FALSE,ntree = ntrees, mtry = ceiling(sqrt(nvar)), replace = TRUE,fraction = 0.632, trace = FALSE))
    
  }
  
  if(analysismode=="classfication"){
    set.seed(290875)
    varimp_res<-varimp(rf1,conditional=TRUE) #varimpAUC(rf1,conditional=TRUE)
  }else{
    varimp_res<-varimp(rf1,conditional=TRUE)
  }
  print(varimp_res)
  rm(dataA) 
  
  #return(varimp_res)
  return(list("rf_model"=rf1,"rf_varimp"=varimp_res))
  
}

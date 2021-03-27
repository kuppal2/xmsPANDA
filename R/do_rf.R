do_rf <-
function(X,classlabels, ntrees=1000, analysismode="classification"){
  dataA<-t(X)
  dataA<-as.data.frame(dataA)
  classlabels<-as.data.frame(classlabels)
  dataA<-cbind(classlabels,dataA)
  dataA<-as.data.frame(dataA)
  
  colnames(dataA)<-c("Targetvar",paste("mz",seq(1,dim(X)[1]),sep=""))
  dataA<-as.data.frame(dataA)
  
  attach(dataA,warn.conflicts=FALSE)
  
  #	return(dataA)	
  
  set.seed(290875)
  if(analysismode=="classfication"){
    set.seed(290875)
    rf2 <- randomForest(as.factor(Targetvar) ~ .,data=dataA, importance=TRUE,proximity=FALSE,ntree=ntrees,keep.forest=FALSE)
  }else{
    set.seed(290875)
    rf2 <- randomForest(Targetvar ~ .,data=dataA, importance=TRUE,proximity=FALSE,ntree=ntrees,keep.forest=FALSE)
  }
  
  #based on permutation importance;
  varimp_res<-randomForest::importance(rf2,type=1,scale=FALSE)
  varimp_res_scaled<-randomForest::importance(rf2,type=1,scale=TRUE)
  rm(dataA)
  
  #return(varimp_res)
  return(list("rf_model"=rf2,"rf_varimp"=varimp_res,"rf_varimp_scaled"=varimp_res_scaled))
}

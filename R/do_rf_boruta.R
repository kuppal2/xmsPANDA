do_rf_boruta <-
function(X,classlabels,maxRuns=1000){
  
  temp1<-cbind(classlabels,t(X))
  temp1<-as.data.frame(temp1)
  set.seed(290875)
  Boruta(classlabels~.,data=temp1,maxRuns=maxRuns)->Bor.test
  mean_imphistory<-apply(Bor.test$ImpHistory,2,mean)
  varimp_res2<-mean_imphistory[1:(length(mean_imphistory)-3)] #rep(0,length(Bor.test$finalDecision))
  temp_decision_vec<-as.character(Bor.test$finalDecision)
  if(length(which(temp_decision_vec=="Rejected"))>0){
    temp_decision_vec<-replace(temp_decision_vec,which(temp_decision_vec=="Rejected"),0)
  }
  if(length(which(temp_decision_vec=="Tentative"))>0){
    temp_decision_vec<-replace(temp_decision_vec,which(temp_decision_vec=="Tentative"),1)
  }
  if(length(which(temp_decision_vec=="Confirmed"))>0){
    temp_decision_vec<-replace(temp_decision_vec,which(temp_decision_vec=="Confirmed"),2)
  }
  temp_decision_vec<-as.numeric(as.character(temp_decision_vec))
  
  if(length(which(temp_decision_vec<2))>0){
    varimp_res2[which(temp_decision_vec<2)]<-0
  }
  return(varimp_res2)
}

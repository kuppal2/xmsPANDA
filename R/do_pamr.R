do_pamr <-
function(X,Y,fdrthresh=0.1,nperms=100,pamr.threshold.select.max=FALSE,kfold=10){
  
  
  library(pamr)
  
  #d1<-list(x=X,y=Y)
  
  d1 <- list(x=as.matrix(X),y=factor(Y[,1]), geneid=as.character(1:nrow(X)),
             genenames=paste("g",as.character(1:nrow(X)),sep=""))
  
  ###save(d1,file="d1.Rda")
  p1<-pamr.train(d1)
  
  set.seed(999)
  p2<-pamr.cv(data=d1,fit=p1,nfold=kfold)
  
  threshold_CVerror_matrix<-cbind(p2$threshold,p2$error)
  threshold_CVerror_matrix<-as.data.frame(threshold_CVerror_matrix)
  colnames(threshold_CVerror_matrix)<-c("Threshold","error")
  
  threshold_value<-max(threshold_CVerror_matrix[which(threshold_CVerror_matrix[,2]==min(threshold_CVerror_matrix[,2])),1])[1]
  
  
  
  set.seed(999)
  p3<-pamr.fdr(data=d1,p1,nperms=nperms)
  
  selected_feature_index<-{}
  
  max.discore<-{}
  
  #use median FDR
  if(length(which(p3$results[,4]<fdrthresh))>0){
    pamr_fdr_filt<-p3$results[which(p3$results[,4]<fdrthresh),]
    
    if(length((pamr_fdr_filt))>0){
      threshold_CVerror_fdrmatrix<-merge(threshold_CVerror_matrix,pamr_fdr_filt,by="Threshold")
      
      if(pamr.threshold.select.max==TRUE){
        threshold_value<-max(threshold_CVerror_fdrmatrix[which(threshold_CVerror_fdrmatrix[,2]==min(threshold_CVerror_fdrmatrix[,2])),1])[1]
      }else{
        threshold_value<-min(threshold_CVerror_fdrmatrix[which(threshold_CVerror_fdrmatrix[,2]==min(threshold_CVerror_fdrmatrix[,2])),1])[1]
        
        
      }
      
    }
  }else{
    
    print("pamr: no features meet the fdr criteria")
  }
  
  p4<-pamr.listgenes(fit=p1,data=d1,threshold=threshold_value)
  p4<-as.data.frame(p4)
  
  selected_feature_index<-p4$id
  selected_feature_index<-as.numeric(as.character(selected_feature_index))
  
  discore_matrix<-p4[,-c(1)]
  if(nrow(discore_matrix)>1){
    discore_matrix<-apply(discore_matrix,2,as.numeric)
    abs.discore_matrix<-abs(discore_matrix)
    max.discore<-apply(abs.discore_matrix,1,max)
    
  }else{
    discore_matrix_1<-unlist(discore_matrix)
    discore_matrix<-as.numeric(as.character(discore_matrix_1))
    abs.discore_matrix<-abs(discore_matrix)
    max.discore<-max(abs.discore_matrix)
    
  }
  
  
  
  
  pall<-pamr.listgenes(fit=p1,data=d1,threshold=0)
  
  ##savepall,file="pall.Rda")
  
  pall<-as.data.frame(pall)
  
  
  
  discore_matrix_all<-pall #[,-c(1)]
  discore_matrix_all<-apply(discore_matrix_all,2,as.numeric)
  discore_matrix_all<-as.data.frame(discore_matrix_all)
  discore_matrix_all<-discore_matrix_all[order(discore_matrix_all$id),]
  
  
  discore_matrix_all<-discore_matrix_all[,-c(1)]
  if(nrow(discore_matrix_all)>1){
    discore_matrix_all<-apply(discore_matrix_all,2,as.numeric)
    abs.discore_matrix_all<-abs(discore_matrix_all)
    max.discore.all<-apply(abs.discore_matrix_all,1,max)
    
    max.discore.all.thresh<-min(max.discore.all[selected_feature_index],na.rm=TRUE)
    
    
  }else{
    discore_matrix_all_1<-unlist(discore_matrix_all)
    discore_matrix_all<-as.numeric(as.character(discore_matrix_all_1))
    abs.discore_matrix_all<-abs(discore_matrix_all)
    max.discore.all<-max(abs.discore_matrix_all)
    
    max.discore.all.thresh<-min(max.discore.all[selected_feature_index],na.rm=TRUE)
    
  }
  
  
  # ###savelist=ls(),file="debug.Rda")
  return(list("feature.list"=selected_feature_index,"max.discore.sigfeats"=max.discore,"pam_train_model"=p1,"pam_toplist"=p4,"max.discore.allfeats"=max.discore.all,"threshold_value"=threshold_value,"max.discore.all.thresh"=max.discore.all.thresh))
}

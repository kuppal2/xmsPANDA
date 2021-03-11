plsda_cv <-
function(v,x,y,ncomp,errortype="total",conflevel=99){
  
  num_samp=dim(x)[1]
  
  num_datasets= floor(num_samp)
  n1<-floor(num_samp/v)
  n2<-num_samp-n1*v
  n3<-v-n2
  
  ind<-rep(c(n1,n1+1),c(n3,n2))
  ind<-diffinv(ind)
  min_err=1
  best_k=1
  
  set.seed(555)
  group<-sample(1:num_samp,num_samp, replace=FALSE)
  
  
  itr=0
  #plsda_error <- matrix(0,v)  # we set K=30 before, it can be changed to any number<100.
  plsda_error<-rep(0,v)
  for ( i in 1:v)
  {
    g<-group[(ind[i]+1):ind[i+1]]
    temptest<-x[g,]
    temptrain <-x[-g,]
    tempclass <-y[-g]
    testclass<-y[g]
    
    #temptest<-as.data.frame(temptest)
    #temptrain<-as.matrix(temptrain)
    
    
    
    #print(dim(temptrain))
    #print(dim(temptest))
    
    #plsda_cv <- plsda(x=temptrain,y=tempclass, type="nu-classification",kernel=kname) 
    
    #opt_comp<-pls.lda.cv(Xtrain=temptrain, Ytrain=tempclass,  ncomp=c(1:10), nruncv=10, alpha=2/3, priors=NULL)
    predfit<-pls.lda(Xtrain=temptrain,Ytrain=tempclass,ncomp=ncomp,nruncv=v,Xtest=temptest)
    #predfit<-predict(plsda_pred,temptest)
    
    #print(length(which(plsda_pred$predclass==testclass)))
    svm_table<-table(predfit$predclass,testclass)
    
    class_names<-rownames(svm_table)
    #print(testclass)
    #print(predfit$predclass)
    predfit<-predfit$predclass
    totacc<-length(which(predfit==testclass))/length(testclass)
    
    beracc<-{}
    auc_acc<-{}
    for(c in 1:dim(svm_table)[1]){
      testclass_ind<-which(testclass==class_names[c])
      beracc<-c(beracc,length(which(predfit[testclass_ind]==testclass[testclass_ind]))/length(testclass_ind))
      
      
      if(errortype=="AUC"){
        pred_acc<-multiclass.roc(testclass,as.numeric(predfit),levels=levels(as.factor(y)))
        pred_acc_orig<-pred_acc$auc[1]
        auc_acc<-c(auc_acc,pred_acc_orig)
      }
      
    }
    
    beracc<-mean(beracc,na.rm=TRUE)
    
    if(errortype=="CV"){
      plsda_error[i]<-(totacc*100)	
    }else{
      if(errortype=="AUC"){
        plsda_error[i]<-(auc_acc*100)
      }else{
        plsda_error[i]<-(beracc*100)
      }
    }
    
    
    
  }
  avgacc <-mean(plsda_error)
  sdacc<-sd(plsda_error)
  
  probval<-(1-(conflevel*0.01))/2
  probval<-1-probval
  
  error <- qnorm(probval)*sdacc/sqrt(length(y))
  
  leftconfint<-avgacc-error
  rightconfint<-avgacc+error
  
  
  
  return(list(mean_acc=avgacc,sd_acc=sdacc, acc_each_fold=plsda_error,confint=c(leftconfint,rightconfint)))
  
  
}

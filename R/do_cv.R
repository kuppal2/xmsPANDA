do_cv <-
function(v,x,y,kname="radial",errortype="CV",conflevel=95,classifier="SVM",seednum=555){
  
  
  num_samp=dim(x)[1]
  
  if(length(num_samp)<1){
    num_samp<-length(x)
    x<-as.data.frame(x)
  }    
  y<-as.data.frame(y)
  
  num_datasets= floor(num_samp)
  n1<-floor(num_samp/v)
  n2<-num_samp-n1*v
  n3<-v-n2
  
  ind<-rep(c(n1,n1+1),c(n3,n2))
  ind<-diffinv(ind)
  min_err=1
  best_k=1
  
  set.seed(seednum)
  group<-sample(1:num_samp,num_samp, replace=FALSE)
  
  
  itr=0
  #svm_acc <- matrix(0,v)  # we set K=30 before, it can be changed to any number<100.
  svm_acc<-rep(0,v)
  for ( i in 1:v)
  {
    g<-group[(ind[i]+1):ind[i+1]]
    temptest<-x[g,]
    temptrain <-x[-g,]
    tempclass <-y[-g,]
    testclass<-y[g,]
    
    
    if(classifier=="SVM"){
      mod_cv <- svm(x=temptrain,y=tempclass, type="nu-classification",kernel=kname)
      
      
      if(v==num_samp){
        
        predfit<-predict(mod_cv,(temptest))
        
      }else{
        predfit<-predict(mod_cv,(temptest))
      }
      
      
    }else{
      if(classifier=="LR"){
        
        Class<-tempclass
        
        dtemp<-cbind(Class,temptrain)
        mod_cv<-glm(as.factor(Class)~.,data=dtemp,family=binomial(logit))
        
        if(v==num_samp){
          
          #predfit<-predict(mod_cv,t(temptest))
          
          predfit<-predict(mod_cv,t(temptest),type="response")
        }else{
          # predfit<-predict(mod_cv,(temptest))
          
          predfit<-predict(mod_cv,temptest,type="response")
        }
        
        
        #mod_cv<-glm.fit(x=temptrain,y=tempclass,family="binomial")
        #predfit<-predict(mod_cv,temptest,type="response")
        
        #print(predfit)
        
        predfit <- ifelse(predfit > 0.5,1,0)
        
        
      }else{
        
        
        if(classifier=="RF"){
          
          
          Class<-tempclass
          
          d1<-cbind(Class,temptrain)
          mod_cv<-randomForest(as.factor(Class)~.,data=d1)
          #mod_cv<-randomForest(x=temptrain,y=tempclass)
          
          # predfit<-predict(mod_cv,temptest)
          if(v==num_samp){
            
            predfit<-predict(mod_cv,t(temptest))
          }else{
            predfit<-predict(mod_cv,(temptest))
          }
          
          
        }
        
      }
      
    }
    
    
    svm_table<-table(predfit,testclass)
    
    class_names<-rownames(svm_table)
    beracc<-{}
    auc_acc<-{}
    totacc<-length(which(predfit==testclass))/length(testclass)
    for(c in 1:dim(svm_table)[1]){
      testclass_ind<-which(testclass==class_names[c])
      beracc<-c(beracc,length(which(predfit[testclass_ind]==testclass[testclass_ind]))/length(testclass_ind))
      
    }
    if(errortype=="AUC"){
      testclass<-as.vector(testclass)
      y1<-as.vector(y[,1])
      pred_acc<-multiclass.roc(testclass,as.numeric(predfit),levels=levels(as.factor(y1)))
      pred_acc_orig<-pred_acc$auc[1]
      auc_acc<-c(auc_acc,pred_acc_orig)
    }
    
    
    
    beracc<-mean(beracc,na.rm=TRUE)
    
    if(errortype=="CV"){
      svm_acc[i]<-(totacc*100)
    }else{
      if(errortype=="AUC"){
        svm_acc[i]<-(auc_acc*100)
      }else{
        svm_acc[i]<-(beracc*100)
      }
    }
    
  }
  avg_acc <-mean(svm_acc,na.rm=TRUE)
  sd_acc<-sd(svm_acc,na.rm=TRUE)
  
  #limit<-avg_acc-(sd.error*(avg_acc) # 1 sd criterion
  #print(avg_acc)
  #print(sd_acc)
  
  #return(list(error=avg_acc,sderror=sd.error))
  probval<-(1-(conflevel*0.01))/2
  probval<-1-probval
  #print(probval)
  error <- qnorm(probval)*sd_acc/sqrt(length(svm_acc))
  
  leftconfint<-avg_acc-error
  rightconfint<-avg_acc+error
  
  options(warn=0)
  #print("done")
  return(list(avg_acc=avg_acc,sd_acc=sd_acc, acc_each_fold=svm_acc,confint=c(leftconfint,rightconfint)))
  #return(list(num=best_k,error=min_err, avg=avg_acc))
}

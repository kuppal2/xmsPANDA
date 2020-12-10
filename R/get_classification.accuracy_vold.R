get_classification.accuracy_vold <-
function(kfold,featuretable,classlabels,kernelname="radial",errortype="AUC",conflevel=95,classifier="svm",seednum=555,testfeaturetable=NA,testclasslabels=NA,match_class_dist=TRUE,plotroc=TRUE,svm.cost=NA,svm.gamma=NA,opt_comp=NA,column.rm.index=c(1,2),svm.type="nu-classification"){
  
  
  
  v=kfold
  kname=kernelname
  
  options(warn=-1)
  
  featuretable<-as.data.frame(featuretable)
  
  
  mz_time<-rownames(featuretable)
  
  if(length(mz_time)<1){
    
    mz_time<-paste(round(featuretable[,1],5),"_",round(featuretable[,2],2),sep="")
  }
  
  if(is.na(column.rm.index)==FALSE){
    x=t(featuretable[,-c(column.rm.index)])
  }else{
    
    x=t(featuretable)
  }
  if(is.na(testfeaturetable)==FALSE){
    
    if(nrow(featuretable)!=nrow(testfeaturetable)){
      
      stop("Number of features/variables should be same in the train and test sets.")
    }
    
    print("Note: the order of features/variables should be same in the train and test sets.")
    
    if(is.na(column.rm.index)==FALSE){
      testfeaturetable<-t(testfeaturetable[,-c(column.rm.index)])
    }else{
      
      testfeaturetable<-t(testfeaturetable)
    }
    colnames(testfeaturetable)<-mz_time
    if(is.na(testclasslabels)==FALSE){
      
      if(is.vector(testclasslabels)==TRUE){
        testclasslabels=as.data.frame(testclasslabels)
      }else{
        testclasslabels=as.data.frame(testclasslabels)
        
        if(dim(testclasslabels)[2]>1){
          testclasslabels<-testclasslabels[,2]
        }else{
          testclasslabels<-testclasslabels[,1]
        }
      }
      
    }
  }
  
  
  
  colnames(x)<-mz_time
  if(is.vector(classlabels)==TRUE){
    y=as.data.frame(classlabels)
  }else{
    classlabels=as.data.frame(classlabels)
    y=classlabels
    if(dim(y)[2]>1){
      y=classlabels[,2]
    }else{
      y=classlabels[,1]
    }
  }
  y=as.data.frame(y)
  
  
  num_samp=dim(x)[1]
  
  if(length(num_samp)<1){
    num_samp<-length(x)
    x<-as.data.frame(x)
  }
  y<-as.data.frame(y)
  x<-as.data.frame(x)
  
  num_datasets= floor(num_samp)
  n1<-floor(num_samp/v)
  n2<-num_samp-n1*v
  n3<-v-n2
  
  ind<-rep(c(n1,n1+1),c(n3,n2))
  ind<-diffinv(ind)
  min_err=1
  best_k=1
  
  suppressMessages(library(CMA))
  if(match_class_dist==TRUE){
    set.seed(seednum)
    #group=GenerateLearningsets(y=y[,1],niter=1,fold=v,method="CV",strat=TRUE)
    
    group=GenerateLearningsets(y=y[,1],niter=v,fold=v,method="MCCV",strat=TRUE,ntrain=(num_samp-floor(num_samp/v)))
    #return(y)
    #return(group)
    
  }else{
    
    set.seed(seednum)
    group<-sample(1:num_samp,num_samp, replace=FALSE)
    
  }
  
  
  itr=0
  
  # ##savev,x,y,classifier,plotroc,errortype,kernelname,group,file="temp.Rda")
  
  
  
  svm_acc<-rep(0,v)
  mod_cv_list<-new("list")
  svm_acc<-lapply(1:v,function(i)
  {
    if(match_class_dist==TRUE){
      g<-group@learnmatrix[i,]
      g<-seq(1,length(y[,1]))[-group@learnmatrix[i,]]
    }else{
      g<-group[(ind[i]+1):ind[i+1]]
    }
    
    temptest<-x[g,]
    temptrain <-x[-g,]
    tempclass <-y[-g,]
    testclass<-y[g,]
    
    # ##save(g,temptrain,temptest,tempclass,testclass,v,classifier,num_samp,errortype,kernelname,svm.type,file="temp.Rda")
    cv_res<-get_classification.accuracy.child(temptrain=temptrain,tempclass=tempclass,kernelname=kernelname,errortype=errortype,classifier=classifier,num_samp=num_samp,temptest=temptest,testclass=testclass,numfolds=v,plotroc=FALSE,rocfeatlist=NA,svm.cost=svm.cost,svm.gamma=svm.gamma,svm.type=svm.type)
    
    #svm_acc[i]<-cv_res$classification_acc
    
    return(cv_res$classification_acc)
  })
  
  
  svm_acc<-unlist(svm_acc)
  
  
  avg_acc <-mean(svm_acc,na.rm=TRUE)
  sd_acc<-sd(svm_acc,na.rm=TRUE)
  
  ##Get confidence interval
  probval<-(1-(conflevel*0.01))/2
  probval<-1-probval
  
  error <- qnorm(probval)*sd_acc/sqrt(length(svm_acc))
  avg_acc<-round(avg_acc,2)
  leftconfint<-avg_acc-error
  rightconfint<-avg_acc+error
  test_acc<-NA
  test_confusion_matrix<-NA
  
  leftconfint<-round(leftconfint,2)
  rightconfint<-round(rightconfint,2)
  print(paste("Training set ", kfold,"-fold CV ",errortype," ",classifier," classification accuracy (%):",avg_acc,sep=""))
  print(paste("Training set ", kfold,"-fold CV ",errortype," ",classifier," classification accuracy ", conflevel,"% confidence interval:(",leftconfint,",",rightconfint,")",sep=""))
  
  
  x<-as.data.frame(x)
  y<-y[,1]
  ###save(x,y,kernelname,errortype,classifier,num_samp,svm.cost,svm.gamma,file="debugclass.Rda")
  train_res<-get_classification.accuracy.child(temptrain=x,tempclass=y,kernelname=kernelname,errortype=errortype,classifier=classifier,num_samp=num_samp,temptest=x,testclass=y,svm.cost=svm.cost,svm.gamma=svm.gamma)
  print(paste("Training set ",errortype," ",classifier," classification accuracy (%):",round(train_res$classification_acc,2),sep=""))
  
  mod_cv<-train_res$classification_model
  
  test_res={}
  #evaluate test set accuracy
  if(is.na(testfeaturetable)==FALSE){
    
    if(is.na(testclasslabels)==FALSE){
      
      testfeaturetable<-as.data.frame(testfeaturetable)
      
      # ##save(list=ls(),file="t2.Rda")
      
      
      
      test_res<-get_classification.accuracy.child(temptrain=x,tempclass=y,kernelname=kernelname,errortype=errortype,classifier=classifier,num_samp=num_samp,temptest=testfeaturetable,testclass=testclasslabels,plotroc=TRUE,svm.cost=svm.cost,svm.gamma=svm.gamma)
      
      test_acc<-test_res$classification_acc
      test_acc<-round(test_acc,2)
      test_pred_table<-test_res$confusion_matrix
      
      print(paste("Test set ", errortype," ",classifier," classification accuracy (%):",test_acc,sep=""))
      print(paste("Test set confusion matrix using ",classifier,sep=""))
      print(test_pred_table)
      
      test_confusion_matrix<-test_pred_table
      
    }
  }
  
  options(warn=0)
  
  if(classifier=="logitreg" | classifier=="LR"){
    
    Class<-y
    
    dtemp<-cbind(Class,x)
    dtemp<-as.data.frame(dtemp)
    
    mod_cv_all<-glm(as.factor(Class)~.,data=dtemp,family=binomial(logit))
    s1<-summary(mod_cv_all)
    
    options(warn=0)
    return(list(avg.train.cv.acc=avg_acc,sd.train.cv.acc=sd_acc, train.cv.acc.each.fold=svm_acc,glm_fit=s1$coefficients,train.cv.acc.confint=c(leftconfint,rightconfint),test.acc=test_acc,test.confusion.matrix=test_confusion_matrix,test_res=test_res))
  }else{
    
    options(warn=0)
    return(list(avg.train.cv.acc=avg_acc,sd.train.cv.acc=sd_acc,train.cv.acc.each.fold=svm_acc,train.cv.acc.confint=c(leftconfint,rightconfint),test.acc=test_acc,test.confusion.matrix=test_confusion_matrix,test_res=test_res))
  }
  
  
}

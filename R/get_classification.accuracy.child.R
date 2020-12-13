get_classification.accuracy.child <-
function(temptrain,tempclass,kernelname="radial",errortype="AUC",classifier="svm",num_samp,temptest=NA,testclass=NA,numfolds=1,plotroc=TRUE,rocfeatlist=NA,svm.cost=NA,svm.gamma=NA,opt_comp=NA,svm.type="nu-classification",confint.auc=TRUE)
{
  
  suppressMessages(library(pROC))
  suppressMessages(library(ROCR))
  
  options(warn=-1)
  v=numfolds
  classifier=tolower(classifier)
  
  auc_val=NA
  
  #  ##save(temptrain,tempclass,svm.type,kernelname,file="temp2.Rda")
  
  if(classifier=="SVM" | classifier=="svm"){
    suppressMessages(library(e1071))
    if(is.na(svm.cost)==TRUE & is.na(svm.gamma)==TRUE){
      mod_cv <- svm(x=temptrain,y=tempclass, type=svm.type,kernel=kernelname)
      
    }else{
      
      mod_cv <- svm(x=temptrain,y=tempclass, type=svm.type,kernel=kernelname,gammma=svm.gamma,cost=svm.cost)
    }
    
    
    
    if(v==num_samp){
      
      predfit<-predict(mod_cv,(temptest))
      
    }else{
      predfit<-predict(mod_cv,(temptest))
    }
    
    #  ##save(predfit,file="predfit.Rda")
    
    
  }else{
    if(classifier=="logitreg" | classifier=="LR" | classifier=="lr"){
      
      Class<-as.numeric(tempclass) #-1
      
      
      labels_1<-levels(as.factor(Class))
      
      Class<-replace(Class,which(Class==labels_1[1]),0)
      Class<-replace(Class,which(Class==labels_1[2]),1)
      
      dtemp<-cbind(Class,temptrain)
      dtemp<-as.data.frame(dtemp)
      
      labels_2<-levels(as.factor(testclass))
      
      testclass<-replace(testclass,which(testclass==labels_2[1]),0)
      testclass<-replace(testclass,which(testclass==labels_2[2]),1)
      
      
      
      mod_cv<-glm(as.factor(Class)~.,data=dtemp,family=binomial(logit))
      
      if(v==num_samp){
        
        predfit<-predict(mod_cv,(temptest),type="response")
      }else{
        
        
        predfit<-predict(mod_cv,temptest,type="response")
      }
      
      
      predfit <- ifelse(predfit > 0.5,1,0)
      
      #      ##save(dtemp,predfit,mod_cv,temptest,file="logit.Rda")
      # testclass<-as.numeric(testclass)-1
      # print(predfit)
      #print(testclass)
      
    }else{
      
      
      if(classifier=="RF" | classifier=="randomforest" | classifier=="rf" | classifier=="cforest"){
        
        
        Class<-tempclass
        
        
        if(classifier=="cforest"){
          suppressMessages(library(party))
          tempdata_c<-cbind(tempclass,temptrain)
          cnames_c<-colnames(tempdata_c)
          
          mod_cv <- cforest(tempclass~.,data=tempdata_c,controls=cforest_unbiased(ntree=5000))
          
          predfit <- predict(mod_cv, newdata=temptest, OOB=TRUE, type="response")
          
          predfit <- ifelse(predfit > 0.5,1,0)
          
          
          
        }else{
          suppressMessages(library(randomForest))
          
          mod_cv<-randomForest(x=(temptrain),y=as.factor(tempclass),ntree=5000)
          
          
          if(v==num_samp){
            
            predfit<-predict(mod_cv,(temptest))
          }else{
            predfit<-predict(mod_cv,(temptest))
          }
          
        }
        
        
        predfit<-as.numeric(as.character(predfit))
        
        
      }else{
        
        if(classifier=="NaiveBayes" | classifier=="naivebayes"){
          
          
          Class<-tempclass
          
          mod_cv<-naiveBayes(x=(temptrain),y=as.factor(tempclass))
          
          # predfit<-predict(mod_cv,temptest)
          if(v==num_samp){
            
            predfit<-predict(mod_cv,(temptest))
          }else{
            predfit<-predict(mod_cv,(temptest))
          }
          
          
        }else{
          
          if(classifier=="plsda" | classifier=="pls" | classifier=="pls.lda"){
            #set.seed(seednum)
            if(is.na(opt_comp)==TRUE){
              
              max_ncomp=min(c(dim(temptrain)[2],10))
              
              opt_comp<-pls.lda.cv(Xtrain=temptrain, Ytrain=tempclass,  ncomp=c(1:max_ncomp), nruncv=v, alpha=2/3, priors=NULL)
            }
            
            if(classifier=="pls.lda"){
              mod_cv<-pls.lda(Xtrain=temptrain,Ytrain=tempclass,ncomp=opt_comp,nruncv=v,Xtest=temptest)
              predfit<-mod_cv$predclass
              predfit<-as.numeric(as.character(predfit))
            }else{
              mod_cv<-mixOmics::plsda(X=temptrain,Y=(tempclass),ncomp = opt_comp)
              
              set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
              
              predfit<-predict(mod_cv,temptest,dist="mahalanobis.dist")
              
              predfit<-as.numeric(as.character(predfit$class$mahalanobis.dist[,1]))
            }
            
          }else{
            if(classifier=="plr" | classifier=="pLR"){
              
              Class<-tempclass
              
              
              param_1<-cv.step.plr(x=(temptrain),y=as.numeric(tempclass),lambda=c(1e-4, 1e-2, 0.1,0.5,1),nfold=v)
              
              mod_cv<-plr(x=(temptrain),y=as.numeric(tempclass)) #,lambda=1e-3)
              
              # predfit<-predict(mod_cv,temptest)
              if(v==num_samp){
                
                predfit<-predict(mod_cv,(temptest),type="class")
              }else{
                predfit<-predict(mod_cv,(temptest),type="class")
              }
              
              
              
              
            } #end
            
          }
          
        }
        
        
      }
      
    }
    
  }
  
  ###savelist=ls(),file="t1.Rda")
  
  
  svm_table<-table(predfit,testclass)
  # print(svm_table)
  # ##save(predfit,file="pred.Rda")
  ###save(testclass,file="testclass.Rda")
  
  conf.auc_acc<-{}
  
  
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
    
    pred_acc<-multiclass.roc(testclass,as.numeric(predfit),levels=levels(as.factor(testclass)))
    pred_acc_orig<-pred_acc$auc[1]
    auc_acc<-c(auc_acc,pred_acc_orig)
    
    if(confint.auc==TRUE){
      set.seed(555)
      conf.auc_acc<-ci.auc(testclass,predfit,method="bootstrap")
    }
  }
  
  
  roc_svm={}
  roc_lr={}
  roc_res<-{}
  
  
  if(dim(svm_table)[1]==2){
    
    if(plotroc==TRUE)
    {
      
      #print(pred)
      
      # ##savelist=ls(),file="t.Rda")
      temptrain<-t(temptrain)
      
      temptest<-t(temptest)
      
      ###save(temptrain,file="temptrain.Rda")
      #    ##save(temptest,file="temptest.Rda")
      ##save(testclass,file="testclass.Rda")
      #   #save(predfit,file="predfit.Rda")
      
      #roc_svm<-get_roc(dataA=temptrain,classlabels=tempclass,classifier="svm",kname="radial",rocfeatlist=seq(1,dim(temptrain)[1]),rocfeatincrement=FALSE,testset=temptest,testclasslabels=testclass,mainlabel="Test set (ROC);\n using SVM",mz_names=rownames(temptrain))
      
      #roc_lr<-get_roc(dataA=temptrain,classlabels=tempclass,classifier="logit",kname="radial",rocfeatlist=seq(1,dim(temptrain)[1]),rocfeatincrement=FALSE,testset=temptest,testclasslabels=testclass,mainlabel="Test set (ROC);\n using LR",mz_names=rownames(temptrain))
      predfit<-as.numeric(factor(predfit)) #as.numeric(as.character(predfit))
      
      
      
      if(length(levels(as.factor(testclass)))==2){
        
        pred1 <- ROCR::prediction(predfit, testclass)
        stats1a <- performance(pred1, 'tpr', 'fpr')
        
        roc_res<-cbind(stats1a@x.values[[1]],stats1a@y.values[[1]])
        colnames(roc_res)<-c(stats1a@x.name,stats1a@y.name)
        fname1=paste("ROC",classifier,".rda",sep="")
        ##saveroc_res,file=fname1)
        x1<-seq(0,1,0.01)
        y1<-x1
        p1<-performance(pred1,"auc")
        
        if(length(conf.auc_acc)>0){
          mod.lab <-c(classifier,paste(' performance, AUC: ',round(p1@y.values[[1]],2),"\n","95% CI:",round(conf.auc_acc[[1]]),"-",round(conf.auc_acc[[2]]),sep=""))
        }else{
          
          mod.lab <-c(classifier,paste(' performance, AUC: ',round(p1@y.values[[1]],2),sep=""))
          
        }
        auc_val=round(p1@y.values[[1]],2)
        
        #if(n==1){
        plot(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', ylab=stats1a@y.name, xlab=stats1a@x.name, col="black", lty=2, main=mod.lab,cex.main=0.7,lwd=2)
        #}else{
        #lines(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', col=col_lab[n], lty=2,lwd=2)
        #}
        
        
        
        
        
        
      }
      
    }
  }
  
  
  
  beracc<-mean(beracc,na.rm=TRUE)
  
  if(errortype=="CV" | errortype=="total"){
    svm_acc<-(totacc*100)
  }else{
    if(errortype=="AUC"){
      svm_acc<-(auc_acc*100)
    }else{
      svm_acc<-(beracc*100)
    }
  }
  options(warn=0)
  return(list(classification_acc=svm_acc,classification_model=mod_cv,confusion_matrix=svm_table,roc_res=roc_res,auc_res=auc_val,test.auc.confint=conf.auc_acc))
}

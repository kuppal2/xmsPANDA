get_roc <-
function(dataA,classlabels,classifier="svm",kname="radial",rocfeatlist=c(1,2,3,4,5),rocfeatincrement=TRUE,
                  testset=NA,testclasslabels=NA,mainlabel=NA,
                  col_lab=c("black"),legend=TRUE,newdevice=FALSE,mz_names=NA,svm.type="C-classification"){
  
  options(warn=-1)
  
  d1<-dataA
  rm(dataA)
  
  
  if(newdevice==TRUE){
    
    if(dir.exists("Figures")){
    
        pdf("Figures/ROC.pdf")
    }else{
      
      pdf("ROC.pdf")
    }
  }
  cnames<-colnames(d1)
  
  
  classlabels<-as.data.frame(classlabels)
  testclasslabels<-as.data.frame(testclasslabels)
  class_inf<-classlabels
  
 # save(classlabels,file="classlabels.Rda")
  #save(d1,mz_names,file="d1.Rda")
  
  if(is.na(rocfeatlist)==TRUE){
    num_select<-dim(d1)[1]
    rocfeatlist=1:num_select
  }else{
    
    num_select<-length(rocfeatlist)
  }
  
  
  if(is.na(mz_names)==TRUE){
    cnames[1]<-"mz"
    cnames[2]<-"time"
    colnames(d1)<-cnames
    d1<-as.data.frame(d1)
    
    #d1<-unique(d1)
    
    d2<-t(d1[,-c(1:2)])
    
    cnames<-colnames(d1)
    
    if(is.na(testset)==TRUE){
      testset<-d1
      testclasslabels<-classlabels
    }
    
    if(is.na(rocfeatlist)==TRUE){
      
      rocfeatlist<-seq(1,dim(d2)[2])
    }
    
    if(max(rocfeatlist,na.rm=TRUE)>ncol(d2))
    {
      rocfeatlist<-rocfeatlist[-which(rocfeatlist>ncol(d2))]
    }
    
    
    
    
    testset<-t(testset[,-c(1:2)])
    testset<-testset[,rocfeatlist]
    
    d2<-d2[,rocfeatlist]
    d1<-d1[rocfeatlist,]
    mz_names<-paste(d1$mz,d1$time,sep="_")
    
    
  }else{
    
    d1<-as.data.frame(d1)
    #d1<-unique(d1)
    
    d2<-t(d1)
    
    cnames<-colnames(d1)
    
    mz_names<-unique(mz_names)
    
    if(is.na(testset)==TRUE){
      testset<-d1
      testclasslabels<-classlabels
      
      #print(classlabels)
    }
    
    if(is.na(rocfeatlist)==TRUE){
      
      rocfeatlist<-seq(1,dim(d2)[2])
    }
    
    if(max(rocfeatlist,na.rm=TRUE)>ncol(d2))
    {
      rocfeatlist<-rocfeatlist[-which(rocfeatlist>ncol(d2))]
    }
    
    # testset<-unique(testset)
    testset<-t(testset)
    testset<-testset[,rocfeatlist]
    d2<-d2[,rocfeatlist]
    d1<-d1[rocfeatlist,]
    mz_names<-mz_names[rocfeatlist]
    
    
    
  }
  
  
  featlist<-rocfeatlist
  featincrement<-rocfeatincrement
  
  
  
  class_vec<-as.character(class_inf[,1])
  
  class_levels<-levels(as.factor(class_vec))
  
  if(length(class_levels)>2){
    
    #stop("More than two classes detected. Please restrict the analysis to only two classes.")
  }
  
	
  #class_vec<-replace(class_vec,which(class_vec==class_levels[2]),1)
  
  if(classifier=="logitreg" || classifier=="logit"){
  class_vec<-as.numeric(as.factor(class_vec))-1
  
  testclasslabels<-as.numeric(as.factor(testclasslabels[,1]))
  }
  
  d3<-cbind(class_vec,d2)
  d3<-as.data.frame(d3)
  
  #featlist<-featlist+1
  
  colnames(d2)<-as.character(mz_names)
  colnames(testset)<-as.character(mz_names)
  mz_names<-c("Class",mz_names)
  
  
  
  colnames(d3)<-as.character(mz_names)
  
  
  
  
  d3<-as.data.frame(d3)
  
 # save(d3,file="d3.Rda")
  
  #featlist<-unique(featlist)
  
  mod.lab<-{}
  
  
  if(is.na(col_lab)==TRUE){
    if(length(featlist)>1){
      col_lab<-palette(rainbow(length(featlist)))
    }else{
      col_lab<-c("blue")
    }
  }
  
  extra_index<-which(featlist>(dim(d1)[1]+1))
  
  if(length(extra_index)>0){
    featlist<-featlist[-extra_index]
  }
  featlist<-unique(featlist)
  
  if(is.na(mainlabel)==TRUE){
    
    
    if(classifier=="logitreg" || classifier=="logit"){
      mainlab<-"ROC curves using logistic regression"
    }else{
      
      mainlab<-"ROC curves using SVM"
    }
    
  }else{
    
    mainlab=mainlabel
  }
  

  roc_res_list<-new("list")
  
  if(featincrement==TRUE){
    
    
    featlist<-featlist+1
    
    alltestset<-testset
    
    for(n in 1:length(featlist)){
    
     # roc_res_list<-lapply(1:length(featlist),function(n){  
      testset<-alltestset
      num_select<-featlist[n]
      
      
      d4<-as.data.frame(d3[,c(1:num_select)])
      
      cnames1<-colnames(d4)
      
      #  ##save(d4,file="d4.Rda")
      testset<-as.data.frame(testset)
      
      ###save(testset,file="testset.Rda")
      
      testset<-testset[,which(colnames(testset)%in%colnames(d4))]
      testset<-as.data.frame(testset)
      
      #testset<-setnames(testset, cnames1[-c(1)])
      colnames(testset)<-cnames1[-c(1)]
      
      
      match_names_check<-match(colnames(testset),cnames1[-c(1)])
      
      if(length(which(is.na(match_names_check))==TRUE)>0){
        
        stop("Column names don't match between training and test sets.")
      }
      
      if(max(abs(diff(match_names_check)))>1){
        
        stop("Column order doesn't match between training and test sets.")
      }
      
      
      if(classifier=="logitreg" || classifier=="logit"){
        
                if(ncol(d4)>2){
                     d4[,-c(1)]<-apply(d4[,-c(1)],2,function(x){as.numeric(as.character(x))})
                }else{
                  
                  d4[,c(2)]<-as.numeric(as.character(d4[,c(2)]))
                }
                model1 <- glm(d4$Class~., data=d4, family="binomial")
                
                
                pred<-predict(model1,testset,type="response")
               # pred<-ifelse(pred>,)
                #save(d4,model1,testset,testclasslabels,pred,file="logit.Rda")
                pred1 <- ROCR::prediction(pred, testclasslabels)
        
                 #pred1<-multiclass.roc(testclasslabels,as.numeric(pred)) #,levels=levels(as.factor(y1)))
        
      }else{
        
                #Use SVM
                #save(d4,d3,testset,testclasslabels,svm.type,featlist,d1,d2,kname,num_select,n,file="d4.Rda")
                #model1 <- svm(as.factor(d4$Class)~., data=d4, type=svm.type,probability=TRUE,kernel=kname)
                
                #d4[,-c(1)]<-apply(d4[,-c(1)],2,function(x){as.numeric(as.character(x))})
                
                if(ncol(d4)>2){
                  d4[,-c(1)]<-apply(d4[,-c(1)],2,function(x){as.numeric(as.character(x))})
                }else{
                  
                  d4[,c(2)]<-as.numeric(as.character(d4[,c(2)]))
                }
                
                model1 <- svm(as.factor(d4$Class)~., data=d4, type="C",probability=TRUE,kernel=kname)
                
                
                
                
                
                if(length(levels(as.factor(d4$Class)))>2){
                  pred<-predict(model1,testset,probability=FALSE,decision.values=FALSE,type="prob")
                  
                  #pred<-as.numeric(as.character(pred))
                  
                  #  testclasslabels<-as.numeric(as.character(testclasslabels))
                  
                  data_1<-cbind(testclasslabels,pred)
                  colnames(data_1)<-c("response","predictor")
                  
                  # ##save(data_1,file="data_1.Rda")
                  
                  pred1<-multiclass.roc(testclasslabels,pred,levels=levels(as.factor(d4$Class)),data=data_1)
                  #   ##save(pred1,file="pred1.Rda")
                  
                }else{
                  
                  pred<-predict(model1,testset,probability=TRUE,decision.values=TRUE,type="prob")
                  
                  pred1 <- ROCR::prediction(attributes(pred)$probabilities[,2], testclasslabels)
                  #pred1<-multiclass.roc(factor(testclasslabels),attributes(pred)$probabilities[,2],levels=levels(as.factor(d4$Class)))
                }
                
        
      }
      
      
      if(length(levels(as.factor(d4$Class)))==2){
        
        stats1a <- performance(pred1, 'tpr', 'fpr')
        
        roc_res<-cbind(stats1a@x.values[[1]],stats1a@y.values[[1]])
        colnames(roc_res)<-c(stats1a@x.name,stats1a@y.name)
        fname1=paste("ROC",classifier,".rda",sep="")
        ##saveroc_res,file=fname1)
        x1<-seq(0,1,0.01)
        y1<-x1
        p1<-performance(pred1,"auc")
        mod.lab <-c(mod.lab,paste('using top ',(num_select-1),' m/z features: AUC ',round(p1@y.values[[1]],2),sep=""))
        if(n==1){
          plot(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', ylab=stats1a@y.name, xlab=stats1a@x.name, col=col_lab[n], lty=2, main=mainlab,cex.main=0.7,lwd=2)
        }else{
          lines(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', col=col_lab[n], lty=2,lwd=2)
        }
        
      }
      
      roc_res_list[[n]]<-roc_res
     # return(roc_res)
    } #)
    
    if(length(levels(as.factor(d4$Class)))==2){
      lines(x1, y1, col="black", type="l",lwd=2)
      if(legend==TRUE){
        
        legend('bottomright', c(mod.lab), col=col_lab, lwd=c(2,1), lty=1:3, cex=0.8, bty='n')
      }
    }
  }else{
    
    
    
    
    
    d4<-as.data.frame(d3)
    testset<-as.data.frame(testset)
    
    #  ##saved4,testset,testclasslabels,file="d4.Rda")
    
    testset<-testset[,which(colnames(testset)%in%colnames(d4))]
    
    match_names_check<-match(colnames(testset),colnames(d4))
    
    
    if(length(which(is.na(match_names_check))==TRUE)>0){
      
      stop("Column names don't match between training and test sets.")
    }
    
    if(max(abs(diff(match_names_check)))>1){
      
      stop("Column order doesn't match between training and test sets.")
    }
    
    
    if(classifier=="logit" || classifier=="logitreg"){
      
      
      
      testset<-as.data.frame(testset)
      
      if(ncol(d4)>2){
        d4[,-c(1)]<-apply(d4[,-c(1)],2,function(x){as.numeric(as.character(x))})
      }else{
        
        d4[,c(2)]<-as.numeric(as.character(d4[,c(2)]))
      }
      
      model1 <- glm(d4$Class~., data=d4, family=binomial)
      
      #pred1 <- prediction(fitted(model1), d3$Class)
      
      pred<-predict(model1,testset)
      
      # pred1 <- ROCR::prediction(pred, testclasslabels)
      
      pred1<-multiclass.roc(testclasslabels,as.numeric(pred))
      
    }else{
      
      if(ncol(d4)>2){
        d4[,-c(1)]<-apply(d4[,-c(1)],2,function(x){as.numeric(as.character(x))})
      }else{
        
        d4[,c(2)]<-as.numeric(as.character(d4[,c(2)]))
      }
      
      
      model1 <- svm(as.factor(d4$Class)~., data=d4, type=svm.type,probability=TRUE,kernel=kname)
      
      
      
      if(length(levels(as.factor(d4$Class)))>2){
        pred<-predict(model1,testset,probability=FALSE,decision.values=FALSE,type="prob")
        
        pred<-as.numeric(as.character(pred))
        
        #  testclasslabels<-as.numeric(as.character(testclasslabels))
        
        data_1<-cbind(testclasslabels,pred)
        colnames(data_1)<-c("response","predictor")
        
        # ##save(data_1,file="data_1.Rda")
        
        pred1<-multiclass.roc(testclasslabels,pred,levels=levels(as.factor(d4$Class)),data=data_1)
        # ##save(pred1,file="pred1.Rda")
        
      }else{
        
        
        pred<-predict(model1,testset,probability=TRUE,decision.values=TRUE,type="prob")
        pred1 <- ROCR::prediction(attributes(pred)$probabilities[,2], testclasslabels)
        
        # pred1<-multiclass.roc(testclasslabels,attributes(pred)$probabilities[,2])
        
        
      }
      
      
    }
    
    n=1
    
    if(length(levels(as.factor(d4$Class)))==2){
      
      stats1a <- performance(pred1, 'tpr', 'fpr')
      roc_res<-cbind(stats1a@x.values[[1]],stats1a@y.values[[1]])
      
      
      colnames(roc_res)<-c(stats1a@x.name,stats1a@y.name)
      
      fname1=paste("ROC",classifier,".rda",sep="")
      ##saveroc_res,file=fname1)
      x1<-seq(0,1,0.01)
      y1<-x1
      p1<-performance(pred1,"auc")
      mod.lab <-c(mod.lab,paste('using selected ',(num_select),' m/z features: AUC ',round(p1@y.values[[1]],2),sep=""))
      if(n==1){
        plot(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', ylab=stats1a@y.name, xlab=stats1a@x.name, col=col_lab[n], lty=2, 
             main=mainlab,cex.main=0.7,lwd=2)
      }else{
        lines(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', col=col_lab[n], lty=2,lwd=2)
      }
      
      roc_res_list[[1]]<-roc_res
    }
    
  }
  
#if(FALSE){
  if(length(levels(as.factor(d4$Class)))==2){
    lines(x1, y1, col="black", type="l",lwd=2)
    if(legend==TRUE){
      legend('bottomright', c(mod.lab), col=col_lab, lwd=c(2,1), lty=1:3, cex=0.8, bty='n')
      
    }
    
    auc_res<-p1@y.values[[1]]
  }
  #else{
    
   # auc_res<-pred1$auc
    
    #roc_res<-pred1
    
  #}

  if(newdevice==TRUE){
    try(dev.off(),silent=TRUE)
  }else{
    
    
  }
  options(warn=0)
  return(list("roc"=roc_res)) #prediction"=pred1, #"auc"=auc_res,
}

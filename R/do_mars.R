do_mars <-
function(X,classlabels, analysismode="classification",kfold=10){
  
  dataA<-t(X)
  dataA<-as.data.frame(dataA)
  rnames1<-rownames(X)
  
  if(analysismode=="classification"){
    
    dataA<-cbind(classlabels,dataA)
    dataA<-as.data.frame(dataA)
    
    cnames<-c("Targetvar",paste("mz",seq(1,dim(X)[1]),sep=""))
    
    colnames(dataA)<-as.character(cnames)
    
    if(dim(dataA)[1]<20){kfold=1}
    
    dataA<-as.data.frame(dataA)
    
    attach(dataA,warn.conflicts=FALSE)
    
    #print(dim(dataA))
    
    
    mars_res<-new("list")
    
    
    mars_res[[1]]<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=1,ncross=kfold,nfold=kfold)
    mars_res[[2]]<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=2,ncross=kfold,nfold=kfold)
    mars_res[[3]]<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=3,ncross=kfold,nfold=kfold)
    
    gcv_list<-c(mars_res[[1]]$gcv,mars_res[[2]]$gcv,mars_res[[3]]$gcv)
    
    max_gcv<-which(gcv_list==max(gcv_list,na.rm=TRUE))
    
    mars_res<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=max_gcv[1],ncross=kfold,nfold=kfold)
    
  }else{
    if(analysismode=="regression"){
      
      mars_res<-new("list")
      
      classlabels<-as.data.frame(classlabels)
      #print(dim(classlabels))
      
      if(dim(classlabels)[2]>1){
        classlabels<-apply(classlabels,2,as.numeric)
        classlabels<-as.numeric(classlabels[,1])
      }else{
        classlabels<-as.numeric(classlabels[,1])
        #classlabels<-apply(classlabels,2,as.numeric)
        
      }
      
      
      dataA<-cbind(classlabels,dataA)
      #dataA<-as.data.frame(dataA)
      
      cnames<-c("Targetvar",paste("mz",seq(1,dim(X)[1]),sep=""))
      
      #cnames<-c(paste("mz",seq(1,dim(X)[1]),sep=""))
      
      colnames(dataA)<-as.character(cnames)
      
      if(dim(dataA)[1]<20){kfold=1}
      
      dataA<-as.data.frame(dataA)
      
      attach(dataA,warn.conflicts=FALSE)
      
     # save(dataA,kfold,file="dataA.rda")
      #mars_res[[1]]<-earth(y=(classlabels),x=dataA,Use.beta.cache=FALSE,degree=1,ncross=kfold,nfold=kfold)
      #		mars_res[[2]]<-earth(y=(classlabels),x=dataA,Use.beta.cache=FALSE,degree=2,ncross=kfold,nfold=kfold)
      #		mars_res[[3]]<-earth(y=(classlabels),x=dataA,Use.beta.cache=FALSE,degree=3,ncross=kfold,nfold=kfold)
      
      mars_res[[1]]<-earth(formula=as.numeric(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=1,ncross=kfold,nfold=kfold)
      mars_res[[2]]<-earth(formula=as.numeric(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=2,ncross=kfold,nfold=kfold)
      mars_res[[3]]<-earth(formula=as.numeric(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=3,ncross=kfold,nfold=kfold)
      gcv_list<-c(mars_res[[1]]$gcv,mars_res[[2]]$gcv,mars_res[[3]]$gcv)
      max_gcv<-which(gcv_list==max(gcv_list,na.rm=TRUE))
      
      mars_res<-earth(formula=as.numeric(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=max_gcv[1],ncross=kfold,nfold=kfold)
      #mars_res<-earth(formula=(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=3,ncross=kfold,nfold=kfold)
      #mars_res<-earth(y=(classlabels),x=dataA,Use.beta.cache=FALSE,degree=min_gcv[1],ncross=kfold,nfold=kfold)
      
      
    }else{
      stop("Invalid analysismode entered. Please use classification or regression.")
    }
  }
  dataA<-NULL
  detach(dataA)
  rm(dataA)
  
  ###savemars_res,file="mars_res.Rda")
  #varimp_res <- evimp(mars_res)
  
  #if(FALSE)
  {
    varimp_res <- evimp(mars_res,trim=FALSE)

    #save(varimp_res,file="varimp_res.Rda")    
    #varimp_res<-as.data.frame(varimp_res)
    varimp_res<-varimp_res[order(varimp_res[,1]),]
    rownames(varimp_res)<-rnames1
    
    #print(varimp_res[1:10,])
    varimp_marsres1<-varimp_res #[order(varimp_res[,2],decreasing=TRUE),]
    
    mars_mznames<-rownames(varimp_marsres1)
    
    
    g1<-grep(pattern="NA",x=mars_mznames)
    if(length(g1)>0){
      varimp_marsres1<-varimp_marsres1[-g1,]
    }
  }
  
  
  return(list("mars_model"=mars_res,"mars_varimp"=varimp_marsres1))
  
}

do_plsda_rand <-
function(X,Y,oscmode="pls",numcomp=3,kfold=10,evalmethod="CV",keepX=15,sparseselect=FALSE,analysismode="classification",
                        vip.thresh=1,sample.col.opt="default",sample.col.vec=c("red","green","blue","purple"),scoreplot_legend=TRUE,
                        feat_names=NA,pairedanalysis=FALSE,optselect=FALSE,class_labels_levels_main=NA,
                        legendlocation="bottomleft",plotindiv=TRUE,alphabetical.order=FALSE)
{
  repeatmeasures=pairedanalysis
  
  
  num_var<-dim(X)[1]
  
  if(keepX>num_var){
    keepX=num_var
  }
  
  X<-t(X)
  Y<-as.data.frame(Y)
  
  classlabels<-Y
  
  if(sample.col.opt=="default"){
    
    col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
               "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
               "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
               "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
    
  }else{
    if(sample.col.opt=="topo"){
      
      col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
    }else{
      if(sample.col.opt=="heat"){
        
        col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
      }else{
        if(sample.col.opt=="rainbow"){
          
          col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
          
          
        }else{
          
          if(sample.col.opt=="terrain"){
            
            
            col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
          }
          
          
        }
        
      }
      
    }
  }

  if(FALSE){  
  if(pairedanalysis==FALSE){
    Yclass<-Y[,1]
    Y<-as.numeric(Y[,1])
  }else{
    
    if(dim(Y)[2]>2){
      if(analysismode=="classification"){
        
        Yclass<-as.factor(Y[,2]):as.factor(Y[,3])
      }else{
        Yclass<-Y[,3] #:Y[,3]
      }
      Y<-as.numeric(Yclass)
    }else{
      Yclass<-Y[,2]
      Y<-as.numeric(Y[,2])
    }
    
    
    
  }
  }
  
  
  if(pairedanalysis==FALSE){
    Yclass<-Y[,1]
    
    if(analysismode=="regression"){
      Y<-as.numeric(Y[,1])
    }else{
      
      if(alphabetical.order==FALSE){
        
        Y[,1]<-factor(Y[,1],levels=unique(Y[,1]))
        Yclass<-Y[,1]
      }
      Y<-as.numeric(as.factor(Y[,1]))
    }
    #Y<-as.factor(Y[,1])
    #Yclass<-as.factor(Y[,1])
  }else{
    
    #repeat measures
    if(dim(Y)[2]>2){
      if(analysismode=="classification"){
        
        if(alphabetical.order==FALSE){
          
          Y[,2]<-factor(Y[,2],levels=unique(Y[,2]))
          Y[,3]<-factor(Y[,3],levels=unique(Y[,3])) 
        }
        
        Yclass<-as.factor(Y[,2]):as.factor(Y[,3])
        
        Y<-as.numeric(as.factor(Yclass))
      }else{
        Yclass<-Y[,3] #:Y[,3]
        
        Y<-as.numeric(Yclass)
      }
      
      
    }else{
      if(analysismode=="classification"){
        
        if(alphabetical.order==FALSE){
          
          Y[,2]<-factor(Y[,2],levels=unique(Y[,2]))
          
        }
      }
      
      Yclass<-Y[,2]
      Y<-as.numeric(Y[,2])
      # Y<-as.factor(Y[,2])
    }
    
  }
  
  
  
  #Y<-as.numeric(as.factor(Y[,1]))
  #Y<-as.vector(Y)
  
  #  print(dim(X))
  #print(dim(Y))
  
  #print("starting")
  if(dim(X)[2]>1){
    if(optselect==TRUE){
      if(analysismode=="classification")
      {
        set.seed(123)
        opt_comp<-pls.lda.cv(Xtrain=X, Ytrain=Yclass,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3, priors=NULL)
        
        #cv_res<-plsda_cv(v=kfold,x=X,y=Y,ncomp=opt_comp,errortype=evalmethod)
        
        #print(paste(kfold," CV evaluation using plsda",sep=""))
        #print(cv_res)
      }else{
        if(analysismode=="regression")
        {
          
          set.seed(123)
          opt_comp<-pls.regression.cv(Xtrain=X, Ytrain=Y,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3)
          
          
        }
        
      }
    }else{
      
      opt_comp<-numcomp
      keep_x_vec<-rep(keepX,opt_comp)
    }
    
  }
  
  if(opt_comp<2){
    opt_comp<-2
  }
  if(oscmode=="o1pls"){
    leukemia.pls <- plsr(Y ~ X, ncomp = opt_comp, validation = "LOO")
    ww <- leukemia.pls$loading.weights[,1]
    pp <- leukemia.pls$loadings[,1]
    w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
    t.ortho <- X %*% w.ortho
    
    p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
    Xcorr <- X - tcrossprod(t.ortho, p.ortho)
    
    if(analysismode=="classification")
    {
      cv_res<-plsda_cv(v=kfold,x=Xcorr,y=Yclass,ncomp=opt_comp,errortype=evalmethod)
      #print(paste(kfold," CV evaluation using o1plsda",sep=""))
      #print(cv_res)
    }
    
    X<-Xcorr
  }
  
  if(oscmode=="o2pls"){
    leukemia.pls <- plsr(Y ~ X, ncomp = opt_comp, validation = "LOO")
    ww <- leukemia.pls$loading.weights[,1]
    pp <- leukemia.pls$loadings[,1]
    w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
    t.ortho <- X %*% w.ortho
    
    p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
    Xcorr <- X - tcrossprod(t.ortho, p.ortho)
    
    
    if(analysismode=="classification")
    {
      cv_res<-plsda_cv(v=kfold,x=Xcorr,y=Yclass,ncomp=opt_comp,errortype=evalmethod)
      
      
     # print(paste(kfold," CV evaluation using o2plsda",sep=""))
      #print(cv_res)
    }
    
    X<-Xcorr
  }
  
  
  
  bad_variables<-{}
  
  if(sparseselect==TRUE)
  {
    
    if(analysismode=="classification")
    {
      if(optselect==TRUE){
        keepx_seq<-seq(5,keepX,5)
        best_cv_res<-c(0)
        best_kvec<-c(5)
        for(kvec in keepx_seq){
          keep_x_vec<-rep(kvec,opt_comp)
          
          if(repeatmeasures==TRUE){
            
            #print("spls classlabels")
            #   print(classlabels)
            #print(dim(X))
            
            #linn.pls <- multilevel(X=X, design=classlabels,ncomp = opt_comp,
            #keepX = keep_x_vec, method = 'splsda')
            
            linn.pls <- try(mixOmics::multilevel(X=X, design=classlabels,ncomp = opt_comp,
                                                 keepX = keep_x_vec, method = 'splsda'),silent=TRUE)
            
            if(is(linn.pls,"try-error")){
              
              
              
              
              linn.pls <- mixOmics::splsda(X=X,Y=classlabels[,-c(1)],ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,1])
              
            }
            
            #print(linn.pls)
            
          }else{
            
            linn.pls <- mixOmics::splsda(X, Yclass,ncomp=opt_comp,keepX=keep_x_vec)
          }
          
          
          
          #linn.vip<-vip(linn.pls)
          linn.vip<-linn.pls$loadings$X
          #linn.vip<-linn.vip[,1]
          
          bad_variables<-linn.pls$nzv$Position
          
          good_feats<-{}
          for(c1 in 1:opt_comp){
            good_feats<-c(good_feats,which(linn.vip[,c1]!=0))
          }
          
          good_feats<-unique(good_feats)
          
          
          
          if(length(good_feats)>1){
            
            cv_res<-plsda_cv(v=kfold,x=X[,good_feats],y=Yclass,ncomp=opt_comp,errortype=evalmethod)
            
            #cv_res<-pls.lda.cv(Xtrain=X, Ytrain=Y,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3, priors=NULL)
            
            #cv_res<-cv_res$cv
        #    print(paste(kfold," CV evaluation using spls top ",kvec," features per component",sep=""))
         #   print(cv_res$mean_acc)
            if(cv_res$mean_acc>best_cv_res){
              
              best_cv_res<-cv_res$mean_acc
              best_kvec<-kvec
              
              #best_ncomp<-cv_res$ncomp
            }
          }else{
            print("Too few variables to perfrom CV.")
          }
        }
        
        keep_x_vec<-rep(best_kvec,opt_comp)
      }else{
        
        keep_x_vec<-rep(keepX,opt_comp)
      }
    }
    else{
      #keep_x_vec<-rep(dim(X)[2],opt_comp)
      keep_x_vec<-rep(keepX,opt_comp)
      
    }
    if(analysismode=="regression"){
      
      if(repeatmeasures==TRUE){
        
        #print("spls classlabels")
        #print(classlabels)
        #print(dim(X))
        #print(keep_x_vec)
        
        
        # linn.pls <- multilevel(X=X, design=classlabels,ncomp = opt_comp,
        #keepX = keep_x_vec, method = 'spls')
        
        linn.pls <- try(mixOmics::multilevel(X=X,Y=Y,design=classlabels[,2],ncomp = opt_comp,
                                             keepX = keep_x_vec, method = 'spls'),silent=TRUE)
        
        if(is(linn.pls,"try-error")){
          
          
          
          
          linn.pls <- mixOmics::spls(X=X,Y=Y,ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,2])
          
        }
        
        
        
        
      }else{
        linn.pls <- mixOmics::spls(X, Y,ncomp=opt_comp,keepX=keep_x_vec,mode="regression")
      }
    }else{
      
      if(repeatmeasures==TRUE){
        #print("spls classlabels")
        #print(classlabels)
        #print(dim(X))
        
        
        #linn.pls <- multilevel(X=X, design=classlabels,ncomp = opt_comp,
        #keepX = keep_x_vec, method = 'splsda')
        
        linn.pls <- try(multilevel(X=X, design=classlabels,ncomp = opt_comp,
                                   keepX = keep_x_vec, method = 'splsda'),silent=TRUE)
        
        if(is(linn.pls,"try-error")){
          
          
          
          
          linn.pls <- mixOmics::splsda(X=X,Y=classlabels[,-c(1)],ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,1])
          
        }
        
        if(FALSE){
          stimu.time<-data.frame(cbind(as.character(classlabels[,2]),
                                       as.character(classlabels[,3])))
          repeat.stimu2<-classlabels[,1]
          
          res.2level <- multilevel(X, cond = stimu.time,
                                   sample = repeat.stimu2, ncomp = 3,
                                   keepX = keep_x_vec, tab.prob.gene = NULL, method = 'splsda')
        }
        
      }else{
        linn.pls <- mixOmics::splsda(X, Y,ncomp=opt_comp,keepX=keep_x_vec)
      }
    }
    
    linn.vip<-linn.pls$loadings$X
    
    
    
    
    
  }else{
    
    # print("opt comp")
    #print(opt_comp)
    if(analysismode=="regression"){
      
      if(repeatmeasures==TRUE){
        linn.pls <- mixOmics::pls(X, Y,ncomp=opt_comp,multilevel=classlabels[,1])
      }else{
        linn.pls <- mixOmics::pls(X, Y,ncomp=opt_comp)
      }
      
    }else{
      
      if(repeatmeasures==TRUE){
        linn.pls <- mixOmics::plsda(X, Yclass,ncomp=opt_comp,multilevel=classlabels[,1])
      }else{
        linn.pls <- mixOmics::plsda(X, Yclass,ncomp=opt_comp)
      }
    }
    
    linn.vip<-mixOmics::vip(linn.pls)
    
    
    
    
    
  }
  
  
  v1<-{}
  cv_res<-{}
  
  
  
  
  
  return(list("model"=linn.pls,"vip_res"=linn.vip))
  
}

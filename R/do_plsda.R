do_plsda <-
function(X,Y,oscmode="pls",numcomp=3,kfold=10,evalmethod="CV",keepX=15,sparseselect=FALSE,analysismode="classification",
                   vip.thresh=1,sample.col.opt="default",sample.col.vec=c("red","green","blue","purple"),
                   scoreplot_legend=TRUE,feat_names=NA,pairedanalysis=FALSE,optselect=FALSE,class_labels_levels_main=NA,
                   legendlocation="bottomleft",plotindiv=TRUE,pls.vip.selection="max",output.device.type="pdf",
                   plots.res=600,plots.width=8,plots.height=8,plots.type="cairo",pls.ellipse=TRUE,alphabetical.order=FALSE)
{
  repeatmeasures=pairedanalysis
  
  #print("Starting here")
  if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/PLS_performance_plots.pdf"
    #pdf(temp_filename_1)
    
   # pdf(temp_filename_1,width=plots.width,height=plots.height)
  }
  
  
  
  num_var<-dim(X)[1]
  
  if(keepX>num_var){
    keepX=num_var
  }
  
  
  
  X<-t(X)
  
  
  Y<-as.data.frame(Y)
  
  classlabels<-Y    
  
  
# save(X,Y,pairedanalysis,classlabels,file="plspaireddebug.Rda")
  
  #only one column for classlabels if mode=classificaion unpaired
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
  
  class_labels_levels<-levels(as.factor(Yclass))
  alphacol=0.3
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
        #col_vec<-heat.colors(256) #length(class_labels_levels))
        
        col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
      }else{
        if(sample.col.opt=="rainbow"){
          #col_vec<-heat.colors(256) #length(class_labels_levels))
          col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
          
          #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
        }else{
          
          if(sample.col.opt=="terrain"){
            #col_vec<-heat.colors(256) #length(class_labels_levels))
            #col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
            
            col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
          }else{
            
            if(sample.col.opt=="colorblind"){
              #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
              # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
              
              if(length(class_labels_levels)<9){
                
                col_vec <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "#E64B35FF", "grey57")
                
              }else{
                
                #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                           "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                
              }
              
              
            }else{
              
              check_brewer<-grep(pattern="brewer",x=sample.col.opt)
              
              if(length(check_brewer)>0){
                sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                
              }else{
                
                if(sample.col.opt=="journal"){
                  
                  col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                             "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                             "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                             
                             "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                             "#E64B3519","#4DBBD519","#631879E5","grey75")
                  if(length(class_labels_levels)<8){
                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75")
                    
                    #col_vec2<-brewer.pal(n = 8, name = "Dark2")
                    
                  }else{
                    if(length(class_labels_levels)<=28){
                      # col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "grey75","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666","#1B9E77", "#7570B3", "#E7298A", "#A6761D", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
                      
                      col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                 "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                 "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                 
                                 "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF", "#8BD76BFF",
                                 "#E64B3519","#9DBBD0FF","#631879E5","#666666","grey75")
                      
                    }else{
                      
                      
                      
                      
                      colfunc <-colorRampPalette(c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75"));col_vec<-colfunc(length(class_labels_levels))
                      
                      col_vec<-col_vec[sample(col_vec)]
                      
                      
                    }
                  }
                }else{
                  #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                #  if(length(sample.col.opt)==1){
                 #   col_vec <-rep(sample.col.opt,length(class_labels_levels))
                #  }else{
                    
                 #   colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                    
                #  }
                  
                  if(length(sample.col.opt)==1){
                    col_vec <-rep(sample.col.opt,length(class_labels_levels))
                  }else{
                    
                    if(length(sample.col.opt)>=length(class_labels_levels)){
                      
                      col_vec <-sample.col.opt
                      col_vec <- rep(col_vec,length(class_labels_levels))
                      
                      
                    }else{
                      colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                    }
                    
                  }
                }
                
              }
              
            }
          }
          
          
        }
        
      }
      
    }
  }
  
  
  
  
  #print("starting")
  if(dim(X)[2]>1){
    if(optselect==TRUE){
      if(analysismode=="classification")
      {
        set.seed(123)
        opt_comp<-plsgenomics::pls.lda.cv(Xtrain=X, Ytrain=Yclass,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3, priors=NULL)
        
        
      }else{
        if(analysismode=="regression")
        {
          
          set.seed(123)
          opt_comp<-plsgenomics::pls.regression.cv(Xtrain=X, Ytrain=Y,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3)
          
          
        }
        
      }
    }else{
      
      opt_comp<-numcomp
      keep_x_vec<-rep(keepX,opt_comp)
    }
    
  }
  
  suppressWarnings(dir.create("Tables",showWarnings = FALSE))
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
      
      
      #print(paste(kfold," CV evaluation using o2plsda",sep=""))
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
            
            
            
            linn.pls <- try(multilevel(X=X, design=classlabels,ncomp = opt_comp,
                                       keepX = keep_x_vec, method = 'splsda'),silent=TRUE)
            
            if(is(linn.pls,"try-error")){
              
              
              
              
              linn.pls <- mixOmics::splsda(X=X,Y=classlabels[,-c(1)],ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,1])
              
            }
            
          }else{
            
            linn.pls <- mixOmics::splsda(X, Yclass,ncomp=opt_comp,keepX=keep_x_vec)
          }
          
          
          
          linn.vip<-linn.pls$loadings$X
          
          
          bad_variables<-linn.pls$nzv$Position
          
          good_feats<-{}
          for(c1 in 1:opt_comp){
            good_feats<-c(good_feats,which(linn.vip[,c1]!=0))
          }
          
          good_feats<-unique(good_feats)
          
          
          
          if(length(good_feats)>1){
            
            cv_res<-try(plsda_cv(v=kfold,x=X[,good_feats],y=Yclass,ncomp=opt_comp,errortype=evalmethod),silent=TRUE)
            
            
           # print(paste(kfold," CV evaluation using spls top ",kvec," features per component",sep=""))
          #  print(cv_res$mean_acc)
            if(cv_res$mean_acc>=best_cv_res){
              
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
      
      keep_x_vec<-rep(keepX,opt_comp)
      
    }
    
    
    if(analysismode=="regression"){
      
      #   #save(X,Y,opt_comp,classlabels,file="splsreg.Rda")
      if(repeatmeasures==TRUE){
        
        #  ###savelist=ls(),file="debugspls.Rda")
        
        
        
        
        linn.pls <- try(multilevel(X=X, design=classlabels[,2],ncomp = opt_comp,
                                   keepX = keep_x_vec, method = 'spls'),silent=TRUE)
        
        
        if(is(linn.pls,"try-error")){
          
          
          
          
          linn.pls <- mixOmics::spls(X=X,Y=Y,ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,2])
          
        }
        
        
        
        
      }else{
        linn.pls <- mixOmics::spls(X, Y,ncomp=opt_comp,keepX=keep_x_vec,mode="regression")
      }
    }else{
      #analysismode classification
      if(repeatmeasures==TRUE){
        
        
        linn.pls <- try(multilevel(X=X, design=classlabels,ncomp = opt_comp,
                                   keepX = keep_x_vec, method = 'splsda'),silent=TRUE)
        
        if(is(linn.pls,"try-error")){
          
          
          linn.pls <- mixOmics::splsda(X=X,Y=classlabels[,-c(1)],ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,1])
          
        }
        
        #  ###savelinn.pls,file="linnpls.Rda")
        
        
        
      }else{
        linn.pls <- mixOmics::splsda(X, Y,ncomp=opt_comp,keepX=keep_x_vec)
      }
    }
    
    linn.vip<-linn.pls$loadings$X
    
    
    bad_variables<-linn.pls$nzv$Position
    
    #VIP based feature selection
    good_feats<-{}
    for(c1 in 1:opt_comp){
      good_feats<-c(good_feats,which(linn.vip[,c1]!=0))
    }
    
    good_feats<-unique(good_feats)
    
    
    
  }else{
    #PLS
    if(analysismode=="regression"){
      
      
      if(repeatmeasures==TRUE){
        
        ##save(X,Y,opt_comp,classlabels,file="plsreg.Rda")
        linn.pls <- mixOmics::pls(X=X,Y=Y,ncomp=opt_comp,multilevel=classlabels[,1])
        
        #  print("done here")
        
      }else{
        
        linn.pls <- mixOmics::pls(X, Y,ncomp=opt_comp)
      }
      
      
    }else{
      
      if(analysismode=="classification"){
        if(repeatmeasures==TRUE){
          
          # print("multilevel PLS classlabels")
          #print(head(classlabels))
          
          linn.pls <- mixOmics::plsda(X, Yclass,ncomp=opt_comp,multilevel=classlabels[,1])
        }else{
          linn.pls <- mixOmics::plsda(X, Yclass,ncomp=opt_comp)
        }
        
      }
    }
    
    
    linn.vip<-mixOmics::vip(linn.pls)
    
    
    bad_variables<-linn.pls$nzv$Position
    
    
    
    good_feats<-{}
    c1<-1
    
    
    good_feats<-{}
    
    if(pls.vip.selection=="max"){
      if(opt_comp>1){
        for(c1 in 1:opt_comp){
          good_feats<-c(good_feats,which(linn.vip[,c1]>vip.thresh))
        }
        
      }else{
        good_feats<-which(linn.vip>vip.thresh)
        
      }
      
    }else{
      
      if(opt_comp>1){
        linn.vip.mean<-apply(linn.vip,1,mean)
        
        good_feats<-which(linn.vip.mean>vip.thresh)
      }else{
        good_feats<-which(linn.vip>vip.thresh)
        
      }
      
    }
    good_feats<-unique(good_feats)
    
  }
  
  
  v1<-{}
  cv_res<-{}
  if(length(good_feats)>1)
  {
    
    #print(paste(oscmode," PLS evaluation using selected variables",sep=""))
    
    
    if(analysismode=="classification"){
      
      if(FALSE){
        if(length(Y)>10){ 
          cv_res<-try(plsda_cv(v=kfold,x=X,y=Y,ncomp=opt_comp,errortype=evalmethod),silent=TRUE)
          
          print(paste(kfold," CV ", evalmethod, " using all features: ",cv_res,sep=""))
          
          
          cv_res<-try(plsda_cv(v=kfold,x=X[,good_feats],y=Y,ncomp=opt_comp,errortype=evalmethod),silent=TRUE)
          
          print(paste(kfold," CV ", evalmethod, " using top features: ",cv_res,sep=""))
        }
        
        if(length(Y)>30){
        set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
        perf.plsda <- perf(linn.pls, validation = "Mfold", folds = 5, 
                           progressBar = FALSE, auc = TRUE, nrepeat = 10) 
        }else{
          set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
          perf.plsda <- perf(linn.pls, validation = "loo", 
                             progressBar = FALSE, auc = TRUE) 
          
        }
        # perf.plsda.srbct$error.rate  # error rates
        plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")
        
        
      }
    }else{
      
      
      print(paste(oscmode," PLS evaluation using all variables",sep=""))
      
      temp_d<-cbind(Y,X[,good_feats])
      temp_d<-as.data.frame(temp_d)
      
     # linn.pls3 <- mixOmics::pls(X[,good_feats], Y)
      
      linn.pls2 <- mixOmics::pls(X, Y,ncomp=opt_comp)
      r2_q2_valid_res<-{}
      
      if(length(Y)>30){
        
        if(opt_comp>1){
          
          #linn.pls2 <- pls(X, Y,ncomp=opt_comp) #pls(X, Y,ncomp=opt_comp)
          #v1<-try(perf(linn.pls2,validation="loo"),silent=TRUE)
          
          v1<-try(mixOmics::perf(linn.pls2,validation="Mfold",folds=kfold),silent=TRUE)
          
          
          
          if(is(v1,"try-error")){
            
          }else{
            
            r2_q2_valid_res<-rbind(v1$R2,v1$Q2,v1$MSEP)
            
            if(nrow(r2_q2_valid_res)>0){
              rownames(r2_q2_valid_res)<-c("R2","Q2","MSEP")
            }
            
            cnames_vres<-paste("PLScomp",seq(1,opt_comp),sep="")
            colnames(r2_q2_valid_res)<-cnames_vres
            
            
            
            write.table(r2_q2_valid_res,file="Tables/pls_r2_q2_res_usingallfeatures.txt",sep="\t",row.names=TRUE)
            if(plotindiv==TRUE){
              w <- 0.1 #grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS leave-one-out validation diagnostics using all features",ylab="Variation",col=c("darkgrey","lightgrey"))
              #legend("topright",c("R2","Q2"),col=c("darkgrey","lightgrey"),pch=c(20))
              
              le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("R2","Q2"), col=c("darkgrey","lightgrey"),pch = c(19), pt.cex = 0.6, title = "",cex=0.8))
              
            }
            
            #linn.pls2 <- pls(X, Y,ncomp=opt_comp) #pls(X, Y,ncomp=opt_comp)
            
          }
          
        }
        
      }else{
        r2_q2_valid_res<-{}
        if(opt_comp>1){
          print("PLS loo validation diagnostics using all features")
          
          v1<-try(mixOmics::perf(linn.pls2,validation="loo"),silent=TRUE)
          
          if(is(v1,"try-error")){
            
          }else{
            
            r2_q2_valid_res<-rbind(v1$R2,v1$Q2,v1$MSEP)
            
            if(nrow(r2_q2_valid_res)>0){
              rownames(r2_q2_valid_res)<-c("R2","Q2","MSEP")
            }
            
            cnames_vres<-paste("PLScomp",seq(1,opt_comp),sep="")
            colnames(r2_q2_valid_res)<-cnames_vres
            
            if(plotindiv==TRUE){
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS loo validation diagnostics \n using all features",ylab="Variation",col=c("darkgrey","lightgrey"))
              
              le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("R2","Q2"), col=c("darkgrey","lightgrey"),pch = c(19), pt.cex = 0.6, title = "",cex=0.8))
              
            }
          }
          v1<-try(perf(linn.pls3,validation="loo"),silent=TRUE)
          if(is(v1,"try-error")){
          }else{
            
            print("PLS loo validation diagnostics using selected features")
            
            
            r2_q2_valid_res<-rbind(v1$R2,v1$Q2,v1$MSEP)
            
            if(nrow(r2_q2_valid_res)>0){
              rownames(r2_q2_valid_res)<-c("R2","Q2","MSEP")
            }
            
            cnames_vres<-paste("PLScomp",seq(1,dim(r2_q2_valid_res)[2]),sep="")
            colnames(r2_q2_valid_res)<-cnames_vres
            
            write.table(r2_q2_valid_res,file="Tables/pls_r2_q2_res_usingselectfeats.txt",sep="\t",row.names=TRUE)
            
            if(plotindiv==TRUE){
              w <- 0.1 #grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS loo validation diagnostics \n using selected features",ylab="Variation",col=c("darkgrey","lightgrey"))
              # legend("topright",c("R2","Q2"),col=c("darkgrey","lightgrey"),pch=c(20))
              
              le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("R2","Q2"), col=c("darkgrey","lightgrey"),pch = c(19), pt.cex = 0.6, title = "",cex=0.8))
            }
            
          }
          
        }
        
      }
      
    }
    
    
  }else{
    print("No variables selected.")
  }
  
  
  # print("Done with plsda")
  
  #linn.pls2 <- pls(X, Y,ncomp=opt_comp)
  
  SS<-get_plscompvar(linn.pls,nvar=dim(X)[2],opt_comp)
  
  # SS<-c(linn.pls$loadings$X)^2*colSums(linn.pls$variates$X^2)
  
  pls_var<-100*SS/(sum(SS))
  
  pls_var<-round(pls_var,2)
  
  
  
  if(output.device.type!="pdf"){
    try(dev.off(),silent=TRUE)
    
  }
  
  
  
  #barplot(pls_var,main="PLS %variation per component",cex.main=0.7)
  
  if(analysismode=="classification")
  {
    # color for plotIndiv
    col.stimu = as.numeric(Y)
    
    #print("plotting PLS")
    #print(opt_comp)
    class_labels_levels<-levels(as.factor(Yclass))
    color_vec<-col_vec #sample.col.vec #rainbow(length(class_labels_levels), start = 0, end = 0.1) #c("green","purple")
    col.stimu<-color_vec[col.stimu]
    
    class_labels_levels2<-class_labels_levels
    # pch for plots
    pch.time = rep(15, length(class_labels_levels))
    #pch.time[time == 't2'] = 4
    
    pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21)
    
    Yclass2=Yclass
    
    
    samplelabels<-as.data.frame(Yclass)
    samplelabels<-as.factor(samplelabels[,1])
    l2<-levels(as.factor(samplelabels))
    col_all=topo.colors(256)
    
    t1<-table(samplelabels)
    if(is.na(class_labels_levels)==TRUE){
      
      l1<-levels(as.factor(samplelabels))
    }else{
      l1<-class_labels_levels
      
      
    }
    
    class_labels_levels<-l1
    
    col <- rep(col_vec[1:length(t1)], t1)
    #col<-rep(col_all[1:length(l1)],t1)
    ## Choose different size of points
    cex <- rep(2, length(Yclass))
    
    pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
    pch <- rep(15,length(Yclass))
    cex <- rep(2, length(Yclass))
    for(p1 in 1:length(l2)){
      
      pch[which(samplelabels==l2[p1])]=pch_vec[p1]
    }
    
    
    if(pairedanalysis==TRUE){
      
      if(ncol(classlabels)>2){
        
        class_labels_levels2<-levels(as.factor(classlabels[,2]):as.factor(classlabels[,3]))
        Yclass2<-as.factor(classlabels[,2]):as.factor(classlabels[,3])
      }else{
        class_labels_levels2<-levels(as.factor(classlabels[,2]))
        Yclass2=classlabels[,2]
      }
    }
    
    #  col.stimu = as.numeric(Yclass2)
    #color_vec<-col_vec #sample.col.vec #rainbow(length(class_labels_levels), start = 0, end = 0.1) #c("green","purple")
    col.stimu<-col #color_vec[col.stimu]
    
    #pch_vec<-seq(1,length(Yclass2))
    
    for(p1 in 1:length(class_labels_levels2)){
      
      pch.time[which(Yclass2==class_labels_levels2[p1])]=pch_vec[p1]
    }
    pch.time=pch
    
    # ###savelist=ls(),file="debug.Rda")
    if(plotindiv==TRUE){
      
      
      if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/PLS_pairwise_component_plots.pdf"
        
        #pdf(temp_filename_1)
        pdf(temp_filename_1,width=plots.width,height=plots.height)
        #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      }
      
   #  save(X,linn.pls,pls_var,Yclass,opt_comp,col.stimu,sample.col.opt,class_labels_levels,pls.ellipse,file="plsplots.Rda")
     
    # print(Sys.time())
      get_plsplots(X,plsres=linn.pls,plsvar=pls_var,samplelabels=Yclass,filename=NA,ncomp=opt_comp,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=col.stimu,
                   sample.col.opt=sample.col.opt,alphacol=0.3,legendlocation="topright",class_levels=class_labels_levels,pls.ellipse=pls.ellipse,alphabetical.order=alphabetical.order)
      #,silent=TRUE)
     # print(Sys.time())
      if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
      }
      
      
      
    }
    #legend = c(class_labels_levels), cex = 0.55)
  }
  
 # print("Done with PLSDA")
  write.table(linn.pls$variates$X,file="Tables/pls_scores.txt",sep="\t")
  write.table(linn.pls$loadings$X,file="Tables/pls_loadings.txt",sep="\t")
  ####savelinn.pls,file="pls_res.Rda")
  return(list("model"=linn.pls,"vip_res"=linn.vip,"valid_res"=v1,"cv_res"=cv_res,"opt_comp"=opt_comp,"selected_variables"=good_feats,"bad_variables"=bad_variables))
  
}

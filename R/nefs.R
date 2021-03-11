nefs <-
function(X=NA,Y=NA,feature_table_file=NA,class_labels_file=NA,feat.sel.methods=c("rf","rfe","pls","limma","lasso"),num.var.sel=10,prop.select.thresh=0.7,split.train.test=FALSE,train.pct=0.7,outloc,kfold=10,pca.ellipse=TRUE,Xtest=NA,Ytest=NA,rsdthresh=1,pls_vip_thresh=2,seedvalue=27611,learningsetmethod="CV",confounder.matrix=NA,fdrmethod="BH",fdrthresh=0.05,num.methods.sel=1,globalcor=FALSE,cor.method="spearman",networktype="complete",abs.cor.thresh=0.4,cor.fdrthresh=0.05,max.cor.num=100,net_node_colors=c("green","red"), net_legend=TRUE,niter=10,output.device.type="pdf",heatmap.col.opt="RdBu",boxplot.col.opt=c("white"),barplot.col.opt=c("grey57","grey90"),sample.col.opt="rainbow",mz.thresh=1,time.thresh=10,svm_kernel="radial",good_feats_index=NA,pvalue.thresh=0.05,plots.width=8,plots.height=8,plots.res=600, plots.type="cairo",num_nodes=2,ylabel="Intensity",cex.plots=0.7,tune_classifiers=FALSE,find.common.features=FALSE,aggregation.method="consensus",aggregation.max.iter=1000,add.pvalues=TRUE,add.jitter=TRUE,alphabetical.order=FALSE,boxplot.type="ggplot",balance.classes=FALSE,deeplearning=FALSE)
{
  
  options(warn=-1)
  
  suppressMessages(library(CMA))
  suppressMessages(library(h2o))
  suppressMessages(library(RankAggreg))
  suppressMessages(library(Boruta))
  suppressMessages(library(stepPlr))
  suppressMessages(library(glmnet))
  
  match_class_dist=TRUE
  analysistype="oneway"
  iter.quantile.thresh=prop.select.thresh
  #library(randomForest)
  #library(CMA)
  # importance=randomForest::importance
  #c("rfe","rf","limma","lasso","elasticnet","wilcox.test","pls","spls","ospls","opls","f.test","t.test")
  
  confounderfdrmethod=fdrmethod
  confounderfdrthresh=fdrthresh
  
  
  try(unlockBinding("importance", as.environment("package:ranger")),silent=TRUE)
  #assign("importance", importance, "package:randomForest")
  try(assign("importance", randomForest::importance, "package:ranger"),silent=TRUE)
  errortype="BER"
  tune_scda<-NA
  if(typeof(X)=="logical"){
    X<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  }
  
  if(typeof(Y)=="logical"){
    Y<-read.table(class_labels_file,sep="\t",header=TRUE)
  }
  
  
  X<-as.data.frame(X)
  
  
  if(is.na(Xtest)==TRUE){
    
    test_data_check<-NA
  }else{
    test_data_check<-1
  }
  
  if(FALSE){
    cl<-makeCluster(num_nodes)
    
    clusterExport(cl,"do_rsd")
    
    feat_rsds<-parApply(cl,X[,-c(1:2)],1,do_rsd)
    
    stopCluster(cl)
    
    abs_feat_rsds<-abs(feat_rsds)
    
    good_metabs<-which(abs_feat_rsds>rsdthresh)
    if(length(good_metabs)>0){
      
      X<-X[good_metabs,]
    }else{
      
      stop("No features meet the defined rsd threshold.")
    }
    
  }
  
  Xorig<-X
  cnames1<-colnames(X)
  
  cnames1<-tolower(cnames1)
  
  check_name1<-grep(cnames1,pattern="name|Name")
  if(length(check_name1)>0){
    
    mzrt<-X$Name
    sampnames<-colnames(X[,-c(1)])
    mzrt_id<-X$Name
    #Transpose X and Xtrain; rows are samples; columns are metabs
    X<-t(X[,-c(1)])
    
  }else{
    mzrt<-X[,c(1:2)]
    
    sampnames<-colnames(X[,-c(1:2)])
    
    mzrt_id<-paste(X$mz,"_",X$time,sep="")
    
    #Transpose X and Xtrain; rows are samples; columns are metabs
    X<-t(X[,-c(1:2)])
  }
  
  
  if(num.var.sel>dim(Xorig)[1]){
    
    num.var.sel<-dim(Xorig)[1]
  }
  
  
  
  suppressWarnings(dir.create(outloc,showWarnings = FALSE))
  setwd(outloc)
  
  outloc=getwd()
  
  if(output.device.type!="pdf"){
    dir.create("Figures",showWarnings = FALSE)
  }
  fname<-"InputParameters.csv"
  
  
  nsamp<-dim(X)[1]
  
  Xtrain<-X
  Ytrain_mat<-Y
  Ytest_mat<-NA
  
  if(split.train.test==TRUE){
    if(nsamp<30){
      print("N is too small for train/test analysis. Continuing using all data as training set.")
      split.train.test=FALSE
      Xtrain<-X
      Y<-as.factor(Y[,2])
      Ytrain_mat<-Y
      Ytest_mat<-NA
    }else{
      
      
      test_data_check=1
      numtrain<-round(train.pct*nsamp)
      suppressMessages(library(CMA))
      set.seed(seedvalue)
      train_test_sets<-GenerateLearningsets(y=Y[,2],method="MCCV",ntrain=numtrain,niter=1,strat=TRUE,fold=kfold)
      allindex<-1:nsamp
      set.seed(seedvalue)
      train_index<-train_test_sets@learnmatrix[1,]
      check_index<-which(train_index==0)
      if(length(check_index)>0){
        train_index<-train_index[-check_index]
      }
      test_index<-allindex[-train_index]
      Ytrain_mat<-Y[train_index,]
      Ytest_mat<-Y[test_index,]
      Xtest<-X[test_index,]
      Ytest<-as.factor(Y[test_index,2])
      Y<-as.factor(Y[train_index,2])
      Ytrain<-Y
      Xtrain<-X[train_index,]
      print("Dim of train set")
      print(dim(Xtrain))
      print("Dim of test set")
      print(dim(Xtest))
      
      if(balance.classes==TRUE){
        
        data1=cbind(Ytrain,Xtrain)
        data1<-as.data.frame(data1)
        
        # colnames(data1)<-c("Ytrain",mzrt)
        ###save(data1,file="data.Rda")
        
        data1$Ytrain<-as.factor(data1$Ytrain)
        
        ## now using SMOTE to create a more "balanced problem"
        #newData <- SMOTE(Ytrain ~ ., data1, perc.over = 600,perc.under=100)
        
        colnames(data1)<-c("Ytrain",paste("var",seq(1,ncol(data1)-1),sep=""))
        newData <- ROSE(Ytrain ~ ., data1, seed = 1,N=nrow(data1)*2)$data
        print(table(newData$Ytrain))
        
        Xtrain<-newData[,-c(1)]
        Xtrain<-as.matrix(Xtrain)
        Ytrain<-newData[,c(1)]
        Y<-Ytrain
        Ytrain_mat<-cbind((rownames(Xtrain)),(Ytrain))
        Ytrain_mat<-as.data.frame(Ytrain_mat)
        print("new data")
        print(dim(Xtrain))
        print(dim(Ytrain_mat))
        print(length(Ytrain))
      }
      
      
      
      
      if(length(check_name1)>0){
        
        sampnames<-colnames(Xorig[,c(train_index+1)])
        sampnames_test<-colnames(Xorig[,c(test_index+1)])
      }else{
        sampnames<-colnames(Xorig[,c(train_index+2)])
        sampnames_test<-colnames(Xorig[,c(test_index+2)])
      }
      
      if(is.na(confounder.matrix)==FALSE){
        confounder.matrix<-confounder.matrix[train_index,]
      }
    }
    
  }else{
    
    Xtrain<-X
    Y<-as.factor(Y[,2])
    
    if(is.na(Xtest)==FALSE){
      
      Xtestorig<-Xtest
      
      mzrt_test<-Xtest[,c(1:2)]
      
      
      sampnames_test<-colnames(Xtest[,-c(1:2)])
      Xtest<-t(Xtest[,-c(1:2)])
      Ytest<-as.factor(Ytest[,2])
      
      #if(find.common.features==TRUE)
      {
        res1<-getVenn(dataA=mzrt,dataB=mzrt_test,name_a="train",name_b="test",time.thresh=time.thresh,mz.thresh=mz.thresh,xMSanalyzer.outloc=getwd())
        
        if(nrow(res1$common)<1){
          stop("No common features found.")
        }else{
          
          if(nrow(res1$common)!=nrow(mzrt)){
            #    print("Not all features were common between the train and test sets. Only using the common features for further analysis.")
            
          }
          
          #  print("Number of common features:")
          #print(nrow(res1$common))
          
          #print(head(res1$common))
          
          Xtrain<-Xtrain[,unique(res1$common$index.A)]
          
          #matching_train_data<-matching_train_data[order(matching_train_data$mz),]
          
          Xtest<-Xtest[,unique(res1$common$index.B)]
          
        }
        
      }
      
      
    }else{
      
      print("Using the training set as the test set")
      Xtest=Xtrain
      Ytest=Y
      sampnames_test<-colnames(Xorig[,-c(1:2)])
    }
    
  }
  
  X<-(Xtrain)
  
  svm_acc<-{}
  class_levels<-levels(as.factor(Y))
  
  
  if(learningsetmethod=="bootstrap"){
    
    niter=1000
  }else{
    niter=niter
  }
  
  if(is.na(Xtest)==FALSE)
  {
    #save(Xtrain,file="Xtrain.Rda")
    #save(Xtest,file="Xtest.Rda")
    #save(Y,file="Y.Rda")
    #save(Ytest,file="Ytest.Rda")
    
    if(ncol(Xtrain)!=ncol(Xtest)){
      
      stop("The train and test sets should have same number of variables.")
      
    }
  }
  
  suppressMessages(library(CMA))
  importance<-randomForest::importance
  set.seed(seedvalue)
  fiveCV10iter<-GenerateLearningsets(y=Y,method=learningsetmethod,fold=kfold,niter=niter,strat=TRUE)
  
  ###save(fiveCV10iter,file="fiveCV10iter.Rda")
  
  feat_sel_matrix<-matrix(nrow=dim(X)[2],ncol=length(feat.sel.methods),0)
  
  #  ##save(feat_sel_matrix,file="feat_sel_matrix.Rda")
  
  #  ##save(feat.sel.methods,file="feat.sel.methods.Rda")
  if(is.na(good_feats_index)==TRUE){
    
    if(is.na(feat.sel.methods)==FALSE){
      
      X1orig=X
      Y1orig=Y
      
      for(m in 1:length(feat.sel.methods)){
        
        X=X1orig
        Y=Y1orig
        method=feat.sel.methods[m]
        
        v1<-matrix(ncol=dim(fiveCV10iter@learnmatrix)[1],nrow=dim(X)[2],0)
        
        rank_matrix<-matrix(nrow=dim(fiveCV10iter@learnmatrix)[1],ncol=dim(X)[2],0)
        
        #rank_matrix<-{}
        ranked_list_2<-{}
        
        if(method=="pls" || method=="spls" || method=="opls" || method=="ospls" || method=="rfboruta" || method=="rferadial" || method=="lmreg" || method=="lmregrobust" || method=="logitreg" || method=="logitregrobust" || method=="poissonregrobust" || method=="poissonreg"){
          
          #changing here
          for(i in 1:dim(fiveCV10iter@learnmatrix)[1]){
            
            
            X=X1orig
            Y=Y1orig
            temptrain<-t(X[fiveCV10iter@learnmatrix[i,],])
            tempclass<-Y[fiveCV10iter@learnmatrix[i,]]
            
            classlabels_sub<-cbind(paste("S",rep(1,length(tempclass)),sep=""),tempclass)
            
            classlabels_sub<-as.data.frame(classlabels_sub)
            
            if(method=="rfboruta"){
              
              varimp_res<-do_rf_boruta(X=temptrain,classlabels=tempclass,maxRuns=1000)
              
              varindex<-which(varimp_res>0)
              
              if(length(varindex)>0){
                v1[varindex,i]<-1
                
              }
              
              
            }else{
              
              if(method=="rferadial"){
                
                ranked_var_list<-diffexpsvmrfe(x=temptrain,y=tempclass,svmkernel="radial")
                
                
                varindex<-which(varimp_res<num.var.sel)
                if(length(varindex)>0){
                  v1[varindex,i]<-1
                  
                }
                
              }else{
                
                if(method=="lmreg"){
                  
                  lmregres<-runlmreg(X=temptrain,Y=tempclass,pvalue.thresh=pvalue.thresh,fdrmethod=fdrmethod,fdrthresh=fdrthresh)
                  
                  varindex<-lmregres$selected.index
                  if(length(varindex)>0){
                    v1[varindex,i]<-1
                    
                  }
                  
                  
                }else{
                  
                  if(method=="logitreg"){
                    
                    lmregres<-runlmreg(X=temptrain,Y=tempclass,pvalue.thresh=pvalue.thresh,fdrmethod=fdrmethod,fdrthresh=fdrthresh,logistic_reg=TRUE)
                    
                    varindex<-lmregres$selected.index
                    if(length(varindex)>0){
                      v1[varindex,i]<-1
                      
                    }
                    
                    
                  }else{
                    
                    if(method=="lmregrobust"){
                      
                      lmregres<-runlmreg(X=temptrain,Y=tempclass,pvalue.thresh=pvalue.thresh,fdrmethod=fdrmethod,fdrthresh=fdrthresh,robust.estimate=TRUE)
                      
                      varindex<-lmregres$selected.index
                      if(length(varindex)>0){
                        v1[varindex,i]<-1
                        
                      }
                      
                      
                    }else{
                      if(method=="logitregrobust"){
                        
                        lmregres<-runlmreg(X=temptrain,Y=tempclass,pvalue.thresh=pvalue.thresh,fdrmethod=fdrmethod,fdrthresh=fdrthresh,robust.estimate=TRUE,logistic_reg=TRUE)
                        
                        varindex<-lmregres$selected.index
                        if(length(varindex)>0){
                          v1[varindex,i]<-1
                          
                        }
                        
                        
                      }
                      
                    }
                    
                    
                  }
                  
                  
                }
              }
            }
            
            
            if(method=="spls"){
              sparseselect=TRUE
            }else{
              sparseselect=FALSE
            }
            
            var_rsd<-apply(temptrain,1,do_rsd)
            
            
            temptrain<-temptrain #[-which(var_rsd<1),]
            
            X=t(temptrain)
            Y=tempclass
            rownames(X)<-seq(1,nrow(X)) #colnames(temptrain)
            #Y=as.vector(Y)
            
            #Y<-t(Y)
            
            
            if(ncol(X)>0){
              numcomp<-5
              #return(list(X=X,Y=Y,numcomp=numcomp,kfold=kfold))
              
              set.seed(123)
              opt_comp<-pls.lda.cv(Xtrain=X, Ytrain=Y,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3, priors=NULL)
              
              if(method=="ospls"){
                
                Ytemp<-as.numeric(Y)
                leukemia.pls <- plsr(Ytemp ~ X, ncomp = opt_comp, validation = "LOO")
                ww <- leukemia.pls$loading.weights[,1]
                pp <- leukemia.pls$loadings[,1]
                w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
                t.ortho <- X %*% w.ortho
                
                p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
                Xcorr <- X - tcrossprod(t.ortho, p.ortho)
                
                
                
                X<-Xcorr
                method="spls"
              }
              
              if(method=="opls"){
                
                Ytemp<-as.numeric(Y)
                leukemia.pls <- plsr(Ytemp ~ X, ncomp = opt_comp, validation = "LOO")
                ww <- leukemia.pls$loading.weights[,1]
                pp <- leukemia.pls$loadings[,1]
                w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
                t.ortho <- X %*% w.ortho
                
                p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
                Xcorr <- X - tcrossprod(t.ortho, p.ortho)
                
                
                
                X<-Xcorr
                method="pls"
              }
              
              
              
              
              #opt_comp<-plsres1$opt_comp
              max_comp_sel<-opt_comp
              if(method=="spls"){
                
                keep_X_vec=rep(num.var.sel,opt_comp)
                
                linn.pls <- splsda(X, Y,ncomp=opt_comp,keepX=keep_X_vec)
                
                linn.vip<-linn.pls$loadings$X
                
                
                if(opt_comp>1){
                  
                  #abs
                  vip_res1<-abs(linn.vip)
                  
                  if(max_comp_sel>1){
                    vip_res1<-apply(vip_res1,1,mean)
                    
                  }else{
                    
                    vip_res1<-vip_res1[,c(1)]
                  }
                }else{
                  
                  vip_res1<-abs(linn.vip)
                }
                
                pls_vip<-vip_res1 #(plsres1$vip_res)
                
                
                #based on loadings for sPLS
                #feat_sel_matrix[which(pls_vip!=0),i]<-1 #pls_vip!=0 & rand_pls_sel_fdr<fdrthresh
                varindex<-which(pls_vip!=0)
                if(length(varindex)>0){
                  v1[varindex,i]<-1
                }
                
                
              }else{
                
                
                
                linn.pls <- plsda(X, Y,ncomp=opt_comp)
                
                linn.vip<-vip(linn.pls)
                #write.table(linn.vip,file="linn.vip.txt",sep="\t",row.names=TRUE)
                
                
                if(opt_comp>1){
                  vip_res1<-(linn.vip)
                  if(max_comp_sel>1){
                    vip_res1<-apply(vip_res1,1,mean)
                  }else{
                    
                    vip_res1<-vip_res1[,c(1)]
                  }
                }else{
                  
                  vip_res1<-linn.vip
                }
                
                
                
                #vip_res1<-plsres1$vip_res
                pls_vip<-vip_res1
                
                pls_vip_order<-pls_vip[order(pls_vip,decreasing=TRUE)]
                
                pls_vip_thresh<-min(pls_vip_order[1:num.var.sel])[1]
                
                rank_matrix[i,]<-t(order(pls_vip,decreasing=TRUE))
                
                
                ###save(rank_matrix,file="rank_matrix.Rda")
                # ##save(pls_vip,file="pls_vip.Rda")
                ###save(ranked_list,file="ranked_list.Rda")
                
                ###save(ranked_list_2,file="ranked_list_2.Rda")
                
                
                
                #print(pls_vip_thresh)
                #pls
                varindex<-which(pls_vip>=pls_vip_thresh)
                if(length(varindex)>0){
                  v1[varindex,i]<-1
                }
              }
            }
          }
          
          ranked_list<-rank_matrix
          for(rnum in 1:nrow(ranked_list)){
            
            ranked_list_2<-rbind(ranked_list_2,t(mzrt_id[ranked_list[rnum,]]))
            
          }
          
          
        }else{
          
          
          #set.seed(27611)
          set.seed(seedvalue)
          if(method=="rf"){
            g1<-GeneSelection(X=X,y=Y,learningsets=fiveCV10iter,method=method,trace=FALSE,seed = 100)
          }else{
            
            g1<-GeneSelection(X=X,y=Y,learningsets=fiveCV10iter,method=method,trace=FALSE)
          }
          gmatrix<-{}
          
          # ##save(g1,file="g1.Rda")
          ranked_list<-{}
          rank_matrix<-g1@rankings[[1]]
          
          v1<-matrix(nrow=dim(X)[2],ncol=dim(rank_matrix)[1],0)
          
          #    ##save(rank_matrix,file="rank_matrix.Rda")
          ###save(mzrt,file="mzrt.Rda")
          
          X<-Xtrain
          
          ranked_list<-g1@rankings[[1]]
          
          for(rnum in 1:nrow(ranked_list)){
            
            ranked_list_2<-rbind(ranked_list_2,t(mzrt_id[ranked_list[rnum,]]))
            
          }
          for(i in 1:dim(rank_matrix)[1]){
            
            varindex<-{}
            varindex1<-toplist(g1,iter=i,k=num.var.sel,show=FALSE)
            
            if(length(g1@rankings)>1){
              
              varindex<-c(varindex,varindex1[g1@rankings][,1])
              
            }else{
              varindex<-varindex1[,1]
            }
            
            varindex<-unique(varindex)
            
            
            v1[varindex,i]<-1
            
          }
        } #end else
        
        
        
        
        ###save(v1,file="v1.Rda")
        
        #hist(svm_acc,main="Inner test set accuracy distribution",col="brown")
        
        #iter.quantile.thresh: means that value is 1 in (1-iter.quantile.thresh)% or more sets;
        stability_measure<-apply(v1,1,function(x){length(which(x==1))/length(x)})  #quantile(x,iter.quantile.thresh)})
        
        stability_matrix<-stability_measure
        if(m==1){
          stability_matrix_1<-cbind(mzrt,stability_measure)
        }else{
          
          stability_matrix_1<-cbind(stability_matrix_1,stability_measure)
        }
        
        max_varsel=num.var.sel
        
        
        
        if(aggregation.method=="consensus"){
          
          feat_sel_matrix[which(stability_matrix>=iter.quantile.thresh),m]<-1
        }else{
          
          if(aggregation.method=="RankAggreg"){
            r1<-RankAggreg(x=ranked_list_2,k=max_varsel,verbose=TRUE,distance="Spearman",method="CE",maxIter=aggregation.max.iter)
          }else{
            
            if(aggregation.method=="RankAggregGA"){
              r1<-RankAggreg(x=ranked_list_2,k=max_varsel,verbose=TRUE,distance="Spearman",method="GA",maxIter=aggregation.max.iter)
            }else{
              
              
              feat_sel_matrix[which(stability_matrix==1),m]<-1
              
            }
          }
          # ##save(stability_matrix_1,file="stability_matrix1.Rda")
          
          stability_matrix_1<-as.data.frame(stability_matrix_1)
          #  ##save(r1,file="r1.Rda")
          
          if(length(check_name1)>0){
            
            colnames(stability_matrix_1)<-c("Name","stability_measure")
            mz_rt_all<-stability_matrix_1$Name
          }else{
            mz_rt_all<-paste(stability_matrix_1$mz,"_",stability_matrix_1$time,sep="")
          }
          
          common_row_index<-which(mz_rt_all%in%r1$top.list)
          
          
          feat_sel_matrix[common_row_index,m]<-1
          
          
          
        }
        
      }
      
      X=X1orig
      Y=Y1orig
      
      
      
    }else{
      
      feat_sel_matrix<-matrix(nrow=dim(Xtrain)[1],ncol=length(feat.sel.methods),1)
      
    }
    
    
    if(length(feat.sel.methods)>1){
      feat_sel_matrix<-apply(feat_sel_matrix,1,sum)
      
      
    }
    
    if(num.methods.sel>length(feat.sel.methods)){
      
      num.methods.sel=length(feat.sel.methods)
    }
    
    if(length(check_name1)>0){
      
      colnames(stability_matrix_1)<-c("mzrt",feat.sel.methods)
      
    }else{
      
      colnames(stability_matrix_1)<-c("mz","time",feat.sel.methods)
    }
    
    #pdf("Results.pdf")
    
    #hist(stability_measure,main="Stability measure distribution",col="brown")
    write.table(stability_matrix_1,file="stability_matrix.txt",sep="\t",row.names=FALSE)
    
    good_feats_index<-which(feat_sel_matrix>=num.methods.sel)
    
    if(FALSE){
      if(is.na(pvalue.thresh)==FALSE){
        
        
        
        
        numcores<-num_nodes #round(detectCores()*0.5)
        
        cl <- parallel::makeCluster(getOption("cl.cores", numcores))
        
        clusterExport(cl,"diffexponewayanova",envir = .GlobalEnv)
        
        clusterExport(cl,"anova",envir = .GlobalEnv)
        
        
        clusterExport(cl,"TukeyHSD",envir = .GlobalEnv)
        
        clusterExport(cl,"aov",envir = .GlobalEnv)
        
        
        #res1<-apply(data_m_fc,1,function(x){
        res1<-parApply(cl,X,2,function(x,classlabels_response_mat){
          xvec<-x
          
          
          data_mat_anova<-cbind(xvec,classlabels_response_mat)
          
          data_mat_anova<-as.data.frame(data_mat_anova)
          cnames<-colnames(data_mat_anova)
          
          cnames[1]<-"Response"
          
          colnames(data_mat_anova)<-c("Response","Factor1")
          
          #print(data_mat_anova)
          
          data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
          
          anova_res<-diffexponewayanova(dataA=data_mat_anova)
          
          
          
          return(anova_res)
        },Y)
        
        stopCluster(cl)
        main_pval_mat<-{}
        
        posthoc_pval_mat<-{}
        pvalues<-{}
        
        
        
        for(i in 1:length(res1)){
          
          
          main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
          pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
          posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthocfactor1)
          
        }
        
        pvalues<-unlist(pvalues)
        
        good_feats_index<-which(feat_sel_matrix>=num.methods.sel & pvalues<pvalue.thresh)
        if(length(good_feats_index)<1){
          stop("No features selected.")
        }
        
      }
      
    }
    good_feats<-Xtrain[,good_feats_index]
    
    
    if(length(check_name1)>0){
      mzrt_sub<-mzrt[good_feats_index]
    }else{
      mzrt_sub<-mzrt[good_feats_index,]
    }
    good_feats<-t(good_feats)
    
    
    
    cnames_1<-colnames(good_feats)
    
    
    #colnames(good_feats)<-sampnames
    
    
    good_feats<-cbind(mzrt_sub,good_feats)
    
    
    
    ###save(sampnames_test,mzrt_sub,good_feats,Xtest,X,Y,Xtrain,mzrt,good_feats_index,file="good_feats_index.Rda")
    
    if(is.na(Xtest)==FALSE){
      good_feats_test<-Xtest[,good_feats_index]
      good_feats_test<-t(good_feats_test)
      colnames(good_feats_test)<-sampnames_test
      good_feats_test<-cbind(mzrt_sub,good_feats_test)
    }else{
      good_feats_test<-NA
      Ytest_mat<-NA
    }
    
    
  }else{
    pdf("Results.pdf")
    
    good_feats<-Xtrain[,good_feats_index]
    
    
    
    mzrt_sub<-mzrt[good_feats_index,]
    
    good_feats<-t(good_feats)
    
    
    
    cnames_1<-colnames(good_feats)
    
    
    colnames(good_feats)<-sampnames
    
    
    good_feats<-cbind(mzrt_sub,good_feats)
    
    if(is.na(Xtest)==FALSE){
      good_feats_test<-Xtest[,good_feats_index]
      good_feats_test<-t(good_feats_test)
      colnames(good_feats_test)<-sampnames_test
      good_feats_test<-cbind(mzrt_sub,good_feats_test)
    }else{
      good_feats_test<-NA
      Ytest_mat<-NA
    }
    
  }
  
  print("Number of features selected")
  print(length(good_feats_index))
  
  #    #save(good_feats,Xtrain,good_feats_test,Xtest,Ytest,Y,good_feats_index,stability_matrix_1,file="feat.sel.res.Rda")
  
  
  if(length(good_feats_index)>1){
    
    X<-X[,good_feats_index]
    
    if(is.na(Xtest)==FALSE){
      
      Xtest<-Xtest[,good_feats_index]
    }
    
    #  #save(X,Y,fiveCV10iter,Xtest,Ytest,good_feats_test,Ytest,file="acc.Rda")
    
    return(list("selected.features.index"=good_feats_index,"train.select"=X,"test.select"=Xtest,"stability_matrix"=stability_matrix_1,"feat.sel.matrix"=feat_sel_matrix,train.class=Ytrain_mat,test.class=Ytest_mat))
  }
}

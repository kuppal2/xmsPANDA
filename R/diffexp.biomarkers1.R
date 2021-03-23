diffexp.biomarkers1 <-
function(X=NA,Y=NA,feature_table_file=NA,class_labels_file=NA,
                              feat.sel.methods=c("rf","rfe","pls","limma","lasso"),num.var.sel=10,prop.select.thresh=0.7,
                              split.train.test=FALSE,train.pct=0.7,outloc,kfold=10,pca.ellipse=TRUE,Xtest=NA,Ytest=NA,
                              rsdthresh=1,pls_vip_thresh=2,seedvalue=27611,learningsetmethod="CV",confounder.matrix=NA,
                              fdrmethod="BH",fdrthresh=0.05,num.methods.sel=1,globalcor=FALSE,
                              cor.method="spearman",networktype="complete",abs.cor.thresh=0.4,
                              cor.fdrthresh=0.05,max.cor.num=100,net_node_colors=c("green","red"),
                              net_legend=TRUE,niter=10,output.device.type="pdf",
                              heatmap.col.opt="RdBu",boxplot.col.opt=c("white"),
                              barplot.col.opt=c("grey57","grey90"),sample.col.opt="rainbow",
                              mz.thresh=1,time.thresh=10,svm_kernel="radial",
                              good_feats_index=NA,pvalue.thresh=0.05,plots.width=8,
                              plots.height=8,plots.res=600, plots.type="cairo",
                              num_nodes=2,ylabel="Intensity",cex.plots=0.7,tune_classifiers=FALSE,
                              find.common.features=FALSE,aggregation.method="consensus",
                              aggregation.max.iter=1000,add.pvalues=TRUE,add.jitter=TRUE,
                              alphabetical.order=FALSE,boxplot.type="ggplot",balance.classes=FALSE,deeplearning=FALSE,alpha.col=1,ggplot.type1=NA)
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
      
      suppressMessages(library(CMA))
      test_data_check=1 
      numtrain<-round(train.pct*nsamp)
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
    #   #save(Xtrain,file="Xtrain.Rda")
    ##save(Xtest,file="Xtest.Rda")
    ##save(Y,file="Y.Rda")
    ##save(Ytest,file="Ytest.Rda")
    
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
    
    pdf("Results.pdf")
    
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
  
  #  #save(good_feats,Xtrain,good_feats_test,Xtest,Ytest,Y,good_feats_index,stability_matrix_1,file="feat.sel.res.Rda")
  
  
  if(length(good_feats_index)>1){
    
    X<-X[,good_feats_index]
    
    if(is.na(Xtest)==FALSE){
      
      Xtest<-Xtest[,good_feats_index]
    }
    
    #    #save(X,Y,fiveCV10iter,Xtest,Ytest,good_feats_test,Ytest,file="acc.Rda")
    
    return(good_feats_index)
    
    if(tune_classifiers==TRUE)
    {
      
      print("Building classification model using training set without tuning")
      
      if(length(good_feats_index)>=3)
      {
        
        
        set.seed(seedvalue)
        tune_plslda <- suppressWarnings(CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, grids = list(comp = 1:10),trace=FALSE))
        
        #     ##save(tune_plslda,file="tune_plslda.Rda")
        
        set.seed(seedvalue)
        tune_plsrf <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, grids = list(comp = 1:10),trace=FALSE)
        
        ##save(tune_plsrf,file="tune_plsrf.Rda")
      }
      
      ##saveX,Y,fiveCV10iter,file="scda.Rda")
      set.seed(seedvalue)
      tune_scda <-try(CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, grids = list( ),trace=FALSE),silent=TRUE)
      
      ##save(tune_scda,file="tune_scda.Rda")
      
      if(is(tune_scda,"try-error")){
        
        #t1<-new("list")
        #t1<-as.list(rep(0.5,nrow(fiveCV10iter@learnmatrix)))
        #  tune_scda<-new("tuningresult",tuneres=t1,method="scDA")
        
        
      }
      
      set.seed(seedvalue)
      tune_svm <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, grids = list( ), kernel = "radial",trace=FALSE)
      
      ##save(tune_svm,file="tune_svm.Rda")
      set.seed(seedvalue)
      #  tune_plr <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA, grids = list( ), trace=FALSE)
      
      
      set.seed(seedvalue)
      
      #if(FALSE)
      {
        if(dim(X)[2]>100){
          tune_nnet <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA, trace=FALSE) #grids = list(size=1:5,decay=c(0,0.0001,0.001,0.005, 0.01,0.05, 0.1)), trace=FALSE)
        }else{
          tune_nnet <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA, trace=FALSE) # grids = list(size=1:5,decay=c(0,0.0001,0.001,0.005, 0.01,0.05, 0.1)), trace=FALSE)
          
        }
        
        ##save(tune_nnet,file="tune_nnet.Rda")
      }
      
      # ##save(X,Y,fiveCV10iter,nnetCMA,tune_nnet,learnmatrix,file="Debug.rda")
      
      set.seed(seedvalue)
      class_nnet<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE,tuneres=tune_nnet)
      
      
      
      #  set.seed(seedvalue)
      #tune_rf <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA, grids = list( ), trace=FALSE)
      
      if(length(good_feats_index)>2){
        set.seed(seedvalue)
        class_plslda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, tuneres = tune_plslda,trace=FALSE)
        set.seed(seedvalue)
        class_plsrf <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, tuneres = tune_plsrf,trace=FALSE)
      }
      
      set.seed(seedvalue)
      
      if(is(tune_scda,"try-error")){
        
        class_scda <- NA
      }else{
        class_scda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, tuneres = tune_scda,trace=FALSE)
      }
      
      
      set.seed(seedvalue)
      class_svm <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres = tune_svm,kernel = "radial",trace=FALSE,probability=TRUE)
    }else{
      
      print("Building classification model using training set without tuning")
      if(length(good_feats_index)>2){
        set.seed(seedvalue)
        class_plslda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, trace=FALSE)
        set.seed(seedvalue)
        class_plsrf <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, trace=FALSE)
      }
      
      set.seed(seedvalue)
      class_scda <- try(classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, trace=FALSE),silent=TRUE)
      
      if(is(class_scda,"try-error")){
        
        class_scda <- NA
        
      }
      
      
      set.seed(seedvalue)
      class_svm <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, kernel = "radial",trace=FALSE,probability=TRUE)
      
      set.seed(seedvalue)
      class_nnet<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE)
    }
    
    
    
    set.seed(seedvalue)
    class_rf <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE) #tuneres = tune_rf,
    #size=3,decay=0.1) #
    set.seed(seedvalue)
    class_plr<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE) #tuneres=tune_plr,
    if(length(good_feats_index)>2){ 
      ##save(class_plslda,file="class_plslda.Rda")
      
      ##save(class_plsrf,file="class_plsrf.Rda")
    }
    ##save(class_scda,file="class_scda.Rda")
    ##save(class_svm,file="class_svm.Rda")
    ##save(class_rf,file="class_rf.Rda")
    ##save(class_nnet,file="class_nnet.Rda")
    ##save(class_plr,file="class_plr.Rda")
    
    
    
    if(length(class_levels)==2){
      
      
      eval_measure="auc"
      
      
    }else{
      
      eval_measure="misclassification"
    }
    
    
    
    
    if(length(good_feats_index)>=2){
      
      
      
      
      
      if(is.na(class_scda)==FALSE){
        
        
        eval_scda<-try(evaluation(class_scda,measure=eval_measure),silent=TRUE)
        if(is(eval_scda,"try-error")){
          
          eval_scda<-NA
          class_scda<-NA
          
        }
        
      }else{
        eval_scda<-NA
        class_scda<-NA
      }
      
      
      ###save(class_scda,eval_scda,class_svm,class_rf,class_nnet,class_plr,eval_measure,X,Y,fiveCV10iter,file="scda_debug1.Rda")
      
      eval_svm<-try(evaluation(class_svm,measure=eval_measure),silent=TRUE)
      eval_rf<-try(evaluation(class_rf,measure=eval_measure),silent=TRUE)
      eval_nnet<-try(evaluation(class_nnet,measure=eval_measure),silent=TRUE)
      eval_plr<-evaluation(class_plr,measure=eval_measure)
      
      if(is(eval_nnet,"try-error")){
        
        eval_nnet<-new("evaloutput",score=1)
        
        
      }
      
      if(is(eval_rf,"try-error")){
        
        eval_rf<-new("evaloutput",score=1)
        
        
      }
      if(is(eval_svm,"try-error")){
        
        eval_svm<-new("evaloutput",score=1)
        
        
      }
      
      if(is(eval_plr,"try-error")){
        
        eval_plr<-new("evaloutput",score=1)
        
        
      }
      
      
      
      
      if(length(good_feats_index)>2){
        eval_plslda<-evaluation(class_plslda,measure=eval_measure)
        eval_plsrf<-evaluation(class_plsrf,measure=eval_measure)
        
        
      }else{
        eval_plslda<-eval_svm
        eval_plsrf<-eval_svm
        eval_plslda@score<-(1)
        eval_plsrf@score<-(1)
      }
      
      if(eval_measure=="auc")
      {
        eval_plslda@score=1-mean(eval_plslda@score)
        eval_plsrf@score=1-mean(eval_plsrf@score)
        
        if(is.na(class_scda)==FALSE){
          eval_scda@score=1-mean(eval_scda@score)
        }else{
          
          eval_scda=eval_svm
          eval_scda@score=(1)
        }
        eval_svm@score=1-mean(eval_svm@score)
        eval_rf@score=1-mean(eval_rf@score)
        eval_nnet@score=1-mean(eval_nnet@score)
        eval_plr@score=1-mean(eval_plr@score)
        
      }
    }
    
    ####saveeval_svm,file="eval_svm.Rda")
    ####saveclass_svm,file="class_svm.Rda")
    
    
    text2<-paste(dim(v1)[2], " learning sets using training data",sep="")
    
    
    if(length(class_levels)==0){
      
      eval_mat1<-cbind(text2,100*(1-mean(eval_plslda@score)),100*(1-mean(eval_plslr@score)),100*(1-mean(eval_plsrf@score)),100*(1-mean(eval_scda@score)),100*(1-mean(eval_svm@score)),100*(1-mean(eval_rf@score)),100*(1-mean(eval_nnet@score)),100*(1-mean(eval_plr@score)),100*(1-mean(eval_lassoplr@score)),100*(1-mean(eval_elasticnetplr@score)))
      
      best_classifier<-which(eval_mat1==max(eval_mat1))
      classifier_names<-c("PLSLDA","PLSLR","PLSRF","SCDA","SVM","RF","NNet","pLR","pLRlasso","pLRelasticnet")
      
      best_classifier_name<-classifier_names[best_classifier]
      
      eval_mat1<-cbind(text2,100*(1-mean(eval_plslda@score)),100*(1-mean(eval_plslr@score)),100*(1-mean(eval_plsrf@score)),100*(1-mean(eval_scda@score)),100*(1-mean(eval_svm@score)),100*(1-mean(eval_rf@score)),100*(1-mean(eval_nnet@score)),100*(1-mean(eval_plr@score)),100*(1-mean(eval_lassoplr@score)),100*(1-mean(eval_elasticnetplr@score)))
      
      colnames(eval_mat1)<-c("Dataset","PLSLDA","PLSLR","PLSRF","SCDA","SVM","RF","NNet","pLR","pLRlasso","pLRelasticnet")
    }else{
      
      if(length(good_feats_index)>2){
        eval_mat1<-cbind(100*(1-mean(eval_plslda@score)),100*(1-mean(eval_plsrf@score)),100*(1-mean(eval_scda@score)),100*(1-mean(eval_svm@score)),100*(1-mean(eval_rf@score)),100*(1-mean(eval_nnet@score)),100*(1-mean(eval_plr@score)))
      }else{
        if(is.na(class_scda)==FALSE){
          eval_scda@score=1-mean(eval_scda@score)
        }else{
          
          eval_scda=eval_svm
          eval_scda@score=(1)
        }
        eval_mat1<-cbind(100*(1-1),100*(1-1),100*(1-mean(eval_scda@score)),100*(1-mean(eval_svm@score)),100*(1-mean(eval_rf@score)),100*(1-mean(eval_nnet@score)),100*(1-mean(eval_plr@score)))
        
      }
      
      best_classifier<-which(eval_mat1==max(eval_mat1)[1])
      
      classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
      
      best_classifier_name<-classifier_names[best_classifier]
      
      best_innerkfold_acc<-max(eval_mat1)[1]
      
      eval_mat_1<-round(eval_mat1,2)
      eval_mat1<-cbind(text2,eval_mat_1)
      
      colnames(eval_mat1)<-c("Dataset","PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
      
    }
    
    if(output.device.type!="pdf"){
      
      temp_filename_1<-"Figures/Barplot_classifier_comparison_CVaccuracy.png"
      
      png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    }
    
    
    w <- 0.1 #grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
    par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
    barplot(eval_mat_1,beside=TRUE,main="Comparison of CV accuracy using different classifiers\n based on learning sets",
            xlab="Classifier",ylab="kfold classification accuracy(%)",col=barplot.col.opt[1],type="p",ylim=c(0,100),xpd=FALSE,cex.axis=0.7,cex.names=0.7)
    
    
    
    if(output.device.type!="pdf"){
      
      try(dev.off(),silent=TRUE)
    }
    
    if(length(check_name1)>0){
      mzrt_1<-good_feats[,1]
    }else{
      mzrt_1<-paste(round(good_feats[,1],5),round(good_feats[,2],1),sep="_")
      
    }
    rownames(good_feats)<-mzrt_1
    
    
    Y1<-cbind(sampnames,as.character(Y))
    Y1<-as.data.frame(Y1)
    
    if(length(good_feats_index)>=3){
      
      if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/HCA_selectedfeats.png"
        
        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      }
      
      good_feats<-as.data.frame(good_feats)
      
      if(FALSE){
        if(length(check_name1)>0){
          
          cnames1<-colnames(good_feats)
          cnames1[1]<-c("Name")
          colnames(good_feats)<-cnames1
          
          good_feats[,-c(1)]<-apply(good_feats[,-c(1)],2,as.numeric)
          
        }else{
          
          cnames1<-colnames(good_feats)
          cnames1[1]<-c("Name")
          colnames(good_feats)<-cnames1
          good_feats[,-c(1:2)]<-apply(good_feats[,-c(1:2)],2,as.numeric)
        }
        
      }
      
      
      
      
      try(get_hca(parentoutput_dir=outloc,X=good_feats,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
                  sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,
                  newdevice=FALSE,labRow.value = labRow.value, labCol.value = labCol.value,similarity.matrix=similarity.matrix,cexLegend=hca.cex.legend,cexRow=cex.plots,cexCol=cex.plots),silent=TRUE)
      
      
      
      if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
      }
      
      
      if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/PCAplots_selectedfeats.pdf"
        
        #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        #pdf(temp_filename_1)
        
        pdf(temp_filename_1,width=plots.width,height=plots.height)
      }
      
      
    }
    
    
    if(output.device.type!="pdf"){
      
      try(dev.off(),silent=TRUE)
    }
    
    
    if(output.device.type!="pdf"){
      
      temp_filename_1<-"Figures/Boxplots_selectedfeats.pdf"
      
      #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      #pdf(temp_filename_1)
      pdf(temp_filename_1,width=plots.width,height=plots.height)
    }
    
    
    
    
    
    par(mfrow=c(1,1),family="sans",cex=cex.plots)
    get_boxplots(X=good_feats,Y=Y1,parentoutput_dir=outloc,boxplot.col.opt=boxplot.col.opt,
                 sample.col.opt=sample.col.opt,
                 newdevice=FALSE,cex.plots=cex.plots,
                 ylabel=ylabel,add.pvalues=add.pvalues,add.jitter=add.jitter,
                 boxplot.type=boxplot.type,study.design=analysistype,multiple.figures.perpanel=multiple.figures.perpanel,alphacol = alpha.col,ggplot.type1=ggplot.type1,facet.nrow=facet.nrow)
    
    
    
    
    
    if(output.device.type!="pdf"){
      
      try(dev.off(),silent=TRUE)
    }
    
    num_levels<-levels(as.factor(Y))
    
    
    
    
    #if(is.na(Xtest)==FALSE)
    {
      
      
      
      best_kfold_acc<-best_innerkfold_acc
      #outerkfold_acc<-100*(1-mean(testeval_res@score))
      
      permkfold_acc<-{}
      permkfold_acc1<-{}
      
      permkfold_acc2<-{}
      
      permkfold_acc3<-{}
      
      permkfold_acc4<-{}
      
      permkfold_acc5<-{}
      
      permkfold_acc6<-{}
      
      permkfold_acc7<-{}
      Yorig<-Y
      
      text2B<-paste("Permuted ",dim(v1)[2], " learning sets using training data",sep="")
      #if(is.na(Xtest)==FALSE)
      
      nperm<-3
      seedvalue_rand_list<-runif(nperm,1,10000)
      #for(p1 in 1:nperm)
      eval_mat_perm<-lapply(1:nperm,function(p1)
      {
        seedvalue_cur=round(seedvalue_rand_list[p1],0)
        #set.seed(27611)
        set.seed(seedvalue_cur)
        Y<-Yorig[sample(1:length(Yorig),size=length(Yorig))]
        #set.seed(27611)
        set.seed(seedvalue_cur)
        suppressMessages(library(CMA))
        fiveCV10iter<-GenerateLearningsets(y=Y,method=learningsetmethod,fold=kfold,niter=1,strat=TRUE)
        
        classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR") #,"PLSLR","pLRlasso","pLRelasticnet")
        
        
        if(length(good_feats_index)>=3){
          
          set.seed(seedvalue)
          class_plslda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA,trace=FALSE) #,tuneres=tune_plslda
          testeval_res1<-evaluation(class_plslda,measure=eval_measure)
          
          set.seed(seedvalue)
          class_plsrf <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA,trace=FALSE) #,tuneres=tune_plsrf
          testeval_res2<-evaluation(class_plsrf,measure=eval_measure)
          
        }
        set.seed(seedvalue)
        
        
        if(is.na(class_scda)==FALSE){
          
          class_scda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA,trace=FALSE) #,tuneres=tune_scda
          testeval_res3<-try(evaluation(class_scda,measure=eval_measure),silent=TRUE)
          
          if(is(testeval_res3,"try-error")){
            
            testeval_res3<-NA
          }
          
        }else{
          testeval_res3<-NA
        }
        
        
        set.seed(seedvalue)
        class_svm <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA,kernel = "radial",trace=FALSE,probability=TRUE) #,tuneres=tune_svm
        testeval_res4<-evaluation(class_svm,measure=eval_measure)
        
        set.seed(seedvalue)
        class_rf <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE) #,tuneres=tune_rf
        testeval_res5<-evaluation(class_rf,measure=eval_measure)
        
        set.seed(seedvalue)
        class_nnet<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE)
        #size=3,decay=0.1) #,tuneres=tune_nnet
        
        testeval_res6<-evaluation(class_nnet,measure=eval_measure)
        
        set.seed(seedvalue)
        class_plr<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE) #,tuneres=tune_plr
        testeval_res7<-evaluation(class_plr,measure=eval_measure)
        
        
        if(eval_measure=="auc")
        {
          eval_plslda@score=1-mean(eval_plslda@score)
          eval_plsrf@score=1-mean(eval_plsrf@score)
          
          if(is.na(class_scda)==FALSE){
            eval_scda@score=1-mean(eval_scda@score)
          }
          eval_svm@score=1-mean(eval_svm@score)
          eval_rf@score=1-mean(eval_rf@score)
          eval_nnet@score=1-mean(eval_nnet@score)
          eval_plr@score=1-mean(eval_plr@score)
          
        }
        
        
        
        if(length(good_feats_index)>=3){
          permkfold_acc1<-c(permkfold_acc1,100*(1-mean(testeval_res1@score)))
          
          permkfold_acc2<-c(permkfold_acc2,100*(1-mean(testeval_res2@score)))
          
        }else{
          permkfold_acc1<-c(permkfold_acc1,0)
          permkfold_acc2<-c(permkfold_acc2,0)
          testeval_res1<-testeval_res4
          testeval_res2<-testeval_res4
          testeval_res1@score<-1
          testeval_res2@score<-1
        }
        
        if(is.na(class_scda)==FALSE){
          
          if(is.na(testeval_res3)==FALSE){
            
            permkfold_acc3<-100*(1-mean(testeval_res3@score))
          }else{
            
            permkfold_acc3<-(0)
          }
          
          
        }else{
          
          permkfold_acc3<-(0)
        }
        
        permkfold_acc4<-c(permkfold_acc4,100*(1-mean(testeval_res4@score)))
        
        permkfold_acc5<-c(permkfold_acc5,100*(1-mean(testeval_res5@score)))
        
        permkfold_acc6<-c(permkfold_acc6,100*(1-mean(testeval_res6@score)))
        
        permkfold_acc7<-c(permkfold_acc7,100*(1-mean(testeval_res7@score)))
        
        temp_res<-c(100*(1-mean(testeval_res1@score)),100*(1-mean(testeval_res2@score)),permkfold_acc3,
                    100*(1-mean(testeval_res4@score)),100*(1-mean(testeval_res5@score)),100*(1-mean(testeval_res6@score)),
                    100*(1-mean(testeval_res7@score)))
        
        return(temp_res)
      })
      
      
      
      eval_mat_perm<-do.call(rbind,eval_mat_perm)
      eval_mat_perm<-apply(eval_mat_perm,2,mean)
      
      
      eval_mat_perm<-t(eval_mat_perm)
      eval_mat_perm<-round(eval_mat_perm,2)
      
      eval_mat_perm_final<-cbind(text2B,eval_mat_perm)
      
      
      colnames(eval_mat_perm_final)<-c("Dataset","PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
      
      eval_mat_actual<-as.data.frame(eval_mat_1)
      eval_mat_perm1<-as.data.frame(eval_mat_perm)
      eval_mat_perm1<-round(eval_mat_perm1,2)
      
      
      colnames(eval_mat_perm1)<-colnames(eval_mat_actual)
      emat1<-rbind(eval_mat_actual,eval_mat_perm1)
      
      emat1<-t(emat1)
      classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
      
      rownames(emat1)<-classifier_names
      
      
      eval_mat_actual<-((emat1[,1]-emat1[,2])*0.5)+(0.5*(emat1[,1]))
      
      emat1<-cbind(emat1,eval_mat_actual)
      #[((Actual-Permuted)*0.5)+(0.5*Actual)]
      colnames(emat1)<-c("Actual.accuracy","Permuted.accuracy","Score")
      
      
      
      
      
      test_acc<-NA
      test_acc_mat<-{}
      Y<-Yorig
      if(is.na(Xtest)==FALSE){
        
        
        if(num_levels==2){
          if(split.train.test==TRUE){
            
            res<-get_classification.accuracy(kfold=kfold,featuretable=good_feats,classlabels=Y,classifier="logit",testfeaturetable=good_feats_test,testclasslabels=Ytest,errortype="BAR",kernelname=svm_kernel,svm.cost=svm.cost,svm.gamma=svm.gamma)
            
            
          }else{
            #get_roc(dataA=good_feats,classlabels=Y,classifier="svm",kname=svm_kernel,rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,testset=X,testclasslabels=Y,mainlabel="Using training set based on SVM")
            
            #get_roc(dataA=good_feats,classlabels=Y,classifier="svm",kname=svm_kernel,rocfeatlist=c(dim(good_feats)[2]),rocfeatincrement=FALSE,testset=good_feats_test,testclasslabels=Ytest,mainlabel="Test set")
            
            res<-get_classification.accuracy(kfold=kfold,featuretable=good_feats,classlabels=Y,classifier="logit",testfeaturetable=good_feats_test,testclasslabels=Ytest,errortype="BAR",kernelname=svm_kernel,svm.cost=svm.cost,svm.gamma=svm.gamma)
            
            
          }
        }
        
        test_acc_mat<-{}
        class_levels_vec<-levels(as.factor(Y))
        
        
        learnmatrix <- matrix(seq(1,nrow(X)), nrow = 1)
        fiveCV10iter<-new("learningsets", learnmatrix = learnmatrix, method = "none",ntrain = ncol(learnmatrix), iter = nrow(learnmatrix))
        X<-rbind(X,Xtest)
        Y<-c(Y,Ytest)
        classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR","H2O.deeplearning","PLSLR","pLRlasso","pLRelasticnet")
        
        result_cma_list<-new("list")
        
        ##save(fiveCV10iter,X,Y,learnmatrix,file="debug1.Rda")
        ###savelist=ls(),file="debug3.Rda")
        
        confusion_matrix_list<-new("list")
        
        if(length(good_feats_index)>2)
        {
          
          if(tune_classifiers==TRUE){
            s1<-ldply(tune_plslda@tuneres,rbind)
            
            s2<-apply(s1,2,median)
            
            
            t1<-new("list")
            confusion_matrix_list<-new("list")
            t1[[1]]<-s2
            tune_plslda1<-tune_plslda
            tune_plslda1@tuneres<-t1
            
            set.seed(seedvalue)
            class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, tuneres = tune_plslda1,trace=FALSE)
            
            b1<-best(tune_plslda1)
            
            
            
            learnmatrix<-as.numeric(learnmatrix)
            
            set.seed(seedvalue)
            class_res2<-pls_ldaCMA(X = X, y = Y, learnind = learnmatrix, comp = median(unlist(b1)))
            
            
          }else{
            
            set.seed(seedvalue)
            class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, trace=FALSE)
            learnmatrix<-as.numeric(learnmatrix)
            
            set.seed(seedvalue)
            class_res2<-pls_ldaCMA(X = X, y = Y, learnind = learnmatrix)
            
            
          }
          ##saveclass_res2,file="pls_ldaCMA.Rda")
          result_cma_list[[1]]<-class_res2
          
          ###saveconfusion_matrix_1,file="confusion_matrix_plslda.Rda")
          confusion_matrix_list[[1]]<-table(class_res2@y,class_res2@yhat)
          
          #debughere
          if(length(class_levels_vec)==2){
            testeval_res_auc<-evaluation(class_res,measure = "auc")
            ###savetesteval_res_auc,file="testeval_res_auc.Rda")
            
            test_auc<-100*(mean(testeval_res_auc@score))
            test_acc_mat<-c(test_acc_mat,test_auc)
            
          }else{
            
            test_acc<-evaluation(class_res)
            
            test_acc<-100*(1-mean(testeval_res_auc@score))
            test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            
          }
          
          
          
          
          #if(best_classifier_name==classifier_names[2]){
          
          if(tune_classifiers==TRUE){
            s1<-ldply(tune_plsrf@tuneres,rbind)
            s2<-apply(s1,2,median)
            t1<-new("list")
            t1[[1]]<-s2
            
            tune_plsrf1<-tune_plsrf
            tune_plsrf1@tuneres<-t1
            set.seed(seedvalue)
            class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, tuneres = tune_plsrf1,trace=FALSE)
            
            b1<-best(tune_plsrf1)
            # ###saveb1,file="b1_pls_rfCMA.Rda")
            
            learnmatrix<-as.numeric(learnmatrix)
            
            set.seed(seedvalue)
            class_res2<-pls_rfCMA(X = X, y = Y, learnind = learnmatrix, comp = median(unlist(b1)))
          }else{
            set.seed(seedvalue)
            class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, trace=FALSE)
            learnmatrix<-as.numeric(learnmatrix)
            
            set.seed(seedvalue)
            class_res2<-pls_rfCMA(X = X, y = Y, learnind = learnmatrix)
            
            
          }
          
          ##saveclass_res2,file="pls_rfCMA.Rda")
          result_cma_list[[2]]<-class_res2
          #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
          ###saveconfusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
          
          confusion_matrix_list[[2]]<-table(class_res2@y,class_res2@yhat)
          
          if(length(class_levels_vec)==2){
            testeval_res_auc<-evaluation(class_res,measure = "auc")
            ###savetesteval_res_auc,file="testeval_res_auc.Rda")
            
            test_auc<-100*(mean(testeval_res_auc@score))
            test_acc_mat<-c(test_acc_mat,test_auc)
            
          }else{
            
            test_acc<-evaluation(class_res)
            
            test_acc<-100*(1-mean(testeval_res_auc@score))
            test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            
          }
          
          
          
          if(tune_classifiers==TRUE){
            
            s1<-ldply(tune_scda@tuneres,rbind)
            s2<-apply(s1,2,median)
            
            t1<-new("list")
            t1[[1]]<-s2
            
            tune_scda1<-tune_scda
            tune_scda1@tuneres<-t1
            
            #tune_scda<-new("tuningresult",tuneres=t1)
            
            set.seed(seedvalue)
            class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, tuneres = tune_scda1,trace=FALSE)
            #scdaCMA(X, y, f, learnind, delta = 0.5, models=FALSE,...)
            
            b1<-best(tune_scda1)
            
            
            learnmatrix<-as.numeric(learnmatrix)
            set.seed(seedvalue)
            class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix, delta = median(unlist(b1)))
            
          }else{
            
            if(is.na(class_scda)==FALSE){
              set.seed(seedvalue)
              class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA,trace=FALSE)
              set.seed(seedvalue)
              class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix)
              
            }
          }
          
          if(is.na(class_scda)==FALSE){
            
            ##saveclass_res2,file="scdaCMA.Rda")
            result_cma_list[[3]]<-class_res2
            #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_scdaCMA.Rda")
            confusion_matrix_list[[3]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ###savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              test_auc<-100*(mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_auc)
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
            
          }else{
            result_cma_list[[3]]<-NA
            confusion_matrix_list[[3]]<-NA
            test_acc_mat<-c(test_acc_mat,0)
            
          }
          
          
          
          #if(best_classifier_name==classifier_names[4])
          {
            
            if(tune_classifiers==TRUE){
              s1<-ldply(tune_svm@tuneres,rbind)
              
              s2<-apply(s1,2,median)
              
              t1<-new("list")
              t1[[1]]<-s2
              tune_svm1<-tune_svm
              tune_svm1@tuneres<-t1
              
              b1<-best(tune_svm1)
              
              set.seed(seedvalue)
              class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
              
              ##save(b1,file="b1.Rda")
              
              
              b2=unlist(unlist(b1))
              
              learnmatrix<-as.numeric(learnmatrix)
              
              set.seed(seedvalue)
              
              
              class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix,probability=TRUE,cost=median(b2[seq(1,length(b2),2)]),gamma=median(b2[seq(2,length(b2),2)])) #,gamma=gamma_1,cost=cost_1)
              
              
              
            }else{
              
              set.seed(seedvalue)
              class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, kernel ="radial",trace=FALSE,probability=TRUE)
              set.seed(seedvalue)
              class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix,probability=TRUE)
              
            }
            
            ##saveclass_res2,file="svmCMA.Rda")
            result_cma_list[[4]]<-class_res2
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_svmCMA.Rda")
            confusion_matrix_list[[4]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ###savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              test_auc<-100*(mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_auc)
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
            
            
          }
          
          #if(best_classifier_name==classifier_names[5])
          {
            
            if(FALSE){
              s1<-ldply(tune_rf@tuneres,rbind)
              s2<-apply(s1,2,median)
              
              t1<-new("list")
              t1[[1]]<-s2
              
              tune_rf1<-tune_rf
              tune_rf1@tuneres<-t1
            }
            
            set.seed(seedvalue)
            class_res <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE) #,tuneres=tune_rf1
            
            learnmatrix<-as.numeric(learnmatrix)
            
            
            set.seed(seedvalue)
            class_res2 <- rfCMA(X =X, y = Y, learnind=learnmatrix, varimp = FALSE) #mtry=tune_rf1$mtry,nodesize=tune_rf1$nodesize
            
            ##saveclass_res2,file="rfCMA.Rda")
            result_cma_list[[5]]<-class_res2
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
            confusion_matrix_list[[5]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ###savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              test_auc<-100*(mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_auc)
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
            
            
          }
          
          #if(best_classifier_name==classifier_names[6])
          {
            
            if(tune_classifiers==TRUE){
              
              
              s1<-ldply(tune_nnet@tuneres,rbind)
              s2<-apply(s1,2,median)
              
              t1<-new("list")
              t1[[1]]<-s2
              
              tune_nnet1<-tune_nnet
              tune_nnet1@tuneres<-t1
              
              
              b1<-best(tune_nnet1)
              b2=unlist(unlist(b1))
              
              
              size_val=b1[[1]]$size
              decay_val=b1[[1]]$decay
              
              ##save(b1,file="b1.Rda")
              
              set.seed(seedvalue)
              class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE,tuneres=tune_nnet1) #size=3,decay=0.1) #,
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              
              set.seed(seedvalue)
              #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
              class_res2 <- nnetCMA(X =X, y = Y, learnind=as.numeric(learnmatrix),size=median(b2[seq(1,length(b2),2)]),decay=median(b2[seq(2,length(b2),2)])) #size=size_val,decay=decay_val)
              
            }else{
              
              set.seed(seedvalue)
              class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE) #size=3,decay=0.1) #,
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              
              set.seed(seedvalue)
              #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
              class_res2 <- nnetCMA(X =X, y = Y, learnind=as.numeric(learnmatrix)) #size=3,decay=0.1) #
            }
            
            ##saveclass_res2,file="nnetCMA.Rda")
            result_cma_list[[6]]<-class_res2
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_nnetCMA.Rda")
            confusion_matrix_list[[6]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ###savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              test_auc<-100*(mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_auc)
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
            
            
          }
          
          #if(best_classifier_name==classifier_names[7])
          {
            
            if(FALSE){
              s1<-ldply(tune_plr@tuneres,rbind)
              s2<-apply(s1,2,median)
              
              t1<-new("list")
              t1[[1]]<-s2
              
              tune_plr1<-tune_plr
              tune_plr1@tuneres<-t1
              
              
              b1<-best(tune_plr1)
            }
            
            set.seed(seedvalue)
            class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE) #,tuneres=tune_plr
            
            set.seed(seedvalue)
            class_res2 <- plrCMA(X =X, y = Y, learnind=learnmatrix)
            
            ##saveclass_res2,file="plrCMA.Rda")
            result_cma_list[[7]]<-class_res2
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_plrCMA.Rda")
            confusion_matrix_list[[7]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ###savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              test_auc<-100*(mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_auc)
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(testeval_res_auc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
            
            
            
            
          }
          ##save(result_cma_list,confusion_matrix_list,test_acc_mat,file="test_acc_mat.Rda")
          
          h2o_res<-NA
          
          #h2o
          if(deeplearning==TRUE){
            #      #save(X,Y,learnmatrix,fiveCV10iter,file="h2o.Rda")
            
            
            #loss="CrossEntropy",
            
            try(h2o.removeAll(),silent=TRUE)
            
            localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size='2g',nthreads=1)
            
            
            train_hex_split<-new("list")
            if(length(learnmatrix)<1000){
              set.seed(555)
              
              learnmatrixa<-learnmatrix[sample(1:length(learnmatrix),size=1000,replace=TRUE)]
            }else{
              learnmatrixa<-learnmatrix
              
            }
            
            stopping_metric="AUC" #logloss
            # stopping_metric="logloss"
            
            # learnmatrixa<-learnmatrix
            
            print("length of learnmatrix")
            print(length(learnmatrix))
            
            train_hex_split<-new("list")
            
            Class<-Y[learnmatrix]
            train_hex_split[[1]]<-cbind(X[learnmatrix,],Class)
            
            
            train_hex_split[[1]]<-as.data.frame(train_hex_split[[1]])
            train_hex_split[[1]]$Class<-as.factor(train_hex_split[[1]]$Class)
            
            
            if(FALSE)
            {
              learnmatrix_1<-learnmatrixa[1:(length(learnmatrixa)*0.7)]
              learnmatrix_2<-learnmatrixa[(length(learnmatrixa)*0.7+1):(length(learnmatrixa))]
              
              Class<-Y[learnmatrix_1]
              train_hex_split[[1]]<-cbind(X[learnmatrix_1,],Class)
              Class<-Y[learnmatrix_2]
              train_hex_split[[2]]<-cbind(X[learnmatrix_2,],Class)
              
              
              train_hex_split[[1]]<-as.data.frame(train_hex_split[[1]])
              train_hex_split[[1]]$Class<-as.factor(train_hex_split[[1]]$Class)
              train_hex_split[[2]]<-as.data.frame(train_hex_split[[2]])
              train_hex_split[[2]]$Class<-as.factor(train_hex_split[[2]]$Class)
            }
            
            hyper_params <- list(
              activation=c("Rectifier","Tanh","Maxout","RectifierWithDropout","TanhWithDropout","MaxoutWithDropout"),
              hidden=list(c(20,20),c(50,50),c(30,30,30),c(25,25,25,25),c(64,64,64),c(128,128),c(128,128,128),rep(128,5),c(100,100,100),c(16,32,64,128,256)),
              input_dropout_ratio=c(0,0.01,0.05),
              l1=seq(0,1e-2,1e-4),
              l2=seq(0,1e-2,1e-4),
              epochs=c(1,5,10),
              rho=c(0.9,0.95,0.99,0.999),
              epsilon=c(1e-10,1e-8,1e-6,1e-4)
            )
            
            
            ## Stop once the top 5 models are within 1% of each other (i.e., the windowed average varies less than 1%)
            search_criteria = list(strategy = "RandomDiscrete", max_runtime_secs = 360, max_models = 100, seed=1234567, stopping_rounds=2, stopping_tolerance=1e-2)
            
            #   set.seed(999)
            dl_random_grid <- h2o.grid(
              algorithm="deeplearning",
              nfolds=10,
              grid_id = "dl_grid_random",
              training_frame=as.h2o(train_hex_split[[1]]),
              validation_frame=as.h2o(train_hex_split[[1]]),
              x=1:(ncol(train_hex_split[[1]])-1),
              y="Class",
              stopping_metric=stopping_metric,reproducible=T,
              
              score_validation_samples=100, ## downsample validation set for faster scoring
              score_duty_cycle=0.025,         ## don't score more than 2.5% of the wall time
              max_w2=10,                      ## can help improve stability for Rectifier
              hyper_params = hyper_params,
              search_criteria = search_criteria
              
            )
            
            if(stopping_metric=="AUC"){ #
              grid <- h2o.getGrid("dl_grid_random",sort_by="AUC",decreasing=TRUE)
            }else{
              
              if(stopping_metric=="logloss"){
                #grid <- h2o.getGrid("dl_grid_random",sort_by="mean_per_class_error",decreasing=FALSE)
                
                grid <- h2o.getGrid("dl_grid_random",sort_by="logloss",decreasing=FALSE)
              }
            }
            #h2o.removeAll()
            
            #grid@summary_table[1,]
            best_model <- h2o.getModel(grid@model_ids[[1]]) ## model with lowest logloss
            
            
            
            train_hex_split<-new("list")
            
            Class<-Y[learnmatrix]
            train_hex_split[[1]]<-cbind(X[learnmatrix,],Class)
            Class<-Y[-learnmatrix]
            train_hex_split[[2]]<-cbind(X[-learnmatrix,],Class)
            
            
            train_hex_split[[1]]<-as.data.frame(train_hex_split[[1]])
            train_hex_split[[1]]$Class<-as.factor(train_hex_split[[1]]$Class)
            train_hex_split[[2]]<-as.data.frame(train_hex_split[[2]])
            train_hex_split[[2]]$Class<-as.factor(train_hex_split[[2]]$Class)
            
            
            amd.dl <- h2o.deeplearning(x=1:(ncol(train_hex_split[[1]])-1),y = "Class", seed=1234567, training_frame = as.h2o(train_hex_split[[1]]), activation=best_model@parameters$activation, epoch=best_model@parameters$epochs, l1=best_model@parameters$l1,l2=best_model@parameters$l2,rho = best_model@parameters$rho,epsilon = best_model@parameters$epsilon,
                                       nfolds=10,balance_classes=TRUE,hidden = best_model@parameters$hidden,input_dropout_ratio=best_model@parameters$input_dropout_ratio,
                                       stopping_metric=stopping_metric,validation_frame=as.h2o(train_hex_split[[2]][,1:(ncol(train_hex_split[[2]]))]),reproducible=T)
            
            #  amd.dl2 <- h2o.deeplearning(x=1:(ncol(train_hex_split[[1]])-1),y = "Class", seed=1234567, training_frame = as.h2o(train_hex_split[[1]]))
            #   amd.dl3 <- h2o.deeplearning(x=1:(ncol(train_hex_split[[1]])-1),y = "Class", seed=1234567, training_frame = as.h2o(train_hex_split[[1]]))
            
            
            predictions <- h2o.predict(amd.dl, as.h2o(train_hex_split[[2]][,1:(ncol(train_hex_split[[2]])-1)]))
            
            p1<-(predictions$predict)
            #  t1<-table(p1,as.vector(train_hex_split[[2]]$Class))
            
            predfit<-as.numeric(as.character(p1))
            
            if(FALSE){
              pred_acc<-multiclass.roc(testclass,as.numeric(predfit),levels=levels(as.factor(testclass)))
              
              pred_acc<-round(pred_acc$auc[1],2)
              
              test_auc<-100*(pred_acc)
              
              test_acc_mat<-c(test_acc_mat,test_auc)
            }
            
            
            h2o_res<-h2o.auc(amd.dl,train=TRUE,valid=TRUE,xval=TRUE)
            
            
            h2o_results<-list("res"=amd.dl,"pred"=predfit)
            
            ##save(amd.dl,h2o_res,file="h2o_res.Rda")
            
            if(is.na(test_data_check)==TRUE){
              h2o_res<-c(100*h2o_res[3],NA,NA,100*h2o_res[1])
              
            }else{
              h2o_res<-c(100*h2o_res[3],NA,NA,100*h2o_res[2])
            }
            
            # pred1 <- ROCR::prediction(as.vector(predfit), (train_hex_split[[2]]$Class))
            #p1<-performance(pred1,"auc")
            #p1=round(p1@y.values[[1]],2)
            
            ###save(amd.dl,predictions,train_hex_split,file="h2ores.Rda")
            
            try(h2o.shutdown(prompt=F),silent=FALSE)
            
          }
          
          
          
          if(length(class_levels_vec)==2){
            
            
            
            acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            acc_mat<-rbind(acc_mat,h2o_res)
            
            
            if(is.na(test_data_check)==TRUE){
              
              colnames(acc_mat)<-c("Training kfold CV Accuracy (AUC)","Training permuted kfold CV accuracy","Score", "Training accuracy(AUC)")
              
              print("Training set evaluation using selected features and the best classifier based on AUC measure")
            }else{ 
              colnames(acc_mat)<-c("Training kfold CV Accuracy (AUC)","Training permuted kfold CV accuracy","Score", "Test accuracy(AUC)")
              
              print("Test set evaluation using selected features and the best classifier based on AUC measure")
            }
            
            
            
            
            #rownames(acc_mat)<-best_classifier_name[1]
            
            write.table(acc_mat,file="Classification_evaluation_results_AUC.txt",sep="\t")
            
            
          }else{
            
            acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            
            if(is.na(test_data_check)==TRUE){
              
              colnames(acc_mat)<-c("Training kfold CV Accuracy (Misclassification)","Training permuted kfold CV Accuracy (Misclassification)","Score", "Training accuracy (Misclassification)")
              print("Training set evaluation using selected features and the best classifier based on misclassification rate measure")
              
            }else{             	
              colnames(acc_mat)<-c("Training kfold CV Accuracy (Misclassification)","Training permuted kfold CV Accuracy (Misclassification)","Score", "Test accuracy (Misclassification)")
              print("Test set evaluation using selected features and the best classifier based on misclassification rate measure")
            }
            
            print(acc_mat[,c(1,4)])
            
            #rownames(acc_mat)<-best_classifier_name[1]
            
            write.table(acc_mat,file="Classification_evaluation_results_misclassification.txt",sep="\t")
            
          }
          
          ###save(acc_mat,file="acc_mat.Rda")
          
          
          acc_mat<-acc_mat[,-c(3)]
          
          
          mainlab<-paste("Performance evaluation using classifiers and selected features",sep="")
          
          if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/Barplot_classification_accuracy.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
          }
          
          
          
          acc_mat1<-t(acc_mat)
          #xaxt="n",
          #
          w <- 0.1
          par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
          if(length(class_levels_vec)==2){
            barplot(acc_mat1[c(1,3),],main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
          }else{
            
            barplot(acc_mat1[c(1,3),],main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
            
          }
          if(FALSE){
            if(length(class_levels_vec)==2){
              axis(side=1,at=seq(1,4),labels=c("kfold CV","Permuted kfold CV","Test set","AUC"),cex.axis=cex.plots)
            }else{
              axis(side=1,at=seq(1,3),labels=c("kfold CV","Permuted kfold CV","Test set"),cex.axis=cex.plots)
            }
          }
          
          
          
          if(is.na(test_data_check)==TRUE){
            le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("Training (kfold)","Training (overall)"), pch = rep(19,2), pt.cex = 0.6, title = "",cex=0.7,col=barplot.col.opt))
            
          }else{
            le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("Training (kfold)","Test"), pch = rep(19,2), pt.cex = 0.6, title = "",cex=0.7,col=barplot.col.opt))
          }
          
          if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
          }
          
          
          try(dev.off(),silent=TRUE)
          
        }else{
          
          test_acc_mat<-c(0,0)
          
          if(tune_classifiers==TRUE){
            
            s1<-ldply(tune_scda@tuneres,rbind)
            s2<-apply(s1,2,median)
            
            t1<-new("list")
            t1[[1]]<-s2
            
            tune_scda1<-tune_scda
            tune_scda1@tuneres<-t1
            
            #tune_scda<-new("tuningresult",tuneres=t1)
            
            set.seed(seedvalue)
            class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, tuneres = tune_scda1,trace=FALSE)
            #scdaCMA(X, y, f, learnind, delta = 0.5, models=FALSE,...)
            
            b1<-best(tune_scda)
            
            
            learnmatrix<-as.numeric(learnmatrix)
            set.seed(seedvalue)
            class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix, delta = median(unlist(b1)))
            
          }else{
            
            if(is.na(class_scda)==FALSE){
              set.seed(seedvalue)
              class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA,trace=FALSE)
              set.seed(seedvalue)
              class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix)
              
            }
          }
          
          if(is.na(class_scda)==FALSE){
            
            ##saveclass_res2,file="scdaCMA.Rda")
            result_cma_list[[3]]<-class_res2
            #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_scdaCMA.Rda")
            confusion_matrix_list[[3]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ###savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              
              test_auc<-100*(mean(testeval_res_auc@score))
              
              test_acc_mat<-c(test_acc_mat,test_auc)
              
              
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(test_acc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
          }else{
            result_cma_list[[3]]<-NA
            confusion_matrix_list[[3]]<-NA
            test_acc_mat<-c(test_acc_mat,0)
            
          }
          
          
          
          #if(best_classifier_name==classifier_names[4])
          {
            
            if(tune_classifiers==TRUE){
              s1<-ldply(tune_svm@tuneres,rbind)
              
              s2<-apply(s1,2,median)
              
              t1<-new("list")
              t1[[1]]<-s2
              
              
              tune_svm1<-tune_svm
              tune_svm1@tuneres<-t1
              
              b1<-best(tune_svm1)
              
              
              #class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
              
              
              #                ##save(X,Y,learnmatrix,seedvalue,fiveCV10iter,tune_svm1,b1,svmCMA,file="svmdebug.Rda")
              set.seed(seedvalue)
              class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
              
              
              cost_1=b1[[1]]$cost
              gamma_1=b1[[1]]$gamma
              
              b2=unlist(unlist(b1))
              
              learnmatrix<-as.numeric(learnmatrix)
              
              
              
              
              set.seed(seedvalue)
              #class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix, cost = cost_1,gamma=gamma_1,probability=TRUE)
              class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix,probability=TRUE,cost=median(b2[seq(1,length(b2),2)]),gamma=median(b2[seq(2,length(b2),2)]))
              
            }else{
              
              set.seed(seedvalue)
              class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, kernel ="radial",trace=FALSE,probability=TRUE)
              set.seed(seedvalue)
              class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix,probability=TRUE)
              
            }
            
            ##saveclass_res2,file="svmCMA.Rda")
            result_cma_list[[4]]<-class_res2
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_svmCMA.Rda")
            confusion_matrix_list[[4]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ####savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              test_auc<-100*(mean(testeval_res_auc@score))
              
              test_acc_mat<-c(test_acc_mat,test_auc)
              
              
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(test_acc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
            
          }
          
          #if(best_classifier_name==classifier_names[5])
          {
            
            if(FALSE){
              s1<-ldply(tune_rf@tuneres,rbind)
              s2<-apply(s1,2,median)
              
              t1<-new("list")
              t1[[1]]<-s2
              
              tune_rf1<-tune_rf
              tune_rf1@tuneres<-t1
            }
            
            set.seed(seedvalue)
            class_res <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE) #,tuneres=tune_rf1
            
            learnmatrix<-as.numeric(learnmatrix)
            
            
            set.seed(seedvalue)
            class_res2 <- rfCMA(X =X, y = Y, learnind=learnmatrix, varimp = FALSE) #mtry=tune_rf1$mtry,nodesize=tune_rf1$nodesize
            
            ##saveclass_res2,file="rfCMA.Rda")
            result_cma_list[[5]]<-class_res2
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
            confusion_matrix_list[[5]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ####savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              test_auc<-100*(mean(testeval_res_auc@score))
              
              test_acc_mat<-c(test_acc_mat,test_auc)
              
              
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(test_acc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
            
          }
          
          #if(best_classifier_name==classifier_names[6])
          {
            
            if(tune_classifiers==TRUE){
              s1<-ldply(tune_nnet@tuneres,rbind)
              s2<-apply(s1,2,median)
              
              t1<-new("list")
              t1[[1]]<-s2
              
              tune_nnet1<-tune_nnet
              tune_nnet1@tuneres<-t1
              
              
              b1<-best(tune_nnet1)
              b2=unlist(unlist(b1))
              
              set.seed(seedvalue)
              class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE,tuneres=tune_nnet1) #size=3,decay=0.1) #
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              
              set.seed(seedvalue)
              #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
              class_res2 <- nnetCMA(X =X, y = Y, learnind=as.numeric(learnmatrix),size=median(b2[seq(1,length(b2),2)]),decay=median(b2[seq(2,length(b2),2)])) #size=3,decay=0.1)
              
            }else{
              
              set.seed(seedvalue)
              class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE) #size=3,decay=0.1) #
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              
              set.seed(seedvalue)
              #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
              class_res2 <- nnetCMA(X =X, y = Y, learnind=as.numeric(learnmatrix)) #size=3,decay=0.1)
              
              
            }
            
            ##saveclass_res2,file="nnetCMA.Rda")
            result_cma_list[[6]]<-class_res2
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_nnetCMA.Rda")
            confusion_matrix_list[[6]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ####savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              test_auc<-100*(mean(testeval_res_auc@score))
              
              test_acc_mat<-c(test_acc_mat,test_auc)
              
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(test_acc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
            
          }
          
          #if(best_classifier_name==classifier_names[7])
          {
            
            if(FALSE){
              s1<-ldply(tune_plr@tuneres,rbind)
              s2<-apply(s1,2,median)
              
              t1<-new("list")
              t1[[1]]<-s2
              
              tune_plr1<-tune_plr
              tune_plr1@tuneres<-t1
              
              
              b1<-best(tune_plr1)
            }
            
            set.seed(seedvalue)
            class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE) #,tuneres=tune_plr
            
            set.seed(seedvalue)
            class_res2 <- plrCMA(X =X, y = Y, learnind=learnmatrix)
            
            ##saveclass_res2,file="plrCMA.Rda")
            result_cma_list[[7]]<-class_res2
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            ###saveconfusion_matrix_1,file="confusion_matrix_plrCMA.Rda")
            confusion_matrix_list[[7]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
              
              testeval_res_auc<-evaluation(class_res,measure = "auc")
              ####savetesteval_res_auc,file="testeval_res_auc.Rda")
              
              test_auc<-100*(mean(testeval_res_auc@score))
              
              test_acc_mat<-c(test_acc_mat,test_auc)
              
              
              
            }else{
              
              test_acc<-evaluation(class_res)
              
              test_acc<-100*(1-mean(test_acc@score))
              test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
              
              
            }
            
          }
          
          
          if(length(class_levels_vec)==2){
            
            
            
            acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            
            if(is.na(test_data_check)==TRUE){     
              colnames(acc_mat)<-c("Training kfold CV Accuracy (AUC)","Training permuted kfold CV accuracy","Score", "Training accuracy(AUC)")
              
              print("Training set evaluation using selected features and the best classifier based on AUC measure")
            }else{
              colnames(acc_mat)<-c("Training kfold CV Accuracy (AUC)","Training permuted kfold CV accuracy","Score", "Test accuracy(AUC)")
              
              print("Test set evaluation using selected features and the best classifier based on AUC measure")
              
            }
            
            
            ###save(acc_mat,file="acc_mat.Rda")
            #    print(acc_mat[,c(1,4)])
            
            #rownames(acc_mat)<-best_classifier_name[1]
            
            write.table(acc_mat,file="Classification_evaluation_results_AUC.txt",sep="\t")
            
            
          }else{
            
            acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            if(is.na(test_data_check)==TRUE){ 
              colnames(acc_mat)<-c("Training kfold CV Accuracy (Misclassification)","Training permuted kfold CV Accuracy (Misclassification)","Score", "Training accuracy (Misclassification)")
              
              print("Training set evaluation using selected features and the best classifier based on misclassification rate measure")
            }else{
              
              colnames(acc_mat)<-c("Training kfold CV Accuracy (Misclassification)","Training permuted kfold CV Accuracy (Misclassification)","Score", "Test accuracy (Misclassification)")
              
              print("Test set evaluation using selected features and the best classifier based on misclassification rate measure")
              
              
            }
            
            print(acc_mat[,c(1,4)])
            
            #rownames(acc_mat)<-best_classifier_name[1]
            
            write.table(acc_mat,file="Classification_evaluation_results_misclassification.txt",sep="\t")
            
          }
          
          
          
          acc_mat<-acc_mat[,-c(3)]
          
          
          mainlab<-paste("Performance evaluation using classifiers and selected features",sep="")
          
          if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/Barplot_classification_accuracy.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
          }
          
          
          acc_mat1<-t(acc_mat)
          #xaxt="n",
          #
          w <- 0.1
          par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
          if(length(class_levels_vec)==2){
            barplot(acc_mat1[c(1,3),],main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
          }else{
            
            barplot(acc_mat1[c(1,3),],main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
            
          }
          if(FALSE){
            if(length(class_levels_vec)==2){
              axis(side=1,at=seq(1,4),labels=c("kfold CV","Permuted kfold CV","Test set","AUC"),cex.axis=cex.plots)
            }else{
              axis(side=1,at=seq(1,3),labels=c("kfold CV","Permuted kfold CV","Test set"),cex.axis=cex.plots)
            }
          }
          #col = col_vec[1:length(t1)],
          le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("Training (kfold)","Test"), pch = rep(19,2), pt.cex = 0.6, title = "",cex=0.7,col=barplot.col.opt))
          
          
          if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
          }
          
          
          try(dev.off(),silent=TRUE)
          print("number of features too small.")
        }
        
      }else{
        
        
        acc_mat<-acc_mat[,-c(3)]
        
        acc_mat1<-t(acc_mat)
        
        colnames(acc_mat)<-c("kfold CV accuracy","Permuted kfold CV accuracy")
        mainlab<-paste("Performance evaluation using ",best_classifier_name," classifier and selected features",sep="")
        
        if(output.device.type!="pdf"){
          
          temp_filename_1<-"Figures/Barplot_classification_accuracy.png"
          
          png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        
        w <- 0.1
        par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
        if(length(class_levels_vec)==2){
          barplot(acc_mat1,main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
        }else{
          
          barplot(acc_mat1,main=mainlab,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE,beside=TRUE,col=barplot.col.opt,cex.axis=0.7,cex.names=0.7)
          
        }
        
        # barplot(acc_mat,xaxt="n",main=mainlab,col=barplot.col.opt,ylab="Classification accuracy(%)",ylim=c(50,100),xpd=FALSE)
        # axis(side=1,at=seq(1,2),labels=c("Best kfold CV","Permuted kfold CV"),cex.axis=cex.plots)
        
        le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("Training (kfold)","Training (permuted)"), pch = rep(19,2), pt.cex = 0.6, title = "",cex=0.7,col=barplot.col.opt))
        
        
        if(output.device.type!="pdf"){
          
          try(dev.off(),silent=TRUE)
        }
        
        
      }
      
      
      
      
    }
    
    
    try(dev.off(),silent=TRUE)
    
    eval_mat1<-rbind(eval_mat1,eval_mat_perm_final)
    
    
    
    
    #emat1
    write.table(emat1,file="Classifier_accuracy_comparison_training_set.txt",sep="\t",row.names=TRUE)
    
    confounder_eval=NA
    if(is.na(confounder.matrix)==FALSE){
      
      response_mat<-confounder.matrix[,-c(1)]
      confounder_eval<-try(runlmreg(X=good_feats,Y=response_mat,fdrmethod=confounderfdrmethod,fdrthresh=confounderfdrthresh),silent=TRUE)
    }
    
    parentoutput_dir=outloc
    setwd(parentoutput_dir)
    
    if(globalcor==TRUE){
      
      print("##############Level 2: Correlation network analysis selected features###########")
      print(paste("Generating metabolome-wide ",cor.method," correlation network",sep=""))
      data_m_fc_withfeats<-as.data.frame(Xorig)
      
      good_feats<-as.data.frame(good_feats)
      #print(goodfeats[1:4,])
      sigfeats_index<-which(data_m_fc_withfeats$mz%in%good_feats$mz)
      sigfeats<-sigfeats_index
      
      #outloc<-paste(parentoutput_dir,"/Allcornetworksigfeats","log2fcthresh",log2.fold.change.thresh,"/",sep="")
      outloc<-paste(parentoutput_dir,"/MWASresults","/",sep="")
      
      suppressWarnings(dir.create(outloc,showWarnings = FALSE))
      setwd(outloc)
      
      #cor.method="spearman",networktype="complete",abs.cor.thresh=0.4,cor.fdrthresh=0.05,max.cor.num=100,net_node_colors=c("green","red"), net_legend=TRUE
      
      if(networktype=="complete"){
        mwan_fdr<-do_cor(data_m_fc_withfeats,subindex=sigfeats_index,targetindex=NA,outloc,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,max.cor.num,net_node_colors,net_legend,cex.plots=cex.plots)
      }else{
        if(networktype=="GGM"){
          mwan_fdr<-get_partial_cornet(data_m_fc_withfeats, sigfeats.index=sigfeats_index,targeted.index=NA,networkscope="global",cor.method,
                                       abs.cor.thresh,cor.fdrthresh,outloc=outloc,net_node_colors,net_legend)
        }else{
          print("Invalid option. Please use complete or GGM.")
        }
      }
      
      print("##############Level 2: processing complete###########")
    }
    
    if(length(good_feats)>0){
      
      write.csv(good_feats,file="train_selected_data.csv",row.names=FALSE)
      try(write.csv(good_feats_test,file="test_selected_data.csv",row.names=FALSE),silent=TRUE)
      try(write.csv(confounder_eval,file="confounder_eval.csv",row.names=FALSE),silent=TRUE)
      
      write.csv(Ytrain_mat,file="train_classlabels_data.csv",row.names=FALSE)
      try(write.csv(Ytest_mat,file="test_classlabels_data.csv",row.names=FALSE),silent=TRUE)
    }
    
    col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF","grey57")
    auc_vec<-{}
    classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
    
    
    
    if(length(class_levels)==2){
      png("ROC_classifier_comparison.png",width=8,height=8,res=300,type="cairo",units="in")
      if(length(good_feats_index)>=3){
        start_ind=1
      }else{
        start_ind=3
      }
      for(i in start_ind:length(result_cma_list))
      {
        
        
        object=result_cma_list[[i]]
        
        
        
        roc1<-try(ROCinternal.panda(test = object@prob[,2], object@y,FALSE),silent=TRUE)
        
        if(is(roc1,"try-error")){
          
        }else{
          
          
          auc_vec<-c(auc_vec,round(roc1$auc,2))
          roc1=roc1$plotcoordinates
          
          
          if(i==start_ind){
            plot(roc1$x,roc1$y,col=col_vec[i],type="s",lty=i,lwd=2,xlab="1-specificity",ylab="Sensitivity",main="Receiver Operating Characteristic curve")
            lines(seq(0,1),seq(0,1),lwd=2)
          }else{
            
            lines(roc1$x,roc1$y,col=col_vec[i],type="s",lty=i,lwd=2)
          }
        }
        
        
      }
      
      legend('bottomright', paste(classifier_names,":",auc_vec,sep=""), col=col_vec[1:7], lty=1:7, cex=0.8, bty='n',lwd=rep(2,7))
      
      dev.off()
      
    }
    
    options(warn=0)
    #print(acc_mat)
    #bestcvacc=best_kfold_acc,permcvacc=permkfold_acc,testacc=test_acc,
    return(list(train.selected.data=good_feats,test.selected.data=good_feats_test,classifier.comparison=emat1,feature.selection.matrix=stability_matrix_1,best.performance.measures=acc_mat,train.class=Ytrain_mat,test.class=Ytest_mat,confounder.eval.res=confounder_eval,test.result.classifiers=result_cma_list,h2o_results=h2o_results))
    
  }
  
  
  
}

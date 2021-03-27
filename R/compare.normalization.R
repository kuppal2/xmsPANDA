compare.normalization <-
function(Xmat=NA,Ymat=NA,Zmat=NA,feature_table_file=NA,parentoutput_dir,class_labels_file=NA,num_replicates=3,feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method="mean",
                                all.missing.thresh=0.5,group.missing.thresh=0.7,
                                missing.val=0,samplermindex=NA, rep.max.missing.thresh=0.5,summary.na.replacement="zeros",pairedanalysis=FALSE,
                                normalization.method=c("all","log2quantilenorm","log2transform","znormtransform","lowess_norm","quantile_norm","rangescaling",
                                                       "paretoscaling","mstus","eigenms_norm","vsn_norm","sva_norm","tic_norm","cubicspline_norm","mad_norm"),input.intensity.scale="raw",
                                abs.cor.thresh=0.4,pvalue.thresh=0.05,cor.fdrthresh=0.2,cex.plots=0.7,plots.width=8,plots.height=8,plots.res=600,
                                plots.type="cairo",heatmap.col.opt="RdBu",cor.method="pearson",pca.ellipse=FALSE,ground_truth_file=NA,cutree.method="dynamic",rsd.filt.thresh=1,alphabetical.order=TRUE,analysistype="classification",lme.modeltype="RI",
                                study.design=c("multiclass","onewayanova","twowayanova","onewayanovarepeat","twowayanovarepeat"),log2.transform.constant=1,
                                featselmethod="limma",similarity.matrix="correlation",...){
  
  #featselmethod=NA
  suppressMessages(require(dynamicTreeCut))
  suppressMessages(require(mclust))
  suppressMessages(require(cluster))
  suppressMessages(require(WGCNA))
  
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  if(is.na(Xmat)==TRUE){
    Xmat<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  }
  
  if(is.na(Ymat)==TRUE){
    Ymat<-read.table(class_labels_file,sep="\t",header=TRUE)
  }
  
  if(is.na(ground_truth_file)==FALSE){
    
    Zmat<-read.table(ground_truth_file,sep="\t",header=TRUE)
  }
  
  
  if(study.design=="multiclass" | study.design=="onewayanova"){
    
    
    print("Treating column B as class labels")
    cnames<-colnames(Ymat)
    cnames[2]<-c("Class")
    colnames(Ymat)<-cnames
  }else{
    if(study.design=="onewayanovarepeat"){
      
      
      print("Treating column B as subject IDs, column C as class labels")
      Ymat<-Ymat[,-c(2)]
      cnames<-colnames(Ymat)
      cnames[2]<-c("Class")
      colnames(Ymat)<-cnames
    }else{
      if(study.design=="twowayanovarepeat"){
        
        
        print("Treating column B as subject IDs, column C as factor 1, and column D as factor 2")
        Ymat<-Ymat[,-c(2)]
        Ymat[,2]<-as.factor(Ymat[,2]):as.factor(Ymat[,3])
        Ymat<-Ymat[,-3]
        cnames<-colnames(Ymat)
        cnames[2]<-c("Class")
        colnames(Ymat)<-cnames
      }else{
        
        if(study.design=="twowayanova"){
          
          
          print("Treating column B as factor 1 and column C as factor 2")
          
          Ymat[,2]<-as.factor(Ymat[,2]):as.factor(Ymat[,3])
          Ymat<-Ymat[,-3]
          cnames<-colnames(Ymat)
          cnames[2]<-c("Class")
          colnames(Ymat)<-cnames
        }
        
        
      }
      
    }
    
    
    
  }
  
  normalization.method=tolower(normalization.method)
  
  
  if(normalization.method[1]==c("all")){
    
    normalization.method=c("log2quantilenorm","log2transform","znormtransform","lowess_norm","quantile_norm","rangescaling","paretoscaling","mstus","eigenms_norm","vsn_norm","sva_norm","tic_norm","cubicspline_norm","mad_norm")
    
  }
  
  if(normalization.method[1]==c("none")){
    print("Using raw data")
  }
  valid_norm_options<-c("log2quantilenorm","log2transform","znormtransform","lowess_norm","quantile_norm","rangescaling","paretoscaling","mstus",
                        "eigenms_norm","vsn_norm","sva_norm","none","tic_norm","cubicspline_norm","mad_norm")
  
  data_matrix_list<-new("list")
  
  rnames_xmat<-colnames(Xmat)
  rnames_xmat<-tolower(rnames_xmat)
  #colnames(Xmat)<-rnames_xmat
  
  rnames_ymat<-colnames(Ymat)
  rnames_ymat<-tolower(rnames_ymat)
  
  check_ylabel<-regexpr(rnames_ymat,pattern="^((class)|(factor))",perl=TRUE)
  
  
  if(length(which(check_ylabel>0))<1){
    
    stop("No column named Class found in the class labels file. Please use the name Class for the column with the main outcome.")
  }
  
  
  colnames(Ymat)<-rnames_ymat
  
  #Xmat1<-Xmat[,c(1,order(rnames_xmat)+1)]
  
  
  
  check_xlabel<-regexpr(rnames_xmat[1],pattern="^(name|Name)",perl=TRUE)
  if(check_xlabel>0){
    
    
    # Xmat<-cbind(Xmat[,c(1)],Xmat[,c(order(Ymat$class)+1)])
    
    cnames1<-colnames(Xmat)
    cnames1[1]<-"Name"
    
    colnames(Xmat)<-cnames1
    
    #    Ymat<-Ymat[order(Ymat$class),]
    
    rnames_xmat<-colnames(Xmat[,-c(1)])
    rnames_ymat<-Ymat[,1]
    
    # rnames_ymat<-tolower(rnames_ymat)
    
  }else{
    
    
    
    
    #Xmat<-cbind(Xmat[,c(1:2)],Xmat[,c(order(Ymat$class)+2)])
    
    
    # Ymat<-Ymat[order(Ymat$class),]
    
    
    rnames_xmat<-colnames(Xmat[,-c(1:2)])
    rnames_ymat<-Ymat[,1]
    # rnames_ymat<-tolower(rnames_ymat)
  }
  
  
  data_matrix_list<-lapply(1:length(normalization.method),function(i)
  {
    
    if(normalization.method[i]%in%valid_norm_options){
      
      print("###################")
      
      print("          ")
      
      dir.create(parentoutput_dir,showWarnings = FALSE)
      setwd(parentoutput_dir)
      dir.create(normalization.method[i],showWarnings = FALSE)
      setwd(normalization.method[i])
      outloc=paste(parentoutput_dir,"/",normalization.method[i],sep="")
      
      print(outloc)
      fname1<-paste(normalization.method[i],".pdf",sep="")
      
      
      pdf(fname1)
      par(mfrow=c(1,1),family="sans",cex=cex.plots)
      
      ######
      
      if(is.na(featselmethod)==FALSE){
        
        if(featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="poissonreg" | featselmethod=="lmregrobust" | featselmethod=="logitregrobust" | featselmethod=="poissonregrobust" | featselmethod=="limma" | featselmethod=="limma1way" | featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls" | featselmethod=="rfesvm" | featselmethod=="pamr" | featselmethod=="poissonreg" | featselmethod=="ttest" | featselmethod=="wilcox" | featselmethod=="lm1wayanova"){
          
          
          Ymat<-Ymat[,c(1:2)]
          
        }
        
      }else{
        
        print("Treating the second column as the class label.")
        Ymat<-Ymat[,c(1:2)]
        
      }
      
      
      
      
      #rnames_ymat[match(rnames_xmat,rnames_ymat)]
      
      if(length(which(duplicated(rnames_ymat)==TRUE))>0){
        
        stop("Duplicate sample IDs are not allowed. Please represent replicates by _1,_2,_3.")
      }
      
      check_ylabel<-regexpr(rnames_ymat[1],pattern="^[0-9]*",perl=TRUE)
      check_xlabel<-regexpr(rnames_xmat[1],pattern="^X[0-9]*",perl=TRUE)
      
      check_ylabel<-regexpr(rnames_ymat[1],pattern="^ ",perl=TRUE)
      check_xlabel<-regexpr(rnames_xmat[1],pattern="^ ",perl=TRUE)
      if(length(check_ylabel)>0 && length(check_xlabel)>0){
        if(attr(check_ylabel,"match.length")>0 && attr(check_xlabel,"match.length")>0){
          
          rnames_ymat<-paste("X",rnames_ymat,sep="")
        }
      }
      
      match_names<-match(rnames_xmat,rnames_ymat)
      
      bad_colnames<-length(which(is.na(match_names)==TRUE))
      
      if(bad_colnames>0){
        print("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names.")
        print("Sample names in feature table")
        print(head(rnames_xmat))
        print("Sample names in classlabels file")
        
        print(head(rnames_ymat))
        #stop("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names. Please try again.")
      }
      
      
      if(is.na(all(diff(match(rnames_xmat,rnames_ymat))))==FALSE){
        if(all(diff(match(rnames_xmat,rnames_ymat)) > 0)==FALSE){
         
          print(head(rnames_xmat))
          print("Sample names in classlabels file")
          
          print(head(rnames_ymat))
          
           stop("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names. Please try again.")
        }
      }else{
        print(head(rnames_xmat))
        print("Sample names in classlabels file")
        
        print(head(rnames_ymat))
        stop("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names. Please try again.")
      }
      
      data_matrix_list<-data_preprocess(Xmat=Xmat,Ymat=Ymat,feature_table_file=feature_table_file,parentoutput_dir=outloc,class_labels_file=class_labels_file,num_replicates=num_replicates,feat.filt.thresh=NA,summarize.replicates=summarize.replicates,
                                        summary.method=summary.method,
                                        all.missing.thresh=all.missing.thresh,group.missing.thresh=group.missing.thresh,log2transform=FALSE,quantile_norm=FALSE,missing.val=missing.val,samplermindex=samplermindex, rep.max.missing.thresh=rep.max.missing.thresh,
                                        summary.na.replacement=summary.na.replacement,
                                        featselmethod=featselmethod,pairedanalysis=pairedanalysis,normalization.method=normalization.method[i],input.intensity.scale=input.intensity.scale,log2.transform.constant=log2.transform.constant) #,silent=TRUE)
      
      if(is(data_matrix_list,"try-error")){
        
        return(list("data_matrix_list"=NA,"Silhouette.sample"=NA,"Silhouette.feature"=NA,"ari_val"=NA,"feat.cor.mat"=NA))
      }else{
        
        
        feature_names<-data_matrix_list$names_with_mz_time
        
        X<-data_matrix_list$data_matrix_afternorm_scaling
        
        sd_check<-apply(X[,-c(1:2)],1,sd)
        
        if(length(which(sd_check==0))>0){
          
          X<-X[-which(sd_check==0),]
        }
        
        feat_rsds<-apply(X[,-c(1:2)],1,do_rsd)
        
        sum_rsd<-summary(feat_rsds,na.rm=TRUE)
        max_rsd<-max(feat_rsds,na.rm=TRUE)
        max_rsd<-round(max_rsd,2)
        
       # print("Summary of RSD across all features:")
        #print(sum_rsd)
        
        abs_feat_rsds<-abs(feat_rsds)
        
        good_metabs<-which(abs_feat_rsds>rsd.filt.thresh)
        
        if(length(good_metabs)>0){
          
          X<-X[good_metabs,]
          
        }else{
          
          
          stop(paste("Please decrease the maximum relative standard deviation (rsd.filt.thresh) threshold to ",max_rsd,sep=""))
          
        }
        
      
        setwd(parentoutput_dir)
        
        setwd(normalization.method[i])
        
        if(is.na(feature_names)==FALSE){
          
          feature_names<-merge(feature_names,X[,c(1:2)],by=c("mz","time"))
          feature_names<-feature_names[order(as.numeric(as.character(feature_names$mz)),feature_names$time),]
          
          dup_name_check<-which(duplicated(feature_names$Name)==TRUE)
          if(length(dup_name_check)>0){
            feature_names<-feature_names[-dup_name_check,]
            
            X<-X[-dup_name_check,]
          }
          rownames(X)<-feature_names$Name
          
        }
        Y<-data_matrix_list$classlabels
        
        size_num<-min(100,dim(X[,-c(1:2)])[2])
        
        data_m_fc<-X[,-c(1:2)]
        
        par(mfrow=c(1,1),family="sans",cex=cex.plots)
        samp_index<-sample(x=1:dim(data_m_fc)[2],size=size_num)
        samp_index<-samp_index[order(samp_index)]
        
        if(normalization.method[i]=="none"){
          
          boxplot(data_m_fc[,samp_index],main=paste("Intensity distribution across samples after no normalization",sep=""),xlab="Samples",ylab="Intensity",col="white",cex=0.8)
          
        }else{
          boxplot(data_m_fc[,samp_index],main=paste("Intensity distribution across samples after ", normalization.method[i],sep=""),xlab="Samples",ylab="Intensity",col="white",cex=0.8)
        }
        
        ##save(Y,file="Y.Rda")
        if(pairedanalysis==TRUE){
          
          Y1=Y[,-c(2)]
        }
        
       # print("Starting PCA")
        
        #  get_pcascoredistplots(X=X,Y=Y,feature_table_file=NA,parentoutput_dir=outloc,class_labels_file=NA,sample.col.opt="journal",plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, alphacol=0.3,col_vec=NA,pairedanalysis=pairedanalysis,
        #pca.cex.val=2,legendlocation="topright",pca.ellipse=pca.ellipse,ellipse.conf.level=0.95,filename=paste(normalization.method[i], " data ",sep=""),paireddesign=NA,error.bar=TRUE,
        #lineplot.col.opt="black",lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),newdevice=FALSE,timeseries.lineplots=FALSE,pcascale=TRUE,pcacenter=TRUE)
        
        
        #print(X[1:3,1:5])
        #print(head(Y))
       # save(X,Y,pairedanalysis,pca.ellipse,normalization.method,i,alphabetical.order,analysistype,lme.modeltype,file="pcad1.Rda")
      
        cat("Starting PCA",sep="\n")  
        pcares<-get_pcascoredistplots(X=X,Y=Y,feature_table_file=NA,parentoutput_dir=outloc,class_labels_file=NA,sample.col.opt="journal",
                                      plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=NA,pairedanalysis=pairedanalysis,pca.cex.val=2,legendlocation="topright",
                                      pca.ellipse=pca.ellipse,ellipse.conf.level=0.95,filename=paste(normalization.method[i], "data ",sep=""),paireddesign=NA,
                                      lineplot.col.opt="black",lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),timeseries.lineplots=FALSE,pcacenter=TRUE,pcascale=TRUE,alphabetical.order=alphabetical.order,study.design=analysistype,lme.modeltype=lme.modeltype) #,silent=TRUE)
        
       #save(pcares,file="pcares.Rda")
      
        setwd(parentoutput_dir)
        
        classgroup=Y[,2]
        classgroup=as.data.frame(classgroup)
        
        
        setwd(normalization.method[i])
        cat("Starting HCA")
        st1=Sys.time()
        
        #cluster features
        x=pcares$X
        #row samples; column: features
        dist1<-as.dist(1-WGCNA::cor(x))
        f1=fastcluster::hclust(d=dist1,method = "complete")
        
        mycl_metabs <-cutreeDynamic(f1,distM=as.matrix((1-WGCNA::cor(x))),method="hybrid",cutHeight = 0.95,
                                    deepSplit = 4,minClusterSize = 2,verbose = FALSE)
        names(mycl_metabs)<-f1$labels
        
        Silhouette.features<- silhouette(mycl_metabs,dmatrix=as.matrix(dist1))
        
        Silhouette.features<-round(mean(Silhouette.features[,3]),3)
        
        #transpose for W.k: row: features; col: samples
        #gap.stat.features<-E.W.k(x=(x), d.power=1,B=100)-W.k(x=t(x), clus=mycl_metabs,d.power=1)
        
        
        #cluster samples
        x=t(pcares$X)
        #row samples; column: features
        dist1<-as.dist(1-WGCNA::cor(x))
        f1=fastcluster::hclust(d=dist1,method = "complete")
        
        mycl_samples <-cutreeDynamic(f1,distM=as.matrix((1-WGCNA::cor(x))),method="hybrid",cutHeight = 0.95,
                                    deepSplit = 2,minClusterSize =2,verbose = FALSE)
        names(mycl_samples)<-f1$labels
        
        Silhouette.samples<- silhouette(mycl_samples,dmatrix=as.matrix(dist1))
        Silhouette.samples<-round(mean(Silhouette.samples[,3]),3)
        
        ari_val<-try(round(adjustedRandIndex(x=classgroup[,1], y=mycl_samples),2),silent=TRUE)
        
       # save(f1,mycl_samples,classgroup,file="hcasamples.Rda")
        
        plotDendroAndColors(f1,mycl_samples,rowText = mycl_samples,groupLabels = c("Cluster"),
                            rowTextAlignment = "center",
                            dendroLabels = rownames(f1),
                            main=paste("Cluster dendrogram for samples\n Adjusted Rand Index:",ari_val,sep=""))
        
        
       
        #transpose for W.k: row: features; col: samples
    #    gap.stat.samples<-E.W.k(x=(x), d.power=1,B=100)-W.k(x=t(x), clus=mycl_samples,d.power=1)
        
        
        #W.k <- function(x, kk,clus,d.power=1) {
      if(FALSE)
      {
        hca_res<-get_hca(parentoutput_dir=outloc,X=X,Y=Y,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                         sample.col.opt="journal",plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, plots.type=plots.type, 
                         alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="",
                         plot.bycluster=FALSE,color.rows=FALSE,deepsplit=2,minclustsize=1,mergeCutHeight=0.1,
                         num_nodes=2,alphabetical.order=alphabetical.order,pairedanalysis=pairedanalysis,cutree.method=cutree.method,
                         study.design=analysistype,similarity.matrix=similarity.matrix,cexRow=cex.plots,cexCol=cex.plots)
        
        classlabels_response_mat=Y[,2]
        st2=Sys.time()
        cat("Done with HCA",sep="\n")
        print(st2-st1)
      }
        
        
        if(FALSE){ 
        if(analysistype=="classification"){
          
          #print("Performing one-way ANOVA analysis")
          
          #numcores<-round(detectCores()*0.6)
          cl <- parallel::makeCluster(getOption("cl.cores", 2))
          clusterExport(cl,"diffexponewayanova",envir = .GlobalEnv)
          clusterExport(cl,"anova",envir = .GlobalEnv)
          clusterExport(cl,"TukeyHSD",envir = .GlobalEnv)
          clusterExport(cl,"aov",envir = .GlobalEnv)
          res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
            #res1<-apply(data_m_fc,1,function(x,classlabels_response_mat){
            xvec<-x
            data_mat_anova<-cbind(xvec,classlabels_response_mat)
            data_mat_anova<-as.data.frame(data_mat_anova)
            cnames<-colnames(data_mat_anova)
            cnames[1]<-"Response"
            colnames(data_mat_anova)<-c("Response","Factor1")
            data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
            #anova_res<-diffexponewayanova(dataA=data_mat_anova)
            data_mat_anova<-as.data.frame(data_mat_anova)
            
            ##save(data_mat_anova,file="data_mat_anova.Rda")
            a1 <- aov(as.numeric(data_mat_anova$Response) ~ .,data=data_mat_anova) # + chocolate$Factor1*chocolate$Factor2)
            
            # posthoc <- TukeyHSD(x=a1, conf.level=0.95,test=univariate())
            anova_res<-anova(a1)
            
            num_rows<-dim(anova_res)
            pvalues_factors<-data.frame(t(anova_res["Pr(>F)"][-c(num_rows),]))
            
            
            
            return(pvalues_factors)
          },classlabels_response_mat)
          
          stopCluster(cl)
          
          pvalues<-{}
          
          
          pvalues<-unlist(res1)
          
        }
        else{
          
          # numcores<-num_nodes #round(detectCores()*0.5)
          
          cl <- parallel::makeCluster(getOption("cl.cores", 2))
          
          clusterExport(cl,"diffexplmreg",envir = .GlobalEnv)
          clusterExport(cl,"lm",envir = .GlobalEnv)
          clusterExport(cl,"glm",envir = .GlobalEnv)
          clusterExport(cl,"summary",envir = .GlobalEnv)
          clusterExport(cl,"anova",envir = .GlobalEnv)
          clusterEvalQ(cl,library(sandwich))
          
          
          #data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
          res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,logistic_reg,poisson_reg,robust.estimate){
            
            xvec<-x
            
            
            
            data_mat_anova<-cbind(xvec,classlabels_response_mat)
            
            cnames<-colnames(data_mat_anova)
            cnames[1]<-"Response"
            
            colnames(data_mat_anova)<-cnames
            
            
            anova_res<-diffexplmreg(dataA=data_mat_anova,logistic_reg,poisson_reg,robust.estimate)
            
            return(anova_res)
          },classlabels_response_mat,FALSE,FALSE,TRUE)
          
          stopCluster(cl)
          main_pval_mat<-{}
          
          posthoc_pval_mat<-{}
          pvalues<-{}
          
          all_inf_mat<-{}
          
          for(i in 1:length(res1)){
            
            main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
            pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
          }
          pvalues<-unlist(pvalues)
          
          #posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthocfactor1)
          
          
        }
        
        
        kstest_res<-ks.test(pvalues,"punif",0,1)
        kstest_res<-round(kstest_res$p.value,3)
        
        hist(as.numeric(pvalues),main=paste("Distribution of pvalues\n","Kolmogorov-Smirnov test for uniform distribution, p=",kstest_res,sep=""),cex.main=0.75,xlab="pvalues")
        #if(FALSE)
        {
          simpleQQPlot = function (observedPValues,mainlab) {
            plot(-log10(1:length(observedPValues)/length(observedPValues)),
                 -log10(sort(observedPValues)),main=mainlab,xlab=paste("Expected -log10pvalue",sep=""),ylab=paste("Observed -logpvalue",sep=""),cex.main=0.75)
            abline(0, 1, col = "brown")
          }
          
          inflation <- function(pvalue) {
            chisq <- qchisq(1 - pvalue, 1)
            lambda <- median(chisq) / qchisq(0.5, 1)
            lambda
          }
          inflation_res<-round(inflation(pvalues),2)
          
          simpleQQPlot(pvalues,mainlab=paste("QQplot pvalues","\np-value inflation factor: ",inflation_res," (no inflation: close to 1; bias: greater than 1)",sep=""))
          
        }
        
        }
        
        
        
        if(FALSE){
        
        Silhouette.sample=hca_res$Silhouette.sample
        
        Silhouette.feature=hca_res$Silhouette.feature
        
        ari_val=hca_res$adj.rand.index
        
        ##save(hca_res,file="hca_res.Rda")
        
        classgroup=hca_res$classgroup[,1]
      }
        #rownames(X)<-NULL
        
        setwd(outloc)
        
        feat.cor.mat=NA
        
        if(is.na(abs.cor.thresh)==FALSE){
          
          
        #  print("Starting pairwise correlation analysis")
          feat.cor.mat=plot.pairwise.correlation(data_matrix=X,newdevice=FALSE,abs.cor.thresh=abs.cor.thresh,pvalue.thresh=pvalue.thresh,cor.fdrthresh=cor.fdrthresh,cex.plots=cex.plots,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, plots.type=plots.type,Y=classgroup,cor.method=cor.method)
          
          
          if(is.na(Zmat)==FALSE){
            
            
            plot.pairwise.correlation.bipartite(data_matrix=X,newdevice=FALSE,abs.cor.thresh=abs.cor.thresh,pvalue.thresh=pvalue.thresh,cor.fdrthresh=cor.fdrthresh,cex.plots=cex.plots,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, plots.type=plots.type,data_matrixB=Zmat,cor.method=cor.method)
          }
        }
        
        dev.off()
        
        
        
        
        return(list("data_matrix_list"=data_matrix_list,"Silhouette.sample"=Silhouette.samples,
                    "Silhouette.feature"=Silhouette.features,"ari_val"=ari_val,
                   #"Gap.feature"=gap.stat.features,
                    #"Gap.sample"=gap.stat.samples,
                    "feat.cor.mat"=feat.cor.mat))
        
      }
      
    }else{
      
      print(paste("Invalid option: ",normalization.method[i],sep=""))
    }
    
    
  })
  
  #   ##save(data_matrix_list,file="data_matrix_list1.Rda")
  if(FALSE){   
  Gap.sample=lapply(1:length(data_matrix_list),function(j){
    
    return(data_matrix_list[[j]]$Gap.sample)
  })
  
  Gap.feature=lapply(1:length(data_matrix_list),function(j){
    
    return(data_matrix_list[[j]]$Gap.feature)
  })
}
  Silhouette.sample=lapply(1:length(data_matrix_list),function(j){
    
    return(data_matrix_list[[j]]$Silhouette.sample)
  })
  
  #Silhouette.feature=data_matrix_list$Silhouette.feature
  Silhouette.feature=lapply(1:length(data_matrix_list),function(j){
    
    return(data_matrix_list[[j]]$Silhouette.feature)
  })
  
  # ari_val=data_matrix_list$ari_val
  
  ari_val=lapply(1:length(data_matrix_list),function(j){
    
    return(data_matrix_list[[j]]$ari_val)
  })
  
  #feat.cor.mat=data_matrix_list$feat.cor.mat
  feat.cor.mat=lapply(1:length(data_matrix_list),function(j){
    
    return(data_matrix_list[[j]]$feat.cor.mat)
  })
  
  # data_matrix_list=data_matrix_list$data_matrix_list
  data_matrix_list=lapply(1:length(data_matrix_list),function(j){
    
    return(data_matrix_list[[j]]$data_matrix_list)
  })
  
  norm_method_names<-(normalization.method)
  norm_method_names<-gsub(norm_method_names,pattern="_norm|norm",replacement="")
  
  setwd(parentoutput_dir)
  
  #save(norm_method_names,Silhouette.feature,Silhouette.sample,ari_val,file="hcaplotdebug.Rda")
  
  pdf("Compare.HCA.performance.pdf",height=8,width=12)
  
  par(mfrow=c(1,3))
  par(mar=c(10,3,4,3))
  #    ##save(varimp_res2,data_m_fc,rf_classlabels,sorted_varimp_res,file="test_rf.Rda")
  #xaxt="n",

  cval1<-unlist(Silhouette.feature)
  names(cval1)<-norm_method_names
  x=barplot(cval1,ylab="",main="mean Silhouette coefficient \n for clustering of features",cex.axis=0.9,cex.names=0.9, xlab="",las=2,ylim=c(0,1))
 # title(ylab = "mean Silhouette coefficient \n for clustering of features", cex.lab = 1.5,line = 4.5)
            
  cval1<-unlist(Silhouette.sample)
  names(cval1)<-norm_method_names
  x=barplot(cval1,ylab="",main="mean Silhouette coefficient \n for clustering of samples",cex.axis=0.9,cex.names=0.9, xlab="",las=2,ylim=c(0,1))
 # title(ylab = "mean Silhouette coefficient \n for clustering of samples", cex.lab = 1.5,line = 4.5)
 if(FALSE){ 
  cval1<-unlist(Gap.feature)
  x=barplot(cval1,xlab="",main="",cex.axis=0.9,cex.names=0.9, ylab="",las=2,ylim=range(pretty(c(min(cval1,na.rm=TRUE),max(cval1,na.rm=TRUE)))))
  title(ylab = "mean Gap statistic \n for clustering of features", cex.lab = 1.5,line = 4.5)
  
  cval1<-unlist(Gap.sample)
  x=barplot(cval1,xlab="",main="",cex.axis=0.9,cex.names=0.9, ylab="",las=2,ylim=range(pretty(c(min(cval1,na.rm=TRUE),max(cval1,na.rm=TRUE)))))
  title(ylab = "mean Gap statistic \n for clustering of samples", cex.lab = 1.5,line = 4.5)
 }  
  
  cval1<-unlist(ari_val)
  print(cval1)
  names(cval1)<-norm_method_names
  x=barplot(cval1,ylab="",main="Adjusted rand index \n for comparing sample clustering \nwith the ground truth (class labels)",cex.axis=0.9,
            cex.names=0.9, xlab="",las=2,ylim=c(0,1))
 # title(ylab = "Adjusted rand index \n for comparing sample clustering \nwith the ground truth (class labels)", cex.lab = 1.5,line = 4.5)
  
  #par(mar=c(1,1))
  
  if(FALSE){
    plot(unlist(Silhouette.feature),type="o",ylim=c(0,1),col="#0072B2",lwd=2,xaxt="n",xlab="",ylab="Value",
       main="Comparison of clustering quality evaluation \n (silhouette, gap, and adjusted rand index) metrics",cex.main=0.8)
  lines(unlist(Silhouette.sample),type="o",lty=2,col="#E69F00",lwd=2)
  lines(unlist(Gap.feature),type="o",lty=5,col="#009E73",lwd=2)
  lines(unlist(Gap.sample),type="o",lty=3,col="#56B4E9",lwd=2)
  lines(unlist(ari_val),type="o",lty=4,col="#D55E00",lwd=2)
  axis(labels=norm_method_names,at=seq(1,length(norm_method_names)),side=1,las=2,cex.axis=0.65)
  legend("topright",legend=c("Silhouette.features","Silhouette.samples","Gap.feautres","Gap.samples","Adjusted Rand Index"),
         lty=c(1,2,5,3,4),col=c("#0072B2","#E69F00",
                                                                                                                      "#009E73","#56B4E9",
                                                                                                                      "#D55E00"))
  }
  
  dev.off()  
  #unlist(Gap.feature),unlist(Gap.sample),
  resmat<-cbind(normalization.method,unlist(Silhouette.feature),unlist(Silhouette.sample),unlist(ari_val))
  colnames(resmat)<-c("Normalization.method","Silhouette.feature","Silhouette.sample","Adjusted_Rand_Index")
  write.csv(resmat,file="comparison.HCA.performance.csv",row.names=FALSE)
  return(list(data_matrix_list=data_matrix_list,normalization.methods=normalization.method,
              Silhouette.sample=Silhouette.sample,Silhouette.feature=Silhouette.feature,
              #Gap.sample=Gap.sample,
              #Gap.feature=Gap.feature,
              ari_val=ari_val,feat.cor.mat=feat.cor.mat))
  
  
  
}

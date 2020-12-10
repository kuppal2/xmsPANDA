diffexp <-
function(Xmat=NA,Ymat=NA,feature_table_file,parentoutput_dir,class_labels_file,num_replicates=3,summarize.replicates=TRUE,summary.method="mean",summary.na.replacement="zeros",missing.val=0,rep.max.missing.thresh=0.3,
                  all.missing.thresh=0.1,group.missing.thresh=0.7,input.intensity.scale="raw",
                  log2transform=TRUE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=TRUE,lowess_norm=FALSE,madscaling=FALSE,TIC_norm=FALSE,rangescaling=FALSE,mstus=FALSE,paretoscaling=FALSE,sva_norm=FALSE,eigenms_norm=FALSE,vsn_norm=FALSE,
                  normalization.method=c("none"),rsd.filt.list=1,
                  pairedanalysis=FALSE,featselmethod=c("limma","pls"),fdrthresh=0.05,fdrmethod="BH",cor.method="spearman",networktype="complete",abs.cor.thresh=0.4,cor.fdrthresh=0.05,kfold=10,
                  pred.eval.method="BER",globalcor=TRUE,
                  target.metab.file=NA,target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=100, 
                  numtrees=20000,analysismode="classification",net_node_colors=c("green","red"), net_legend=TRUE,
                  svm_kernel="radial",heatmap.col.opt="brewer.RdBu",manhattanplot.col.opt=c("darkblue","red3"),boxplot.col.opt=c("journal"),barplot.col.opt=c("journal"),sample.col.opt="journal",lineplot.col.opt="journal",scatterplot.col.opt=c("journal"),hca_type="two-way",pls_vip_thresh=2,num_nodes=2,max_varsel=100,
                  pls_ncomp=5,pca.stage2.eval=FALSE,scoreplot_legend=TRUE,pca.global.eval=TRUE,rocfeatlist=seq(2,6,1),rocfeatincrement=TRUE,rocclassifier="svm",foldchangethresh=1,wgcnarsdthresh=20,WGCNAmodules=FALSE,
                  optselect=TRUE,max_comp_sel=1,saveRda=FALSE,legendlocation="topleft",pcacenter=TRUE,pcascale=TRUE,pca.cex.val=6,
                  pca.ellipse=FALSE,ellipse.conf.level=0.95,pls.permut.count=NA,svm.acc.tolerance=5,limmadecideTests=TRUE,pls.vip.selection="max",globalclustering=FALSE,plots.res=600,plots.width=8,plots.height=8,plots.type="cairo",
                  output.device.type="pdf",pvalue.thresh=0.05,individualsampleplot.col.opt="journal",pamr.threshold.select.max=FALSE,aggregation.method="RankAggreg",aggregation.max.iter=1000,mars.gcv.thresh=1,error.bar=TRUE,cex.plots=1,lme.modeltype="RI",
                  barplot.xaxis="Factor1",lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),match_class_dist=TRUE,timeseries.lineplots=FALSE,alphabetical.order=FALSE,
                  kegg_species_code="hsa",database="pathway",reference_set=NA,target.data.annot=NA,add.pvalues=FALSE,add.jitter=FALSE,fcs.permutation.type=1,
                  fcs.method="zscore",fcs.min.hits=2,names_with_mz_time=NA,
                  ylab_text="Abundance",xlab_text=NA,boxplot.type="ggplot",
                  samplermindex=NA,differential.network.analysis=TRUE,
                  degree.centrality.method="eigenvector",log2.transform.constant=1,
                  balance.classes=FALSE,balance.classes.sizefactor=10,
                  balance.classes.seed=1,cv.perm.count=100,multiple.figures.perpanel=FALSE,...)
{
  
  time_start<-Sys.time()
  options(warn=-1)
  
  
  print("**")
  print("**")
  
  #print(paste("g is ",group.missing.thresh,sep=""))
  modeltype=lme.modeltype
  balance.classes.method="ROSE"
  
  if(differential.network.analysis==TRUE){
    
    
    
    degree_rank_method="diffrank"
    
  }else{
    degree_rank_method="none"
  }
  
  
  
  feat.filt.thresh=NA
  feat_weight=1
  samplermindex=NA
  #pcacenter=TRUE
  #pcascale=TRUE
  alphacol=0.3
  
  
  print(parentoutput_dir)
  dir.create(parentoutput_dir)
  
  
  # print(is(group.missing.thresh<0.8))
  
  if(is.na(group.missing.thresh)==FALSE){
    if(group.missing.thresh<0.8){
      
      
      print("********************************************************************************")
      print("***** WARNING: group.missing.thresh is set to below 0.8. This can lead to false significance for class or group-wise comparisons.****")
      print("********************************************************************************")
      print("**")
      print("**")
      
    }
  }
  options(warn=-1)
  
  if(input.intensity.scale=="raw" || input.intensity.scale=="log2"){
    
    print("##################################################################################")
    print("Note 1: The order of samples should be same in the feature table and classlabels file")
    print(paste("Note 2: Treating input intensities as ",input.intensity.scale," values",sep=""))
    
  }else{
    
    stop("Input intensities should either be at raw or log2 scale")
  }
  
  
  suppressWarnings(suppressWarnings(sink(file=NULL)))
  x<-date()
  x<-strsplit(x,split=" ")
  
  #x<-gsub(x,pattern=":",replacement="_")
  targeted_feat_raw<-{}
  logistic_reg=FALSE
  x1<-unlist(x)
  x1<-gsub(x1,pattern=":",replacement="_")
  #fname<-paste(x1[2:5],collapse="_")
  
  #fname<-paste(x1[2:3],x1[5],x1[4],collapse="_")
  
  fname<-paste(x1[2:3],collapse="")
  
  #fname<-paste(fname,x1[6],sep="")
  x1[4]<-gsub(x1[4],pattern=":",replacement="_")
  fname<-paste(fname,x1[5],sep="")
  fname<-paste(fname,x1[4],sep="_")
  
  
  dir.create(parentoutput_dir)
  setwd(parentoutput_dir)
  
  
  #fname<-paste(parentoutput_dir,"/Log",fname,".txt",sep="")
  
  fname<-paste(parentoutput_dir,"/Log.txt",sep="")
  
  
  if(is.na(foldchangethresh)==FALSE){
    if(log2transform==TRUE && znormtransform==TRUE){
      
      stop("Both log2transform and znormtransform can not be true if foldchangethresh is not equal to NA.")
    }
  }
  
  if(featselmethod=="limma2way")
  {
    
    print("Note 3: lm2wayanova option is recommended for greater than 2x2 designs and this includes post-hoc comparisons")
    
    
  }
  
  
  if(featselmethod=="lm1wayanovarepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="limma1wayrepeat")
  {
    print("Note 3: Class labels format should be: Sample ID, Subject, Time. lm1wayanovarepeat is based on the nlme::lme() function with post-hoc Tukey HSD test.")
  }else{
    
    if(featselmethod=="lm2wayanovarepeat" | featselmethod=="spls2wayrepeat" | featselmethod=="limma2wayrepeat")
    {
      print("Note 3: Class labels format should be: Sample ID, Subject, Factor, Time. lm2wayanovarepeat is based on the nlme::lme() funciton with post-hoc Tukey HSD test. ")
    }
    
  }
  
  
  print("##############################Starting processing now################################")
  print(paste("**Program is running. Please check the logfile for runtime status: ",fname,"**",sep=""))
  
  fname_params<-paste(parentoutput_dir,"/InputParameters.csv",sep="")
  #sink(fname_params)
  # ###savelist=ls(),file="cur.Rda")
  c1<-"feature_table_file:"
  c2<-feature_table_file
  #c2<-rbind(c2,feature_table_file)
  c1<-rbind(c1,"parentoutput_dir:")
  c2<-rbind(c2,parentoutput_dir)
  c1<-rbind(c1,"class_labels_file:")
  c2<-rbind(c2,class_labels_file)
  
  c1<-rbind(c1,"num_replicates:")
  c2<-rbind(c2,num_replicates)
  c1<-rbind(c1,"summarize.replicates:")
  c2<-rbind(c2,summarize.replicates)
  c1<-rbind(c1,"summary.method:")
  c2<-rbind(c2,summary.method)
  c1<-rbind(c1,"summary.na.replacement:")
  c2<-rbind(c2,summary.na.replacement)
  c1<-rbind(c1,"rep.max.missing.thresh:")
  c2<-rbind(c2,rep.max.missing.thresh)
  c1<-rbind(c1,"all.missing.thresh:")
  c2<-rbind(c2,all.missing.thresh)
  c1<-rbind(c1,"group.missing.thresh:")
  c2<-rbind(c2,group.missing.thresh)
  c1<-rbind(c1,"input.intensity.scale:")
  c2<-rbind(c2,input.intensity.scale)
  
  c1<-rbind(c1,"normalization.method:")
  c2<-rbind(c2,normalization.method)
  
  
  c1<-rbind(c1,"log2transform:")
  c2<-rbind(c2,log2transform)
  c1<-rbind(c1,"medcenter:")
  c2<-rbind(c2,medcenter)
  c1<-rbind(c1,"znormtransform:")
  c2<-rbind(c2,znormtransform)
  c1<-rbind(c1,"quantile_norm:")
  c2<-rbind(c2,quantile_norm)
  
  
  c1<-rbind(c1,"TIC_norm:")
  c2<-rbind(c2,TIC_norm)
  c1<-rbind(c1,"lowess_norm:")
  c2<-rbind(c2,lowess_norm)
  c1<-rbind(c1,"madscaling:")
  c2<-rbind(c2,madscaling)
  
  c1<-rbind(c1,"rsd.filt.list:")
  c2<-rbind(c2,rsd.filt.list)
  c1<-rbind(c1,"pairedanalysis:")
  c2<-rbind(c2,pairedanalysis)
  c1<-rbind(c1,"featselmethod:")
  
  c2<-rbind(c2,paste(featselmethod,collapse=";"))
  
  c1<-rbind(c1,"pvalue.thresh:")
  c2<-rbind(c2,pvalue.thresh)
  c1<-rbind(c1,"fdrthresh:")
  c2<-rbind(c2,fdrthresh)
  c1<-rbind(c1,"fdrmethod:")
  c2<-rbind(c2,fdrmethod)
  c1<-rbind(c1,"cor.method:")
  c2<-rbind(c2,cor.method)
  c1<-rbind(c1,"abs.cor.thresh:")
  c2<-rbind(c2,abs.cor.thresh)
  c1<-rbind(c1,"cor.fdrthresh:")
  c2<-rbind(c2,cor.fdrthresh)
  c1<-rbind(c1,"kfold:")
  c2<-rbind(c2,kfold)
  c1<-rbind(c1,"globalcor:")
  c2<-rbind(c2,globalcor)
  c1<-rbind(c1,"target.metab.file:")
  c2<-rbind(c2,target.metab.file)
  c1<-rbind(c1,"target.mzmatch.diff:")
  c2<-rbind(c2,target.mzmatch.diff)
  c1<-rbind(c1,"target.rtmatch.diff:")
  c2<-rbind(c2,target.rtmatch.diff)
  c1<-rbind(c1,"max.cor.num:")
  c2<-rbind(c2,max.cor.num)
  c1<-rbind(c1,"missing.val:")
  c2<-rbind(c2,missing.val)
  c1<-rbind(c1,"networktype:")
  c2<-rbind(c2,networktype)
  c1<-rbind(c1,"samplermindex:")
  c2<-rbind(c2,samplermindex)
  c1<-rbind(c1,"numtrees:")
  c2<-rbind(c2,numtrees)
  c1<-rbind(c1,"analysismode:")
  c2<-rbind(c2,analysismode)
  c1<-rbind(c1,"net_node_colors:")
  c2<-rbind(c2,net_node_colors)
  c1<-rbind(c1,"net_legend:")
  c2<-rbind(c2,net_legend)
  c1<-rbind(c1,"heatmap.col.opt:")
  c2<-rbind(c2,heatmap.col.opt)
  c1<-rbind(c1,"manhattanplot.col.opt:")
  c2<-rbind(c2,paste(manhattanplot.col.opt,collapse=";"))
  c1<-rbind(c1,"boxplot.col.opt:")
  c2<-rbind(c2,boxplot.col.opt)
  
  c1<-rbind(c1,"barplot.col.opt:")
  c2<-rbind(c2,barplot.col.opt)
  
  c1<-rbind(c1,"sample.col.opt:")
  c2<-rbind(c2,sample.col.opt)
  c1<-rbind(c1,"alphacol:")
  c2<-rbind(c2,alphacol)
  c1<-rbind(c1,"pls_vip_thresh:")
  c2<-rbind(c2,pls_vip_thresh)
  c1<-rbind(c1,"alphacol:")
  c2<-rbind(c2,alphacol)
  c1<-rbind(c1,"max_varsel:")
  c2<-rbind(c2,max_varsel)
  c1<-rbind(c1,"pls_ncomp:")
  c2<-rbind(c2,pls_ncomp)
  c1<-rbind(c1,"pcacenter:")
  c2<-rbind(c2,pcacenter)
  c1<-rbind(c1,"pcascale:")
  c2<-rbind(c2,pcascale)
  c1<-rbind(c1,"pred.eval.method:")
  c2<-rbind(c2,pred.eval.method)
  c1<-rbind(c1,"rocfeatlist:")
  c2<-rbind(c2,paste(rocfeatlist,collapse=";"))
  c1<-rbind(c1,"rocfeatincrement:")
  c2<-rbind(c2,rocfeatincrement)
  c1<-rbind(c1,"rocclassifier:")
  c2<-rbind(c2,rocclassifier)
  c1<-rbind(c1,"foldchangethresh:")
  c2<-rbind(c2,foldchangethresh)
  c1<-rbind(c1,"wgcnarsdthresh:")
  c2<-rbind(c2,wgcnarsdthresh)
  c1<-rbind(c1,"WGCNAmodules:")
  c2<-rbind(c2,WGCNAmodules)
  c1<-rbind(c1,"optselect:")
  c2<-rbind(c2,optselect)
  c1<-rbind(c1,"max_comp_sel:")
  c2<-rbind(c2,max_comp_sel)
  c1<-rbind(c1,"saveRda:")
  c2<-rbind(c2,saveRda)
  c1<-rbind(c1,"pca.cex.val:")
  c2<-rbind(c2,pca.cex.val)
  c1<-rbind(c1,"pls.permut.count:")
  c2<-rbind(c2,pls.permut.count)
  c1<-rbind(c1,"pca.ellipse:")
  c2<-rbind(c2,pca.ellipse)
  c1<-rbind(c1,"ellipse.conf.level:")
  c2<-rbind(c2,ellipse.conf.level)
  c1<-rbind(c1,"legendlocation:")
  c2<-rbind(c2,legendlocation)
  c1<-rbind(c1,"svm.acc.tolerance:")
  c2<-rbind(c2,svm.acc.tolerance)
  c1<-rbind(c1,"limmadecideTests:")
  c2<-rbind(c2,limmadecideTests)
  c1<-rbind(c1,"pls.vip.selection:")
  c2<-rbind(c2,pls.vip.selection)
  c1<-rbind(c1,"globalclustering:")
  c2<-rbind(c2,globalclustering)
  c1<-rbind(c1,"plots.res:")
  c2<-rbind(c2,plots.res)
  c1<-rbind(c1,"plots.width:")
  c2<-rbind(c2,plots.width)
  c1<-rbind(c1,"plots.height:")
  c2<-rbind(c2,plots.height)
  c1<-rbind(c1,"plots.type:")
  c2<-rbind(c2,plots.type)
  c1<-rbind(c1,"output.device.type:")
  c2<-rbind(c2,output.device.type)
  c1<-rbind(c1,"pamr.threshold.select.max:")
  c2<-rbind(c2,pamr.threshold.select.max)
  c1<-rbind(c1,"aggregation.method:")
  c2<-rbind(c2,aggregation.method)
  c1<-rbind(c1,"mars.gcv.thresh")
  c2<-rbind(c2,mars.gcv.thresh)
  c1<-rbind(c1,"error.bar")
  c2<-rbind(c2,error.bar)
  c1<-rbind(c1,"timeseries.lineplots")
  c2<-rbind(c2,timeseries.lineplots)
  
  
  c1<-cbind(c1,c2)
  c1<-as.data.frame(c1)
  
  colnames(c1)<-c("InputParameter:","Value")
  write.csv(c1,file=fname_params,row.names=FALSE)
  rm(c1)
  
  
  
  
  sink(fname)
  print(sessionInfo())
  analysistype="oneway"
  
  if(featselmethod=="limma2way" | featselmethod=="lm2wayanova" | featselmethod=="spls2way"){
    analysistype="twowayanova"
  }else{
    
    if(featselmethod=="limma2wayrepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="spls2wayrepeat"){
      analysistype="twowayrepeat"
      pairedanalysis=TRUE
    }else{
      
      if(featselmethod=="limma1wayrepeat" | featselmethod=="lm1wayanovarepeat" | featselmethod=="spls1wayrepeat"){
        analysistype="onewayrepeat"
        pairedanalysis=TRUE
      }
      
    }
    
  }
  
  if(length(featselmethod)>1){
    
    if(length(rsd.filt.list)>1){
      
      print("Warning: only one RSD threshold allowed for multiple feature selection methods. Only the first RSD threshold will be used.")
      rsd.filt.list=rsd.filt.list[1]
      
    }
    
    consensus_res<-{}
    diffexp.res<-new("list")
    consensus_analysis=TRUE
    ranked_list<-{}
    common_feats<-{}
    
    pass_method_list<-{}
    
    for(i in 1:length(featselmethod))
    {
      
      if(featselmethod[i]=="rfesvm" && analysismode=="regression"){
        try(dev.off(),silent=TRUE)
        next;
      }
      
      outloc<-paste(parentoutput_dir,featselmethod[i],sep="/")
      suppressWarnings(diffexp.res[[i]]<-try(diffexp.child(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,summary.na.replacement,missing.val,rep.max.missing.thresh,
                                                           all.missing.thresh,group.missing.thresh,input.intensity.scale,
                                                           log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,TIC_norm,rangescaling,mstus,paretoscaling,sva_norm,eigenms_norm,vsn_norm,
                                                           normalization.method,rsd.filt.list,
                                                           pairedanalysis,featselmethod[i],fdrthresh,fdrmethod,cor.method,networktype,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
                                                           target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num,samplermindex,pcacenter,pcascale,
                                                           numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,barplot.col.opt,sample.col.opt,lineplot.col.opt, scatterplot.col.opt,hca_type,alphacol,pls_vip_thresh,num_nodes,max_varsel, pls_ncomp=pls_ncomp,pca.stage2.eval=pca.stage2.eval,scoreplot_legend=scoreplot_legend,pca.global.eval=pca.global.eval,rocfeatlist=rocfeatlist,rocfeatincrement=rocfeatincrement,rocclassifier=rocclassifier,foldchangethresh=foldchangethresh,wgcnarsdthresh=wgcnarsdthresh,WGCNAmodules=WGCNAmodules,
                                                           optselect=optselect,max_comp_sel=max_comp_sel,saveRda=saveRda,legendlocation=legendlocation,degree_rank_method=degree_rank_method,pca.cex.val=pca.cex.val,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,pls.permut.count=pls.permut.count,
                                                           svm.acc.tolerance=svm.acc.tolerance,limmadecideTests=limmadecideTests,pls.vip.selection=pls.vip.selection,globalclustering=globalclustering,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type,
                                                           output.device.type=output.device.type,pvalue.thresh,individualsampleplot.col.opt,pamr.threshold.select.max,mars.gcv.thresh,error.bar,cex.plots,modeltype,barplot.xaxis,lineplot.lty.option,match_class_dist=match_class_dist,
                                                           timeseries.lineplots=timeseries.lineplots,alphabetical.order=alphabetical.order,kegg_species_code=kegg_species_code,database=database,reference_set=reference_set,target.data.annot=target.data.annot,add.pvalues=add.pvalues,
                                                           add.jitter=add.jitter,fcs.permutation.type=fcs.permutation.type,fcs.method=fcs.method,fcs.min.hits=fcs.min.hits,names_with_mz_time=names_with_mz_time,ylab_text=ylab_text,xlab_text=xlab_text,boxplot.type=boxplot.type,
                                                           degree.centrality.method=degree.centrality.method,log2.transform.constant=log2.transform.constant,
                                                           balance.classes=balance.classes,balance.classes.sizefactor=balance.classes.sizefactor,
                                                           balance.classes.method=balance.classes.method,balance.classes.seed=balance.classes.seed,
                                                           cv.perm.count=cv.perm.count,multiple.figures.perpanel=multiple.figures.perpanel)
                                             
                                             ,silent=TRUE))
      if(is(diffexp.res[[i]],"try-error")){
        print(paste("Error processing option ",featselmethod[i],sep=""))
        print(paste("Error message: ",diffexp.res[[i]],sep=""))
        
        #diffexp.res[[i]]<-
        
        
        
      }else{
        pass_method_list<-c(pass_method_list,i)
        print(paste("Done with ",featselmethod[i],sep=""))
      }
    }
    
    ##save(diffexp.res,featselmethod,file="diffexp.res.Rda")
    
    max_num_feats<-max(c(unlist(lapply(1:length(diffexp.res),function(j){
      
      if(j%in%pass_method_list){
        return(dim(diffexp.res[[j]]$all_metabs)[1])
      }
      
    }
    ))),0)
    
    if(max_num_feats<1){
      
      return("No features selected.")
    }
    
    sel_feat_matrix<-matrix(0,nrow=max_num_feats,ncol=length(featselmethod))
    
    feat_rank_matrix<-matrix(1,nrow=max_num_feats,ncol=length(featselmethod))
    
    ranked_list<-{}
    
    for(i in 1:length(featselmethod))
    {
      if(is.na(aggregation.method)==TRUE){
        consensus_analysis=FALSE
      }else{
        
        if(aggregation.method=="none"){
          consensus_analysis=FALSE
        }else{
          if(aggregation.method=="consensus"){
            
            consensus_analysis=TRUE
            
          }
          
        }
      }
      if(length(diffexp.res[[i]]$all_metabs)>0 && consensus_analysis==TRUE){
        
        tvec=diffexp.res[[i]]$all_metabs$diffexp_rank
        
        if(length(tvec)<1){
          
          try(dev.off(),silent=TRUE)
          next;
        }
        cur_data_res<-diffexp.res[[i]]$all_metabs
        
        diffexp.res[[i]]$all_metabs<-diffexp.res[[i]]$all_metabs[order(diffexp.res[[i]]$all_metabs$mz),]
        mz_rt_all<-paste(diffexp.res[[i]]$all_metabs$mz,"_",diffexp.res[[i]]$all_metabs$time,sep="")
        mz_rt_selected<-paste(diffexp.res[[i]]$diffexp_metabs$mz,"_",diffexp.res[[i]]$diffexp_metabs$time,sep="")
        
        
        sel_feat_matrix[which(mz_rt_all%in%mz_rt_selected),i]<-1
        
        rownames(sel_feat_matrix)<-mz_rt_all
        
        if(i==1){
          
          #sort(rankingCriteria, index.return = TRUE)$ix
          ranked_indices=sort(tvec,index.return=TRUE)$ix
          ranked_list<-mz_rt_all[ranked_indices]
          
        }
        if(i>1){
          
          if(length(ranked_list)>0){
            
            ranked_list<-rbind(ranked_list,mz_rt_all[sort(tvec,index.return=TRUE)$ix])
          }else{
            
            ranked_indices=sort(tvec,index.return=TRUE)$ix
            ranked_list<-mz_rt_all[ranked_indices]
          }
          
          
        }
      }
    }
    
    ###save(sel_feat_matrix,file="feature.selection.different.methods.Rda")
    ###save(ranked_list,file="ranked_list.Rda")
    
    if(length(ranked_list)<1){
      consensus_analysis=FALSE
      
    }
    
    if(consensus_analysis==TRUE){
      
      print("################")
      print(paste("Aggregating results from different methods using ",aggregation.method," aggregation method.",sep=""))
      
      
      
      print("################")
      
      if(aggregation.method=="RankAggreg" | aggregation.method=="RankAggregGA"){
        
        if(max_varsel>dim(ranked_list)[2]){
          max_varsel=round(dim(ranked_list)[2]*0.3)
          
          
        }
        
        if(aggregation.method=="RankAggreg"){
          r1<-RankAggreg(x=ranked_list,k=max_varsel,verbose=TRUE,distance="Spearman",method="CE",maxIter=aggregation.max.iter)
        }else{
          
          r1<-RankAggreg(x=ranked_list,k=max_varsel,verbose=TRUE,distance="Spearman",method="GA",maxIter=aggregation.max.iter)
        }
        
        common_row_index<-which(mz_rt_all%in%r1$top.list)
        
        
        common_feats<-cur_data_res[common_row_index,]
        
        cnamesd1<-colnames(common_feats)
        time_ind<-which(cnamesd1=="time")
        mz_ind<-which(cnamesd1=="mz")
        
        
      }else{
        
        if(length(featselmethod)>=1){
          
          check_sel_status<-apply(sel_feat_matrix,1,sum)
          
          common_row_index<-which(check_sel_status==length(featselmethod))
          if(length(common_row_index)>0){
            common_feats<-cur_data_res[common_row_index,]
          }else{
            
            print("No features selected by all methods.")
          }
        }
        
        
      }
      
      cnames_1<-try(colnames(common_feats),silent=TRUE)
      
      
      if(consensus_analysis==FALSE | is(cnames_1,"try-error") | length(common_feats)<1){
        
        print("Skipping aggregation.")
      }else{
        
        
        bad_colind<-grep(tolower(cnames_1),pattern="rank")
        
        bad_colind_2<-grep(tolower(cnames_1),pattern="max.fold.change.log2")
        
        if(length(bad_colind_2)>1){
          bad_colind_2<-bad_colind_2[-c(1)]
        }
        
        bad_colind<-c(bad_colind,bad_colind_2)
        
        
        if(nrow(common_feats)>0){
          
          if(length(bad_colind)>0){
            common_feats<-common_feats[,-c(bad_colind)]
          }
        }
        
        
        common_feats<-unique(common_feats)
        print("Dimension of aggregated feature table of selected features:")
        print(dim(common_feats))
        
        num_common_feats<-dim(common_feats)[1]
        
        if(num_common_feats<1){
          
          stop("No common features found.")
        }
        ####savecommon_feats,file="common_feats.Rda")
        
        
        cnamesd1<-colnames(common_feats)
        time_ind<-which(cnamesd1=="time")
        mz_ind<-which(cnamesd1=="mz")
        
        
        #Xmat<-common_feats[,-c(1:2)]
        mz<-common_feats[,mz_ind]
        time<-common_feats[,time_ind]
        rnames1<-paste(mz,time,sep="_")
        
        Xmat<-cbind(mz,time,common_feats[,-c(1:time_ind)])
        
        rownames(Xmat)<-rnames1
        
        if(max_varsel>nrow(Xmat)){
          
          max_varsel<-nrow(Xmat)
        }
        
        #Ymat<-cbind(colnames(demetabs_res$norm_data[,-c(1:2)]),diffexp.res[i]$classlabels)
        Ymat<-diffexp.res[[i]]$classlabels
        
        
        outloc<-paste(parentoutput_dir,"AggregatedResults/",sep="/")
        
        if(log2transform==TRUE){
          
          input.intensity.scale="log2"
        }else{
          input.intensity.scale="raw"
        }
        
        dir.create(outloc)
        setwd(outloc)
        
        dir.create("Tables")
        
        write.table(Xmat,file="Tables/Aggregated_selected_features.txt",sep="\t",row.names=FALSE)
        
        if(nrow(Xmat)>1){
          
          dir.create("Figures")
          
          subdata<-t(Xmat[,-c(1:2)])
          classlabels<-Ymat[,2]
          classlabels<-as.data.frame(classlabels)
          if(analysismode=="classification"){
            svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist),silent=TRUE)
            #
            #svm_model<-svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist)
            
            
            if(is(svm_model,"try-error")){
              kfold_acc<-NA
              
            }else{
              kfold_acc<-svm_model$avg_acc
            }
          }else{
            
            svm_model_reg<-try(svm(x=subdata,y=(classlabels[,1]),type="eps",cross=kfold),silent=TRUE)
            
            if(is(svm_model_reg,"try-error")){
              kfold_acc<-NA
              
            }else{
              
              kfold_acc<-svm_model_reg$tot.MSE
            }
          }
          numcores<-num_nodes #round(detectCores()*0.6)
          
          kfold_acc_rand<-{}
          #for(p1 in 1:100){
          cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
          clusterEvalQ(cl,library(e1071))
          clusterEvalQ(cl,library(pROC))
          clusterEvalQ(cl,library(ROCR))
          clusterEvalQ(cl,library(CMA))
          clusterExport(cl,"svm_cv",envir = .GlobalEnv)
          clusterExport(cl,"svm",envir = .GlobalEnv)
          if(analysismode=="classification"){
            kfold_acc_rand<-parLapply(cl,1:100,function(p1){
              
              sample_ind<-sample(1:dim(classlabels)[1],size=dim(classlabels)[1])
              classlabels_rand<-classlabels[sample_ind,]
              classlabels_rand<-as.data.frame(classlabels_rand)
              svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels_rand,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist),silent=TRUE)
              #
              
              if(is(svm_model,"try-error")){
                # kfold_acc_rand<-c(kfold_acc_rand,NA)
                return(NA)
              }else{
                #kfold_acc_rand<-c(kfold_acc_rand,svm_model$avg_acc)
                return(svm_model$avg_acc)
              }
            })
            
          }else{
            
            kfold_acc_rand<-parLapply(cl,1:100,function(p1){
              
              sample_ind<-sample(1:dim(classlabels)[1],size=dim(classlabels)[1])
              classlabels_rand<-classlabels[sample_ind,]
              classlabels_rand<-as.data.frame(classlabels_rand)
              
              svm_model_reg<-try(svm(x=subdata,y=(classlabels[,1]),type="eps",cross=kfold),silent=TRUE)
              
              if(is(svm_model_reg,"try-error")){
                # kfold_acc_rand<-c(kfold_acc_rand,NA)
                return(NA)
              }else{
                #kfold_acc_rand<-c(kfold_acc_rand,svm_model$avg_acc)
                return(svm_model_reg$tot.MSE)
              }
            })
          }
          stopCluster(cl)
          
          kfold_acc_rand<-unlist(kfold_acc_rand)
          
          kfold_acc_rand<-mean(kfold_acc_rand,na.rm=TRUE)
          
          num_common_feats<-dim(common_feats)[1]
          summary_res<-cbind(num_common_feats,kfold_acc,kfold_acc_rand)
          colnames(summary_res)<-c("Number of selected features after aggregation",paste(pred.eval.method,"-accuracy",sep=""),paste(pred.eval.method," permuted accuracy",sep=""))
          
          
          file_name<-paste("../Results_summary_aggregated.txt",sep="")
          write.table(summary_res,file=file_name,sep="\t",row.names=FALSE)
          
          if(output.device.type=="pdf"){
            pdf("Figures/Aggregatedresults.pdf")
          }
          
          if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/HCA_aggregated_selectedfeats.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
          }
          
          #   print("here")
          #print(head(Ymat))
          ####saveXmat,file="Xmat.Rda")
          ####saveYmat,file="Ymat.Rda")
          
          g1<-get_hca(feature_table_file=NA,parentoutput_dir=outloc,class_labels_file=NA,X=Xmat,Y=Ymat,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode=analysismode,
                      sample.col.opt=sample.col.opt,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, plots.type=plots.type, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="",cexRow=0.4,cexCol=0.4,alphabetical.order=alphabetical.order,study.design=analysistype)
          
          if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
          }
          
          
          if(analysismode=="classification")
          {
            best_subset<-{}
            best_acc<-0
            
            xvec<-{}
            yvec<-{}
            for(i in 2:max_varsel){
              
              subdata<-t(Xmat[1:i,-c(1:2)])
              svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95,match_class_dist=match_class_dist),silent=TRUE)
              #svm_model<-svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95
              
              if(is(svm_model,"try-error")){
                
                svm_model<-NA
                
              }else{
                
                xvec<-c(xvec,i)
                yvec<-c(yvec,svm_model$avg_acc)
                if(svm_model$avg_acc>best_acc){
                  
                  best_acc<-svm_model$avg_acc
                  best_subset<-seq(1,i)
                  
                  
                }
                
                if(svm_model$avg_acc<best_acc){
                  
                  diff_acc<-best_acc-svm_model$avg_acc
                  if(diff_acc>50){
                    
                    break;
                    
                  }
                  
                }
              }
              
              
            }
            
            if(pred.eval.method=="CV"){
              ylab_text=paste(pred.eval.method," accuracy (%)",sep="")
              
            }else{
              if(pred.eval.method=="BER"){
                ylab_text=paste("Balanced accuracy"," (%)",sep="")
              }else{
                
                ylab_text=paste("AUC"," (%)",sep="")
              }
            }
            
            
            if(length(yvec)>0){
              
              msg1<-paste("k-fold CV classification accuracy based on forward selection of\n aggregated features ordered by ",featselmethod[1],sep="")
              
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/kfold_forward_selection_aggregated_selectedfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              
              plot(x=xvec,y=yvec,main=msg1,xlab="Feature index",ylab=ylab_text,type="b",col="brown")
              
              
              if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
              }
              
              
              cv_mat<-cbind(xvec,yvec)
              colnames(cv_mat)<-c("Feature Index",ylab_text)
              
              write.table(cv_mat,file="Tables/aggregated_kfold_cv_mat.txt",sep="\t")
            }
            
            
            
            if(output.device.type!="pdf"){
              
              temp_filename_1<-"Figures/Boxplots_aggregated_selectedfeats.png"
              
              #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              pdf(temp_filename_1)
              
            }
            
            if(log2transform==TRUE || input.intensity.scale=="log2"){
              
              if(znormtransform==TRUE){
                ylab_text_2="scale normalized"
              }else{
                if(quantile_norm==TRUE){
                  
                  ylab_text_2="quantile normalized"
                }else{
                  
                  ylab_text_2=""
                }
              }
              ylab_text=paste("log2 intensity ",ylab_text_2,sep="")
            }else{
              if(znormtransform==TRUE){
                ylab_text_2="scale normalized"
              }else{
                if(quantile_norm==TRUE){
                  
                  ylab_text_2="quantile normalized"
                }else{
                  ylab_text_2=""
                }
              }
              ylab_text=paste("Raw intensity ",ylab_text_2,sep="")
            }
            
            
            
            # par(mfrow=c(2,2))
            par(mfrow=c(1,1),family="sans",cex=cex.plots)
            get_boxplots(X=Xmat,Y=Ymat,parentoutput_dir=outloc,sample.col.opt=sample.col.opt,
                         boxplot.col.opt=boxplot.col.opt, alphacol=0.3,newdevice=FALSE,
                         cex.plots=cex.plots,ylabel=ylab_text,alphabetical.order=alphabetical.order,
                         boxplot.type=boxplot.type,study.design=analysistype,multiple.figures.perpanel=multiple.figures.perpanel)
            
            
            
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
            }
            
            
            
            if(output.device.type!="pdf"){
              
              temp_filename_1<-"Figures/Barplots_aggregated_selectedfeats.pdf"
              
              #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              pdf(temp_filename_1)
            }
            # par(mfrow=c(2,2))
            
            
            
            par(mfrow=c(1,1),family="sans",cex=cex.plots)
            
            get_barplots(feature_table_file=NA,class_labels_file=NA,X=Xmat,Y=Ymat,parentoutput_dir=outloc,newdevice=FALSE,ylabel=ylab_text,cex.plots=cex.plots,barplot.col.opt=barplot.col.opt,error.bar=error.bar,barplot.xaxis=barplot.xaxis,study.design=analysistype)
            
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
            }
            
            
            if(output.device.type!="pdf"){
              
              temp_filename_1<-"Figures/Individual_sample_plots_aggregated_selectedfeats.png"
              
              #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              pdf(temp_filename_1)
            }
            
            par(mfrow=c(1,1),family="sans",cex=cex.plots)
            get_individualsampleplots(feature_table_file=NA,class_labels_file=NA,X=Xmat,Y=Ymat,parentoutput_dir=outloc,newdevice=FALSE,
                                      ylabel=ylab_text,cex.plots=cex.plots,sample.col.opt=individualsampleplot.col.opt)
            
            
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
            }
            
            if(pairedanalysis==TRUE || timeseries.lineplots==TRUE)
            {
              
              
              if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/Lineplots_aggregated_selectedfeats.png"
                pdf(temp_filename_1)
                #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
              }
              
              par(mfrow=c(1,1),family="sans",cex=cex.plots)
              classlabels_orig<-read.table("../Stage2/classlabels_orig.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE, quote="")
              classlabels_orig<-classlabels_orig[,-c(1)]
              
              #save(Xmat,classlabels_orig,lineplot.col.opt,col_vec,pairedanalysis,
              #     pca.cex.val,pca.ellipse,ellipse.conf.level,legendlocation,ylab_text,error.bar,
               #    cex.plots,lineplot.lty.option,timeseries.lineplots,analysistype,file="debuga_lineplots.Rda")
              
              get_lineplots(X=Xmat,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),
                                class_labels_file=NA,lineplot.col.opt=lineplot.col.opt, 
                                alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,
                                point.cex.val=pca.cex.val,legendlocation=legendlocation,
                                pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,
                                filename="selected",ylabel=ylab_text,error.bar=error.bar,cex.plots=cex.plots,
                                lineplot.lty.option=lineplot.lty.option,timeseries.lineplots=timeseries.lineplots,
                                study.design=analysistype,multiple.figures.perpanel = multiple.figures.perpanel) #,silent=TRUE)
            }
            
            if(output.device.type!="pdf"){
              
              try(dev.off(),silent=TRUE)
            }
          }
          
          
          if(output.device.type=="pdf"){
            try(dev.off(),silent=TRUE)
          }
        }
      }
      print("#######################")
      print("#")
      print("#")
      print("Program ended successfully.")
      print("#")
      print("#")
      print("#######################")
      
      suppressWarnings(sink(file=NULL))
      
      #print("###############################")
      #print("###############################")
      #print("###############################")
      time_end<-Sys.time()
      
      time_taken_panda<-round(time_end-time_start,2)
      
      
      #print("*********")
      print(paste("**Program ended successfully in ",time_taken_panda," ",units(time_taken_panda),". Please see the ReadMe.txt file for description of output files and folders.**", sep=""))
      
      #print("*********")
      #print(paste("All result files are in the specified output location: ",parentoutput_dir,sep=""))
      #print("*********")
      #print("There will be a sub-folder for each step: Stage 1: pre-processing, Stage 2: statistical analysis (e.g. limma, PLS), and Stages 3 and 4: network analysis (Global and/or Targeted).")
      #print("*********")
      #print("")
      #print("*********")
      # print("Enjoy!")
      print("##############################END################################")
      
      
      
      
      s1<-"Stage 1 results: Preprocessing (Normalization, transformation)"
      s2<-"Stage 2 results: Feature selection & evaluation results for each feature selection method. Description of files and folder: A) *selected.features.final.txt file includes final table of selected features; B) Figures subfolder: includes Manhattan plots, boxplots, barplots, and other graphics, and C) Tables sub-folder: includes data files with feature selection results for all features, PCA (and PLS) scores and loadings, HCA clusters, k-fold CV results, and files for generating boxplots and barplots."
      s3<-"Stage 3 results: Correlation based network analysis"
      s4<-"Stage 4 results: Correlation based targeted network analysis"
      s5<-"Consensus results: HCA, k-fold CV, boxplots, and barplots for aggregated selected features"
      sm<-rbind(s1,s2,s3,s4,s5)
      setwd(parentoutput_dir)
      write.table(sm,file="ReadMe.txt",sep="\t",row.names=FALSE)
      return(list("individual.featsel.res"=diffexp.res,"aggregated.res"=common_feats))
    }else{
      
      print("#######################")
      print("#")
      print("#")
      print("Program ended successfully.")
      print("#")
      print("#")
      print("#######################")
      
      suppressWarnings(sink(file=NULL))
      
      # print("###############################")
      #print("###############################")
      #print("###############################")
      time_end<-Sys.time()
      
      time_taken_panda<-round(time_end-time_start,2)
      
      print(paste("**Program ended successfully in ",time_taken_panda," ",units(time_taken_panda),". Please see the ReadMe.txt file for description of output files and folders.**", sep=""))
      
      
      # print(paste("*******Program ended successfully in ",round(time_taken_panda,2)," minutes*******", sep=""))
      
      #  print("*     *")
      #  print("Consensus analysis could not be performed as not all features were selected by all feature selection methods.")
      #print("*     *")
      #print(paste("All result files are in the specified output location: ",parentoutput_dir,sep=""))
      # print("*     *")
      #print("There will be a sub-folder for each step: Stage 1: pre-processing, Stage 2: statistical analysis (e.g. limma, PLS), and Stages 3 and 4: network analysis (Global and/or Targeted).")
      #  print("*     *")
      
      #print("Please see the ReadMe.txt file for more information.")
      # print("*     *")
      # print("Enjoy!")
      print("##############################END################################")
      
      
      
      s1<-"Stage 1 results: Preprocessing (Normalization, transformation)"
      s2<-"Stage 2 results: Feature selection & evaluation results for each feature selection method. Description of files and folder: A) *selected.features.final.txt file includes final table of selected features; B) Figures subfolder: includes Manhattan plots, boxplots, barplots, and other graphics, and C) Tables sub-folder: includes data files with feature selection results for all features, PCA (and PLS) scores and loadings, HCA clusters, k-fold CV results, and files for generating boxplots and barplots."
      s3<-"Stage 3 results: Correlation based network analysis"
      s4<-"Stage 4 results: Correlation based targeted network analysis"
      
      sm<-rbind(s1,s2,s3,s4)
      setwd(parentoutput_dir)
      write.table(sm,file="ReadMe.txt",sep="\t",row.names=FALSE)
      options(warn=0)
      return(list("individual.featsel.res"=diffexp.res))
      
    }
    
    
    
  }else{
    
    suppressWarnings(
      diffexp.res<-diffexp.child(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,summary.na.replacement,missing.val,rep.max.missing.thresh,
                                 all.missing.thresh,group.missing.thresh,input.intensity.scale,
                                 log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,TIC_norm,rangescaling,mstus,paretoscaling,sva_norm,eigenms_norm,vsn_norm,
                                 normalization.method,rsd.filt.list,
                                 pairedanalysis,featselmethod,fdrthresh,fdrmethod,cor.method,networktype,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
                                 target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num,samplermindex,pcacenter,pcascale,
                                 numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,barplot.col.opt,sample.col.opt,lineplot.col.opt,scatterplot.col.opt,hca_type,alphacol,pls_vip_thresh,num_nodes,max_varsel, pls_ncomp,pca.stage2.eval,scoreplot_legend,pca.global.eval,rocfeatlist,rocfeatincrement,rocclassifier,foldchangethresh,wgcnarsdthresh,WGCNAmodules,
                                 optselect,max_comp_sel,saveRda,legendlocation,degree_rank_method,pca.cex.val,pca.ellipse,ellipse.conf.level,pls.permut.count,svm.acc.tolerance,limmadecideTests,pls.vip.selection,globalclustering,plots.res,plots.width,plots.height,
                                 plots.type,output.device.type,pvalue.thresh,individualsampleplot.col.opt,pamr.threshold.select.max,mars.gcv.thresh,error.bar,cex.plots,modeltype,barplot.xaxis,lineplot.lty.option,match_class_dist=match_class_dist,
                                 timeseries.lineplots=timeseries.lineplots,alphabetical.order=alphabetical.order,kegg_species_code=kegg_species_code,database=database,reference_set=reference_set,target.data.annot=target.data.annot,add.pvalues=add.pvalues,
                                 add.jitter=add.jitter,fcs.permutation.type=fcs.permutation.type,
                                 fcs.method=fcs.method,fcs.min.hits=fcs.min.hits,
                                 names_with_mz_time=names_with_mz_time,ylab_text=ylab_text,
                                 xlab_text=xlab_text,boxplot.type=boxplot.type,
                                 degree.centrality.method=degree.centrality.method,
                                 log2.transform.constant=log2.transform.constant,
                                 balance.classes=balance.classes,
                                 balance.classes.sizefactor=balance.classes.sizefactor,
                                 balance.classes.method=balance.classes.method,
                                 balance.classes.seed=balance.classes.seed,cv.perm.count=cv.perm.count,multiple.figures.perpanel=multiple.figures.perpanel)
    )
    
    time_end<-Sys.time()
    
    time_taken_panda<-round(time_end-time_start,2)
    
    
    print("#######################")
    print("#")
    print("#")
    print("Program ended successfully.")
    print("#")
    print("#")
    print("#######################")
    suppressWarnings(sink(file=NULL))
    
    print(paste("**Program ended successfully in ",time_taken_panda," ",units(time_taken_panda),". Please see the ReadMe.txt file for description of output files and folders.**", sep=""))
    
    #     print("###############################")
    #print("###############################")
    #print("###############################")
    #print("*     *")
    
    # print(paste("***Program ended successfully in ",round(time_taken_panda,2)," minutes***", sep=""))
    #  print(paste("*******Program ended successfully in ",round(time_taken_panda,2)," minutes*******", sep=""))
    
    # print("*     *")
    #    print(paste("All result files are in the specified output location: ",parentoutput_dir,sep=""))
    # print("*     *")
    #   print("There will be a sub-folder for each step: Stage 1: pre-processing, Stage 2: statistical analysis (e.g. limma, PLS), and Stages 3 and 4: network analysis (Global and/or Targeted).")
    # print("*     *")
    
    #   print("Please see the ReadMe.txt file for more information.")
    #   print("*     *")
    # print("Enjoy!")
    print("##############################END################################")
    
    
    s1<-"Stage 1 results: Preprocessing (Normalization, transformation)"
    s2<-"Stage 2 results: Feature selection & evaluation results for each feature selection method. Description of files and folder: A) *selected.features.final.txt file includes final table of selected features; B) Figures subfolder: includes Manhattan plots, boxplots, barplots, and other graphics, and C) Tables sub-folder: includes data files with feature selection results for all features, PCA (and PLS) scores and loadings, HCA clusters, k-fold CV results, and files for generating boxplots and barplots."
    s3<-"Stage 3 results: Correlation based network analysis"
    s4<-"Stage 4 results: Correlation based targeted network analysis"
    sm<-rbind(s1,s2,s3,s4)
    setwd(parentoutput_dir)
    write.table(sm,file="ReadMe.txt",sep="\t",row.names=FALSE)
    options(warn=0)
    return(diffexp.res)
    
    
  }
}

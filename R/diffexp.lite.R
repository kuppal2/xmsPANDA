diffexp.lite <-
function(Xmat=NA,Ymat=NA,outloc=NA,
                  summary.na.replacement="halffeaturemin",missing.val=0,
                  all.missing.thresh=0.1,group.missing.thresh=0.8,
                  input.intensity.scale="raw",
                
                  normalization.method=c("log2transform","znormtransform","lowess_norm","log2quantilenorm","quantile_norm","rangescaling",
                                         "paretoscaling","mstus","eigenms_norm","vsn_norm","sva_norm","tic_norm","cubicspline_norm","mad_norm","none"),
                  rsd.filt.list=1,
                  pca.global.eval=TRUE,
                  
                  pairedanalysis=FALSE,featselmethod=c("limma"),
                  pvalue.thresh=0.05,fdrthresh=0.05,fdrmethod="BH",
                  foldchangethresh=0,
                  analysismode="classification",
                  pls_vip_thresh=2,optselect=TRUE,max_comp_sel=2,
                  
                  hca_type="two-way",
                   
                  timeseries.lineplots=FALSE,alphabetical.order=FALSE,
                  
                  ylab_text="Abundance",boxplot.type="ggplot",
                  
                  color.palette=c("journal","npg","nejm","jco","lancet","custom1","brewer.RdYlBu","brewer.RdBu","brewer.PuOr","brewer.PRGn","brewer.PiYG","brewer.BrBG",
                                  "brewer.Set2","brewer.Paired","brewer.Dark2","brewer.YlGnBu","brewer.YlGn","brewer.YlOrRd","brewer.YlOrBr","brewer.PuBuGn",
                                  "brewer.PuRd","brewer.PuBu",
                                  "brewer.OrRd","brewer.GnBu","brewer.BuPu","brewer.BuGn","brewer.blues","black","grey65","terrain","rainbow","heat","topo"),
                  generate.boxplots=TRUE,
                  hca.cex.legend=0.7, lme.modeltype="lme.RI",globalcor=FALSE,abs.cor.thresh=0.4,cor.fdrthresh=NA,net_legend=TRUE,evaluate.classification.accuracy=FALSE,
                  limma.contrasts.type=c("contr.sum","contr.treatment"),limmadecideTests=FALSE,cex.plots=0.9,parentoutput_dir=NA, hca.labRow.value = TRUE,hca.labCol.value = TRUE,...)
{
  
  options(warn=-1)
  time_start<-Sys.time()
  options(warn=-1)
   
  if(is.na(evaluate.classification.accuracy)==FALSE){
    
    if(evaluate.classification.accuracy==TRUE){
        rocclassifier="svm"
    }else{
      rocclassifier=NA
      
    }
  }else{
    
    rocclassifier=NA
  }
  
  if(is.na(outloc)==TRUE){
    
    if(is.na(parentoutput_dir)==FALSE){
      outloc=parentoutput_dir
    }
  }
  dir.create(outloc,showWarnings = FALSE)
  
  setwd(outloc)
  sink(file="tmp.txt")
  
  suppressMessages(require(tidyverse))
  suppressMessages(require(tidyr))
  suppressMessages(require(MASS))
  suppressMessages(require(RColorBrewer))
  suppressMessages(require(lattice))
  suppressMessages(require(tidyverse))
  if(globalcor==TRUE){
    suppressMessages(require(igraph))
    suppressMessages(require(fdrtool))
  }
  if(is.na(rocclassifier)==FALSE){
    suppressMessages(require(e1071))
    suppressMessages(require(ROCR))
    suppressMessages(require(pROC))
    
  }
  suppressMessages(require(parallel))
  num_nodes=detectCores()*0.5
  
  if(length(which(featselmethod%in%c("limma","limma1way","limma2way","limma2wayrepeat","limma1wayrepeat")))>0){
    
    suppressMessages(require(limma))
  }else{
    if(length(which(featselmethod%in%c("pls","spls","spls2way","spls2wayrepeat","spls1wayrepeat","o1pls","o2pls")))>0){
      
      suppressMessages(require(mixOmics))
    }else{
      
      if(length(which(featselmethod%in%c("lm1wayanova","lm2wayanova","lm2wayanovarepeat","lm1wayanovarepeat","lm1wayrepeat")))>0){
        
        suppressMessages(require(nlme))
        suppressMessages(require(multcomp))
        suppressMessages(require(lsmeans))
        
      }else{
        
        if(length(which(featselmethod%in%c("RF","rf")))>0){
          suppressMessages(require(Boruta))
          
        }else{
          if(length(which(featselmethod%in%c("svmrfe","rfesvm")))>0){
            suppressMessages(require(e1071))
          }else{
            
            if(length(which(featselmethod%in%c("pamr")))>0){
              suppressMessages(require(pamr))
            }else{
              
              if(length(which(featselmethod%in%c("MARS")))>0){
                suppressMessages(require(earth))
              }
            }
          }
          
        }
      }
    }
    
  }
  try(sink(file=NULL),silent=TRUE)
  
  try(unlink("tmp.txt"),silent=TRUE)
  #print(paste("g is ",group.missing.thresh,sep=""))
  
  #if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
  #   Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  #}
  pvalue.dist.plot=FALSE
  lme.modeltype=lme.modeltype[1]
  
  if(lme.modeltype=="RI"){
    
    lme.modeltype="lme.RI"
  }else{
    if(lme.modeltype=="RIRS"){
      
      lme.modeltype="lme.RIRS"
    }
    
  }
  kfold=10
  output.device.type="png"
  feature_table_file=NA
  parentoutput_dir=outloc
  class_labels_file=NA
  num_replicates=1
  summarize.replicates=TRUE
  summary.method="mean"
  rep.max.missing.thresh=0.3
  
  log2transform=FALSE
  medcenter=FALSE
  znormtransform=FALSE
  quantile_norm=FALSE
  lowess_norm=FALSE
  madscaling=FALSE
  TIC_norm=FALSE
  rangescaling=FALSE
  mstus=FALSE
  paretoscaling=FALSE
  sva_norm=FALSE
  eigenms_norm=FALSE
  vsn_norm=FALSE
  normalization.method=normalization.method[1]
  cor.method="spearman"
  networktype=c("complete","GGM")
  
 # cor.fdrthresh=fdrthresh
  #kfold=10
  pred.eval.method="BER"
  
  target.metab.file=NA
  target.mzmatch.diff=10
  target.rtmatch.diff=NA
  max.cor.num=100
  numtrees=20000
  
  net_node_colors=c("green","red")
  
  network.label.cex=0.6
  svm_kernel="radial"
  
  pca.stage2.eval=FALSE
  scoreplot_legend=TRUE
 
  rocfeatlist=seq(1,5,1)
  rocfeatincrement=TRUE
  #rocclassifier=NA
  
  wgcnarsdthresh=20
  WGCNAmodules=FALSE
  

  
  max_varsel=100
  pls_ncomp=5
 
  saveRda=FALSE
  legendlocation="topleft"
  pcacenter=TRUE
  pcascale=TRUE
  pca.cex.val=6
  
  pca.ellipse=FALSE
  ellipse.conf.level=0.95
  pls.permut.count=NA
  svm.acc.tolerance=5
  
  pls.vip.selection="max"
  globalclustering=FALSE
  
  plots.res=600
  plots.width=10
  plots.height=8
  plots.type="cairo"
  xlab_text=NA
  
  
  pamr.threshold.select.max=FALSE
  aggregation.method="RankAggreg"
  aggregation.max.iter=1000
  mars.gcv.thresh=1
  
  error.bar=TRUE
  
 
  barplot.xaxis="Factor1"
  lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  match_class_dist=TRUE
  
  kegg_species_code="hsa"
  database="pathway"
  reference_set=NA
  target.data.annot=NA
  add.pvalues=FALSE
  add.jitter=FALSE
  fcs.permutation.type=1
  fcs.method="zscore"
  fcs.min.hits=2
  names_with_mz_time=NA
  
  samplermindex=NA
  differential.network.analysis=FALSE
  degree.centrality.method="hybrid.DEC"
  log2.transform.constant=1
  balance.classes=FALSE
  balance.classes.sizefactor=10
  balance.classes.seed=1
  cv.perm.count=NA
  multiple.figures.perpanel=FALSE
  
  heatmap.col.opt="brewer.RdBu"
  manhattanplot.col.opt=c("red3","darkblue")
  
  modeltype=lme.modeltype[1]
  balance.classes.method="ROSE"
  
  alpha.col=1
  similarity.matrix="correlation"
  outlier.method=NA #c("pcout","sumtukey","pcatukey","pcachisq")
  removeRda=TRUE
  
  plot_DiNa_graph=FALSE
  
  
  hca.labRow.value = TRUE
  hca.labCol.value = TRUE
  plot.boxplots.raw=FALSE
  vcovHC.type="HC3"
  ggplot.type1=TRUE
  facet.nrow=1
  pairwise.correlation.analysis=globalcor
  


  match_col.opt=match(color.palette,c("journal","npg","nejm","jco","lancet","custom1","brewer.RdYlBu","brewer.RdBu","brewer.PuOr","brewer.PRGn","brewer.PiYG","brewer.BrBG",
                                        "brewer.Set2","brewer.Paired","brewer.Dark2","brewer.YlGnBu","brewer.YlGn","brewer.YlOrRd","brewer.YlOrBr","brewer.PuBuGn",
                                        "brewer.PuRd","brewer.PuBu",
                                        "brewer.OrRd","brewer.GnBu","brewer.BuPu","brewer.BuGn","brewer.blues","black","grey65","terrain","rainbow","heat","topo"))
  
  
  match_col.opt=length(which(is.na(match_col.opt)==TRUE))
  
  #all colors match
  if(length(match_col.opt)<1){
    
    color.palette=color.palette[1]
  }else{
    
    if(length(grep(color.palette,pattern="brewer."))>1){
      
      color.palette=color.palette[1]
    }
  }
  
  color.palette=get_hexcolors_for_palettes(color.palette=color.palette,alpha.col=alpha.col[1])
  
  
  limma.contrasts.type=limma.contrasts.type[1]
  outlier.method=outlier.method[1]
  lineplot.lty.option=lineplot.lty.option[1]
  networktype=networktype[1]
  differential.network.analysis.method=c("degree.centrality")
  differential.network.analysis.method=differential.network.analysis.method[1] #"st5.differential.correlation"
  
  boxplot.col.opt=color.palette
  barplot.col.opt=color.palette
  sample.col.opt=color.palette
  lineplot.col.opt=color.palette
  scatterplot.col.opt=color.palette
  individualsampleplot.col.opt=color.palette
  
  if(length(grep(heatmap.col.opt,pattern = "brewer."))>0){
    
    heatmap.col.opt<-gsub(heatmap.col.opt,pattern="brewer.",replacement="")
  }
  
  labRow.value=hca.labRow.value
  labCol.value=hca.labCol.value
  
  if(is.na(parentoutput_dir)==TRUE){
    
    parentoutput_dir=getwd()
  }
  
  if(is.na(Xmat)==FALSE){
    
    feature_table_file=NA
  }
  if(is.na(Ymat)==FALSE){
    
    class_labels_file=NA
  }
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
  
  
  # print(parentoutput_dir)
  dir.create(parentoutput_dir,showWarnings = FALSE)
  
  
  # print(is(group.missing.thresh<0.8))
  
  if(is.na(group.missing.thresh)==FALSE){
    if(group.missing.thresh<0.8){
      
      
      
      
    }
  }
  options(warn=-1)
  
  if(input.intensity.scale=="raw" || input.intensity.scale=="log2"){
    
   # print("##################################################################################")
    #print("Note 1: The order of samples should be the same in the feature table and classlabels files")
    #print(paste("Note 2: Treating input intensities as ",input.intensity.scale," values",sep=""))
    
  }else{
    
    stop("Input intensities should either be at raw or log2 scale")
  }
  
  suppressMessages(suppressWarnings(try(sink(file=NULL),silent=TRUE)))
  #suppressMessages(suppressWarnings(sink(file=NULL)))
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
  
  
  suppressWarnings(dir.create(parentoutput_dir,showWarnings = FALSE))
  setwd(parentoutput_dir)
  
  
  #fname<-paste(parentoutput_dir,"/Log",fname,".txt",sep="")
  
  fname<-paste(parentoutput_dir,"/Log.txt",sep="")
  
  
  if(is.na(foldchangethresh)==FALSE){
    if(log2transform==TRUE && znormtransform==TRUE){
      
      stop("Both log2transform and znormtransform can not be true if foldchangethresh is not equal to NA.")
    }
  }
  
  if(featselmethod[1]=="limma2way")
  {
    
    #print("Note 3: lm2wayanova option is recommended for greater than 2x2 designs and this includes post-hoc comparisons")
    
    
  }
  
  
  if(featselmethod[1]=="lm1wayanovarepeat" | featselmethod[1]=="spls1wayrepeat" | featselmethod[1]=="limma1wayrepeat")
  {
   # print("Note 3: Class labels format should be: Sample ID, Subject, Time. lm1wayanovarepeat is based on the nlme::lme() function with post-hoc Tukey HSD test.")
  }else{
    
    if(featselmethod[1]=="lm2wayanovarepeat" | featselmethod[1]=="spls2wayrepeat" | featselmethod[1]=="limma2wayrepeat")
    {
    #  print("Note 3: Class labels format should be: Sample ID, Subject, Factor, Time. lm2wayanovarepeat is based on the nlme::lme() funciton with post-hoc Tukey HSD test. ")
    }
    
  }
  
  
  #print("##############################Starting processing now################################")
  #print(paste("**Program is running. Please check the logfile for runtime status: ",fname,"**",sep=""))
  
  cat("","\n")
  cat("Starting processing...")
  cat("","\n")
  
  fname_params<-paste(parentoutput_dir,"/InputParameters.csv",sep="")
  #sink(fname_params)
  # ###savelist=ls(),file="cur.Rda")
  c1={}
  c2={}
  #c2<-rbind(c2,feature_table_file)
  c1<-rbind(c1,"parentoutput_dir:")
  c2<-rbind(c2,parentoutput_dir)

  c1<-rbind(c1,"summary.na.replacement:")
  c2<-rbind(c2,summary.na.replacement)

  c1<-rbind(c1,"all.missing.thresh:")
  c2<-rbind(c2,all.missing.thresh)
  c1<-rbind(c1,"group.missing.thresh:")
  c2<-rbind(c2,group.missing.thresh)
  c1<-rbind(c1,"input.intensity.scale:")
  c2<-rbind(c2,input.intensity.scale)
  
  c1<-rbind(c1,"normalization.method:")
  c2<-rbind(c2,normalization.method[1])
  

  
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

  c1<-rbind(c1,"globalcor:")
  c2<-rbind(c2,globalcor)

  c1<-rbind(c1,"missing.val:")
  c2<-rbind(c2,missing.val)
 
  c1<-rbind(c1,"analysismode:")
  c2<-rbind(c2,analysismode)
  c1<-rbind(c1,"net_node_colors:")
  c2<-rbind(c2,net_node_colors)
  c1<-rbind(c1,"net_legend:")
  c2<-rbind(c2,net_legend)

  c1<-rbind(c1,"pls_vip_thresh:")
  c2<-rbind(c2,pls_vip_thresh)

 
  c1<-rbind(c1,"pls_ncomp:")
  c2<-rbind(c2,pls_ncomp)
  
  c1<-rbind(c1,"foldchangethresh:")
  c2<-rbind(c2,foldchangethresh)
  
 # c1<-rbind(c1,"optselect:")
  #c2<-rbind(c2,optselect)
  
  c1<-rbind(c1,"max_comp_sel:")
  c2<-rbind(c2,max_comp_sel)

  c1<-rbind(c1,"limmadecideTests:")
  c2<-rbind(c2,limmadecideTests)
  c1<-rbind(c1,"pls.vip.selection:")
  c2<-rbind(c2,pls.vip.selection)

  c1<-rbind(c1,"timeseries.lineplots")
  c2<-rbind(c2,timeseries.lineplots)
  
  c1<-rbind(c1,"alphabetical.order")
  c2<-rbind(c2,alphabetical.order)
  
  
  
  c1<-rbind(c1,"ylab_text")
  c2<-rbind(c2,ylab_text)
  
  
  c1<-rbind(c1,"boxplot.type")
  c2<-rbind(c2,boxplot.type)
  
  
  c1<-rbind(c1,"color.palette:")
  c2<-rbind(c2,paste(color.palette,collapse=";"))
  
  
 
  c1<-rbind(c1,"limma.contrasts.type")
  c2<-rbind(c2,limma.contrasts.type)
  
  c1<-rbind(c1,"hca.cex.legend")
  c2<-rbind(c2,hca.cex.legend)
  
 
  c1<-cbind(c1,c2)
  c1<-as.data.frame(c1)
  
  colnames(c1)<-c("InputParameter:","Value")
  write.csv(c1,file=fname_params,row.names=FALSE)
  rm(c1)
  rm(c2)
  
  
  
 # sink(fname)
#  print(sessionInfo())
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
      
     # print("Warning: only one RSD threshold allowed for multiple feature selection methods. Only the first RSD threshold will be used.")
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
      
      #suppressWarnings(
      diffexp.res[[i]]<-diffexp.child(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,summary.na.replacement,missing.val,rep.max.missing.thresh,
                                      all.missing.thresh,group.missing.thresh,input.intensity.scale,
                                      log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,TIC_norm,rangescaling,mstus,paretoscaling,sva_norm,eigenms_norm,vsn_norm,
                                      normalization.method[1],rsd.filt.list,
                                      pairedanalysis,featselmethod[i],fdrthresh,fdrmethod,cor.method,networktype,network.label.cex,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
                                      target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num,samplermindex,pcacenter,pcascale,
                                      numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,barplot.col.opt,sample.col.opt,lineplot.col.opt, scatterplot.col.opt,hca_type,alphacol,pls_vip_thresh,num_nodes,max_varsel,
                                      pls_ncomp=pls_ncomp,pca.stage2.eval=pca.stage2.eval,scoreplot_legend=scoreplot_legend,pca.global.eval=pca.global.eval,rocfeatlist=rocfeatlist,rocfeatincrement=rocfeatincrement,rocclassifier=rocclassifier,foldchangethresh=foldchangethresh,wgcnarsdthresh=wgcnarsdthresh,WGCNAmodules=WGCNAmodules,
                                      optselect=optselect,max_comp_sel=max_comp_sel,saveRda=saveRda,legendlocation=legendlocation,degree_rank_method=degree_rank_method,pca.cex.val=pca.cex.val,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,pls.permut.count=pls.permut.count,
                                      svm.acc.tolerance=svm.acc.tolerance,limmadecideTests=limmadecideTests,pls.vip.selection=pls.vip.selection,globalclustering=globalclustering,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type,
                                      output.device.type=output.device.type,pvalue.thresh,individualsampleplot.col.opt,pamr.threshold.select.max,mars.gcv.thresh,error.bar,cex.plots,modeltype,barplot.xaxis,lineplot.lty.option,match_class_dist=match_class_dist,
                                      timeseries.lineplots=timeseries.lineplots,alphabetical.order=alphabetical.order,kegg_species_code=kegg_species_code,database=database,reference_set=reference_set,target.data.annot=target.data.annot,add.pvalues=add.pvalues,
                                      add.jitter=add.jitter,fcs.permutation.type=fcs.permutation.type,fcs.method=fcs.method,fcs.min.hits=fcs.min.hits,names_with_mz_time=names_with_mz_time,ylab_text=ylab_text,xlab_text=xlab_text,boxplot.type=boxplot.type,
                                      degree.centrality.method=degree.centrality.method,log2.transform.constant=log2.transform.constant,
                                      balance.classes=balance.classes,balance.classes.sizefactor=balance.classes.sizefactor,
                                      balance.classes.method=balance.classes.method,balance.classes.seed=balance.classes.seed,
                                      cv.perm.count=cv.perm.count,
                                      multiple.figures.perpanel=multiple.figures.perpanel,labRow.value = labRow.value, 
                                      labCol.value = labCol.value,alpha.col=alpha.col,
                                      similarity.matrix=similarity.matrix,outlier.method=outlier.method[1],removeRda=removeRda,color.palette=color.palette,plot_DiNa_graph=plot_DiNa_graph,
                                      limma.contrasts.type=limma.contrasts.type,hca.cex.legend=hca.cex.legend,
                                      differential.network.analysis.method=differential.network.analysis.method,plot.boxplots.raw=plot.boxplots.raw,vcovHC.type=vcovHC.type,
                                      ggplot.type1=ggplot.type1,facet.nrow=facet.nrow,
                                      pairwise.correlation.analysis=pairwise.correlation.analysis,generate.boxplots=generate.boxplots,pvalue.dist.plot=pvalue.dist.plot)
      #,silent=TRUE)
      
      
      
      
      if(is(diffexp.res[[i]],"try-error")){
        print(paste("Error processing option ",featselmethod[i],sep=""))
        print(paste("Error message: ",diffexp.res[[i]],sep=""))
        
        #diffexp.res[[i]]<-
        
        
        
      }else{
        pass_method_list<-c(pass_method_list,i)
        #print(paste("Done with ",featselmethod[i],sep=""))
      }
      
      fname_del<-paste(outloc,"/Rplots.pdf",sep="")
      try(unlink(fname_del),silent=TRUE)
      
    }
    
    ##save(diffexp.res,featselmethod,file="diffexp.res.Rda")
    alphacol=alpha.col
    
    max_num_feats<-max(c(unlist(lapply(1:length(diffexp.res),function(j){
      
      if(j%in%pass_method_list){
        return(dim(diffexp.res[[j]]$all_metabs)[1])
      }
      
    }
    ))),0)
    
    if(max_num_feats<1){
      
      return("No features selected.")
    }
    
  
    
    ###save(sel_feat_matrix,file="feature.selection.different.methods.Rda")
    ###save(ranked_list,file="ranked_list.Rda")

    
 {
      
      
      
      time_end<-Sys.time()
      
      time_taken_panda<-round(time_end-time_start,2)
      
      
      cat(paste("**Program ended successfully in ",time_taken_panda," ",units(time_taken_panda),".**", sep=""),sep="\n")
      cat("",sep="\n")
    #  cat("Description of output folders:",sep="\n")
      
      
      
      s1<-"Stage 1 results: Includes raw (original) and preprocessed (filtered,normalized, and imputed) data tables"
      s2<-"Stage 2 results: Feature selection & evaluation results for each method. The following files and folders are generated for each method:: A) *selected.features.final.txt file includes final table of selected features; B) Figures subfolder: includes Manhattan plots, boxplots, HCA heatmap, and other figures, and C) Tables sub-folder: includes data files with feature selection results for all features, PCA (and PLS) scores and loadings, HCA clusters, k-fold CV results, and other tables.."
      
      sm<-rbind(s1,s2)
      
      #cat(sm,sep="\n")
      #cat("",sep="\n")
      #cat("##############################END################################",sep="\n")
      
      setwd(parentoutput_dir)
      write.table(sm,file="ReadMe.txt",sep="\t",row.names=FALSE)
      options(warn=0)
      return(list("individual.featsel.res"=diffexp.res))
      
    }
    
    
    
  }else{
    
    # suppressWarnings(
    diffexp.res<-diffexp.child(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,summary.na.replacement,missing.val,rep.max.missing.thresh,
                               all.missing.thresh,group.missing.thresh,input.intensity.scale,
                               log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,TIC_norm,rangescaling,mstus,paretoscaling,sva_norm,eigenms_norm,vsn_norm,
                               normalization.method,rsd.filt.list,
                               pairedanalysis,featselmethod,fdrthresh,fdrmethod,cor.method,networktype,network.label.cex,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
                               target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num,samplermindex,pcacenter,pcascale,
                               numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,barplot.col.opt,sample.col.opt,lineplot.col.opt,scatterplot.col.opt,hca_type,alphacol,pls_vip_thresh,num_nodes,max_varsel,
                               pls_ncomp,pca.stage2.eval,scoreplot_legend,pca.global.eval,rocfeatlist,rocfeatincrement,rocclassifier=rocclassifier,foldchangethresh,wgcnarsdthresh,WGCNAmodules,
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
                               balance.classes.seed=balance.classes.seed,cv.perm.count=cv.perm.count,
                               multiple.figures.perpanel=multiple.figures.perpanel,
                               labRow.value = labRow.value, labCol.value = labCol.value,alpha.col=alpha.col,
                               similarity.matrix=similarity.matrix,outlier.method=outlier.method[1],removeRda=removeRda,color.palette=color.palette,
                               plot_DiNa_graph=plot_DiNa_graph,limma.contrasts.type=limma.contrasts.type,hca.cex.legend=hca.cex.legend,
                               differential.network.analysis.method=differential.network.analysis.method,plot.boxplots.raw=plot.boxplots.raw,
                               vcovHC.type=vcovHC.type,ggplot.type1=ggplot.type1,facet.nrow=facet.nrow,
                               pairwise.correlation.analysis=pairwise.correlation.analysis,generate.boxplots=generate.boxplots,pvalue.dist.plot=pvalue.dist.plot) #,silent=TRUE)
    #)
    
    time_end<-Sys.time()
    
    time_taken_panda<-round(time_end-time_start,2)
    
    
 
   # print(paste("**Program ended successfully in ",time_taken_panda," ",units(time_taken_panda),".**", sep=""))
    
   
    #print("##############################END################################")
    cat(paste("**Program ended successfully in ",time_taken_panda," ",units(time_taken_panda),".**", sep=""),sep="\n")
    cat("",sep="\n")
    
    
    s1<-"Stage 1 results: Includes raw (original) and preprocessed (filtered,normalized, and imputed) data tables"
    s2<-"Stage 2 results: Feature selection & evaluation results for each method. The following files and folders are generated for each method:: A) *selected.features.final.txt file includes final table of selected features; B) Figures subfolder: includes Manhattan plots, boxplots, HCA heatmap, and other figures, and C) Tables sub-folder: includes data files with feature selection results for all features, PCA (and PLS) scores and loadings, HCA clusters, k-fold CV results, and other tables.."
    
    sm<-rbind(s1,s2)
   
   # cat("Description of output folders:",sep="\n")
    
    
    
    #s1<-"Stage 1 results: Includes raw (original) and preprocessed (filtered,normalized, and imputed) data tables"
    #s2<-"Stage 2 results: Feature selection & evaluation results for each method. The following files and folders are generated for each method:: A) *selected.features.final.txt file includes final table of selected features; B) Figures subfolder: includes Manhattan plots, boxplots, HCA heatmap, and other figures, and C) Tables sub-folder: includes data files with feature selection results for all features, PCA (and PLS) scores and loadings, HCA clusters, k-fold CV results, and other tables.."
    
    #sm<-rbind(s1,s2)
    
  #  cat(sm,sep="\n")
   # cat("",sep="\n")
  #  cat("##############################END################################",sep="\n")
    
    
    setwd(parentoutput_dir)
    write.table(sm,file="ReadMe.txt",sep="\t",row.names=FALSE)
    options(warn=0)
    return(diffexp.res)
    
    
  }
}

metabnet <-
function(feature_table_file,target.metab.file,sig.metab.file,class_labels_file=NA,parentoutput_dir,num_replicates=1,cor.method="spearman",abs.cor.thresh=0.4,cor.fdrthresh=0.05,
                   target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=100,feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method="mean",all.missing.thresh=0.5, group.missing.thresh=0.7,
                   missing.val=0, networktype="complete", samplermindex=NA,
                   rep.max.missing.thresh=0.3,summary.na.replacement="zeros",net_node_colors=c("pink","skyblue"),net_legend=FALSE,netrandseed=555,normalization.method=c("none"),
                   input.intensity.scale="raw",log2.transform.constant=1,...){
  
  options(warn=-1)
  
  
  
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  data_matrix<-data_preprocess(Xmat=NA,Ymat=NA,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh=NA,
                               summarize.replicates=summarize.replicates,summary.method=summary.method,
                               all.missing.thresh=all.missing.thresh,group.missing.thresh=group.missing.thresh,
                               missing.val=missing.val,samplermindex=samplermindex, rep.max.missing.thresh=rep.max.missing.thresh,summary.na.replacement=summary.na.replacement,featselmethod=featselmethod,
                               normalization.method=normalization.method,input.intensity.scale=input.intensity.scale,log2.transform.constant=log2.transform.constant)
  
  data_matrix<-data_matrix$data_matrix_afternorm_scaling
  
  
  
  
  dir.create(parentoutput_dir,showWarnings = FALSE)
  setwd(parentoutput_dir)
  
  data_m<-data_matrix[,-c(1:2)]
  
  #print(data_matrix[1:3,1:10])
  
  
  
  if(is.na(sig.metab.file)==FALSE){
    goodfeats<-read.table(sig.metab.file,sep="\t",header=TRUE)
  }else{
    goodfeats<-data_matrix
  }
  
  
  if(is.na(target.metab.file)==FALSE){
    dataA<-read.table(target.metab.file,sep="\t",header=TRUE)
    
    dataA<-as.data.frame(dataA)
    outloc<-paste(parentoutput_dir,"/metabnet","/",sep="")
    dir.create(outloc,showWarnings = FALSE)
    print(paste("Searching for metabolites matching target list",sep=""))
    
    g1<-getVenn(dataA=dataA,name_a="TargetSet",name_b="ExperimentalSet",dataB=data_matrix[,c(1:2)],mz.thresh=target.mzmatch.diff,time.thresh=target.rtmatch.diff,
                xMSanalyzer.outloc=outloc,alignment.tool=NA)
    #names(g1)
    
    
    
    if(length(g1$common)>1){
      
      if(is.na(sig.metab.file)==FALSE){
        com_mzs<-find.Overlapping.mzs(dataA=data_matrix,dataB=goodfeats,mz.thresh=1,time.thresh=1,alignment.tool=NA)
        
        sigfeats.index<-com_mzs$index.A #which(data_matrix$mz%in%goodfeats$mz)
        print(paste(length(unique(sigfeats.index))," selected features",sep=""))
        
        
      }else{
        sigfeats.index<-NA #seq(1,dim(data_matrix)[1])
        #sigfeats.index<-seq(1,50)
        
      }
      
      print(paste(length(unique(g1$common$index.B))," metabolites matched the target list",sep=""))
      
      print(paste("Generating targeted network",sep=""))
      
      
      if(networktype=="complete"){
        targetedan_fdr<-do_cor(data_matrix,subindex=sigfeats.index,targetindex=g1$common$index.B,outloc,networkscope="targeted",cor.method,
                               abs.cor.thresh,cor.fdrthresh,max.cor.num,net_node_colors,net_legend,netrandseed)
        
        #  pdf("Cornetworkplot.pdf")
        #load("metabnet.Rda")
        #print(plot(net_result))
        #dev.off()
        
      }else{
        if(networktype=="GGM"){
          targetedan_fdr<-get_partial_cornet(data_matrix, sigfeats.index,targeted.index=g1$common$index.B,networkscope="targeted",cor.method,abs.cor.thresh,
                                             cor.fdrthresh,outloc=outloc,net_node_colors,net_legend,netrandseed)
          
          
          
        }else{
          print("Invalid option. Please use complete or GGM.")
        }
      }
    }else{
      print(paste("Targeted metabolites were not found.",sep=""))
    }
  }else{
    outloc<-paste(parentoutput_dir,"/metabnet","/",sep="")
    #outloc<-paste(outloc,"/MWASresults","/",sep="")
    dir.create(outloc,showWarnings = FALSE)
    setwd(outloc)
    
    if(is.na(sig.metab.file)==FALSE){
      com_mzs<-find.Overlapping.mzs(dataA=data_matrix,dataB=goodfeats,mz.thresh=1,time.thresh=1,alignment.tool=NA)
      
      sigfeats.index<-com_mzs$index.A #which(data_matrix$mz%in%goodfeats$mz)
    }else{
      sigfeats.index<-NA #seq(1,dim(data_matrix)[1])
    }
    
    if(networktype=="complete"){
      
      
      #targetedan_fdr<-do_cor(data_matrix,sigfeats.index,outloc,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,max.cor.num)
      targetedan_fdr<-do_cor(data_matrix,subindex=sigfeats.index,targetindex=NA,outloc,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,
                             max.cor.num,net_node_colors,net_legend,netrandseed)
      
      if(FALSE){
        pdf("Cornetworkplot.pdf")
        load("metabnet.Rda")
        print(plot(net_result))
        dev.off()
      }
      
    }else{
      if(networktype=="GGM"){
        targetedan_fdr<-get_partial_cornet(data_matrix, sigfeats.index,targeted.index=NA,networkscope="global",cor.method,abs.cor.thresh,
                                           cor.fdrthresh,outloc=outloc,net_node_colors,net_legend,netrandseed)
        
        if(FALSE){
          pdf("GGMnetworkplot.pdf")
          load("metabnet.Rda")
          print(plot(net_result))
          dev.off()
        }
        
      }else{
        print("Invalid option")
      }
    }
    
  }
  
  
  print("Processing complete.")
  return(targetedan_fdr)
}

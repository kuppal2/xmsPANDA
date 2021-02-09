plot.pairwise.correlation <-
function(data_matrix,newdevice=FALSE,abs.cor.thresh=0.4,pvalue.thresh=0.05,cor.fdrthresh=0.2,cex.plots=0.7,output.device.type="pdf",plots.width=8,plots.height=8,plots.res=600, plots.type="cairo",Y=NA,cor.method="spearman"){
  
  goodfeats_temp=data_matrix
  
  suppressMessages(library(WGCNA))
  
  if(newdevice==TRUE){
    if(output.device.type!="pdf"){
      
      try(dev.off(),silent=TRUE)
      temp_filename_1<-"Pairwise.correlation.plots.pdf"
      pdf(temp_filename_1)
    }else{
      
      temp_filename_1<-"Pairwise.correlation.plots.png"
      
      png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      
    }
  }
  
  par(mfrow=c(1,1),family="sans",cex=cex.plots,cex.main=0.7)
  
  cor1<-WGCNA::cor(t(goodfeats_temp[,-c(1:2)]),method=cor.method)
  
  save(cor1,abs.cor.thresh,pvalue.thresh,cor.fdrthresh,file="corres.Rda")
  
  if(is.na(abs.cor.thresh)==FALSE){
    cor1[(abs(cor1)<abs.cor.thresh)]<-0
  }
  newnet <- cor1
  
  if(is.na(pvalue.thresh)==FALSE){
    
    corpval1=apply(cor1,2,function(x,goodfeats_temp){
      corPvalueStudent(x,n=ncol(goodfeats_temp[,-c(1:2)]))
      
    },goodfeats_temp=goodfeats_temp)
    
    
    if(is.na(cor.fdrthresh)==FALSE){
      fdr_adjust_pvalue<-try(fdrtool(as.vector(cor1[upper.tri(cor1)]),statistic="correlation",verbose=FALSE,plot=FALSE),silent=TRUE)
      
      newnet[upper.tri(newnet)][fdr_adjust_pvalue$qval > cor.fdrthresh] <- 0
    }                
    
    
    newnet[upper.tri(newnet)][as.vector(corpval1[upper.tri(corpval1)]) > pvalue.thresh] <- 0
  }
  
  newnet[lower.tri(newnet)] <- t(newnet)[lower.tri(newnet)]
  newnet <- as.matrix(newnet)
  
  corqval1=newnet
  diag(corqval1)<-0
  upperTriangle<-upper.tri(cor1, diag=F)
  lowerTriangle<-lower.tri(cor1, diag=F)
  
  if(is.na(cor.fdrthresh)==FALSE){
    corqval1[upperTriangle]<-fdr_adjust_pvalue$qval
    corqval1[lowerTriangle]<-corqval1[upperTriangle]
  }
  
  cor1=newnet
  rm(newnet)
  
  #  rownames(cor1)<-paste(goodfeats_temp[,c(1)],goodfeats_temp[,c(2)],sep="_")
  #  colnames(cor1)<-rownames(cor1)
  
  cor_range<-round(range(cor1[upperTriangle],na.rm=TRUE),2)
  
  mainlab1<-paste("Pairwise correlations between features\n correlation range: ",cor_range[1]," to ",cor_range[2],sep="")
  
  
  
  h1<-heatmap.2(cor1,col=brewer.pal(11,"RdBu"),Rowv=TRUE,Colv=TRUE,dendrogram="none",scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none",main=mainlab1,cexRow = 0.5,cexCol = 0.5,cex.main=0.7)
  
  #if(FALSE)
    {
    #plot_sample_correlation_child(X=goodfeats_temp[,-c(1:2)],abs.cor.thresh=abs.cor.thresh,pvalue.thresh=pvalue.thresh,cor.fdrthresh=cor.fdrthresh,groupname="all",cor.method=cor.method)
    
    
    if(is.na(Y)==FALSE){
      
      class_levels=levels(as.factor(Y))
      
      for(i in 1:length(class_levels)){
        plot_sample_correlation_child(X=goodfeats_temp[,which(Y==class_levels[i])+2],abs.cor.thresh=abs.cor.thresh,pvalue.thresh=pvalue.thresh,cor.fdrthresh=cor.fdrthresh,groupname=class_levels[i],cor.method=cor.method)
        
      }
    }
    
  }
  
  return(cor1)
  
}

plot_sample_correlation_child <-
function(X,abs.cor.thresh=0.4,pvalue.thresh=0.05,cor.fdrthresh=0.2,groupname="",cor.method="spearman"){
  
  
  cor1<-WGCNA::cor(X,method=cor.method)
  
  
  
  
  
  
  if(is.na(abs.cor.thresh)==FALSE){
    cor1[(abs(cor1)<abs.cor.thresh)]<-0
  }
  
  newnet <- cor1
  
  if(is.na(cor.fdrthresh)==FALSE){
    fdr_adjust_pvalue<-try(suppressWarnings(fdrtool(as.vector(cor1[upper.tri(cor1)]),statistic="correlation",verbose=FALSE,plot=FALSE)),silent=TRUE)
    
    
    newnet[upper.tri(newnet)][fdr_adjust_pvalue$qval > cor.fdrthresh] <- 0
  }
  
  if(is.na(pvalue.thresh)==FALSE){
    
    corpval1=apply(cor1,2,function(x){
      corPvalueStudent(x,n=ncol(X))
      
    })
    
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
  
  
  cor_range<-round(range(cor1[upperTriangle],na.rm=TRUE),2)
  
  mainlab1<-paste("Pairwise correlations between ",groupname," samples\n correlation range: ",cor_range[1]," to ",cor_range[2],sep="")
  
  
  
  h1<-heatmap.2(cor1,col=brewer.pal(11,"RdBu"),Rowv=TRUE,Colv=TRUE,dendrogram="none",scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none",main=mainlab1,cexRow = 0.5,cexCol = 0.5,cex.main=0.7)
  
}

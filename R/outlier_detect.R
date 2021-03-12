outlier_detect <-
function(data_matrix,ncomp=2,pthresh=0.005,outlier.method="sumtukey",column.rm.index=NA,iqr.multiplier=1.5){
  
  
  cnames<-{}
  if(outlier.method=="sumtukey"){
    
    
    if(is.na(column.rm.index)==FALSE){
      s2=apply(data_matrix[,-c(column.rm.index)],2,function(x){sum(x,na.rm=TRUE)})
      cnames<-colnames(data_matrix[,-c(column.rm.index)])
    }else{
      s2=apply(data_matrix,2,function(x){sum(x,na.rm=TRUE)})
      cnames<-colnames(data_matrix)
    }
    
    iqr_val<-quantile(s2,0.75)-quantile(s2,0.25)
    upper_limit=quantile(s2,0.75)+iqr.multiplier*(iqr_val)
    lower_limit=quantile(s2,0.25)-iqr.multiplier*(iqr_val)
    cnames<-cnames[which(s2>upper_limit | s2<lower_limit)]
    
    
  }
  else{
  if(is.na(column.rm.index)==FALSE){
  p1<-mixOmics::pca(t(data_matrix[,-c(column.rm.index)]),center=TRUE,scale=TRUE,ncomp=ncomp)
  cnames<-colnames(data_matrix[,-c(column.rm.index)])
  }else{
    p1<-mixOmics::pca(t(data_matrix),center=TRUE,scale=TRUE,ncomp=ncomp)
    cnames<-colnames(data_matrix)
    
  }
  U3=p1$variates$X
  
  ##save(p1,U3,file="pcaoutlier.Rda")
  
  if(outlier.method=="pcachisq"){
    suppressMessages(library(robust))
    dist2=covRob(p1$variates$X,estim="pairwiseGK")$dist
    pval <- pchisq(dist2, df = ncomp, lower.tail = FALSE)
    is.out <- (pval < (pthresh)) # / length(dist2)))
    #qplot(U3[, 1], U3[, 2], color =is.out,size = I(2),main="Outlier (green dots) detection using PCA and ChiSq test (p<0.005) ") + coord_equal()
    col=rep("blue",length(is.out))
    col[is.out==TRUE]<-"brown"
    plotIndiv(p1,col=col,cex=2,title=paste("Outlier (brown dots) detection using PCA and ChiSq test (p<",pthresh,")",sep=""),size.title=8) #0.005) )
    cnames<-cnames[which(is.out==TRUE)]
  }else{
    if(outlier.method=="pcout"){
      suppressMessages(library(mvoutlier))
      res<-pcout(U3,makeplot=TRUE)
      cnames<-names(res$wfinal01[which(res$wfinal01==0)])
    }else{
      
      if(outlier.method=="pcatukey"){
        #s2=apply(data_matrix[,-c(1:2)],2,sum)
        s2=U3[,1]
        iqr_val<-quantile(s2,0.75)-quantile(s2,0.25)
        upper_limit=quantile(s2,0.75)+iqr.multiplier*(iqr_val)
        lower_limit=quantile(s2,0.25)-iqr.multiplier*(iqr_val)
        cnames1<-cnames[which(s2>upper_limit | s2<lower_limit)]
        
       # s2=U3[,2]
       # iqr_val<-quantile(s2,0.75)-quantile(s2,0.25)
      #  upper_limit=quantile(s2,0.75)+iqr.multiplier*(iqr_val)
      #  lower_limit=quantile(s2,0.25)-iqr.multiplier*(iqr_val)
       # cnames1<-c(cnames1,cnames[which(s2>upper_limit | s2<lower_limit)])
        
        cnames=cnames1
      }else{
        
       
        
      }
    }
  }
  
  }
  if(length(cnames)>0){
  print(paste("Outliers: ",cnames,sep=""))
  }else{
    
    print("No outliers detected.")
  }
  write.table(cnames,file="Potential.Outliers.txt",sep="\t",row.names=FALSE)
  return(cnames)
}

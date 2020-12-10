plot.pairwise.correlation.bipartite <-
function(data_matrix,newdevice=FALSE,abs.cor.thresh=0.4,pvalue.thresh=0.05,cor.fdrthresh=0.2,cex.plots=0.7,output.device.type="pdf",plots.width=8,plots.height=8,plots.res=600, plots.type="cairo",data_matrixB,cor.method="spearman"){
  
  
  
  cnames<-colnames(data_matrixB)
  cnames<-tolower(cnames)
  
  check_names<-grep(cnames,pattern="^name$")
  
  names_with_mz_time<-NA
  
  X<-data_matrixB
  if(length(check_names)>0){
    
    if(check_names==1){
      
      check_names1<-grep(cnames,pattern="^mz$")
      check_names2<-grep(cnames,pattern="^time$")
      
      
      if(length(check_names1)<1 & length(check_names2)<1){
        mz<-seq(1.00001,nrow(X)+1,1)
        time<-seq(1.01,nrow(X)+1,1.00)
        check_ind<-gregexpr(cnames,pattern="^name$")
        check_ind<-which(check_ind>0)
        X<-as.data.frame(X)
        
        
        Name<-as.character(X[,check_ind])
        
        X<-cbind(mz,time,X[,-check_ind])
        names_with_mz_time=cbind(Name,mz,time)
        
        names_with_mz_time<-as.data.frame(names_with_mz_time)
        X<-as.data.frame(X)
        
        #write.table(names_with_mz_time,file="Stage1/Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        
      }else{
        
        if(length(check_names1)>0 & length(check_names2)>0){
          
          check_ind<-gregexpr(cnames,pattern="^name$")
          check_ind<-which(check_ind>0)
          Name<-as.character(X[,check_ind])
          X<-X[,-check_ind]
          names_with_mz_time=cbind(Name,X$mz,X$time)
          colnames(names_with_mz_time)<-c("Name","mz","time")
          names_with_mz_time<-as.data.frame(names_with_mz_time)
          X<-as.data.frame(X)
          #write.table(names_with_mz_time,file="Stage1/Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        }
      }
      
    }
  }else{
    
    
    
    check_names1<-grep(cnames[1],pattern="^mz$")
    check_names2<-grep(cnames[2],pattern="^time$")
    if(length(check_names1)<1 || length(check_names2)<1){
      stop("Invalid feature table format. First two columns should be mz and time. Please check example files.")
    }
  }
  
  data_matrixB<-X
  
  
  goodfeats_temp=rbind(data_matrix,data_matrixB)
  
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
  
  #save(goodfeats_temp,data_matrix,file="pcormat.Rda")
  
  cor1<-WGCNA::cor(t(goodfeats_temp[,-c(1:2)]),method=cor.method)
  
  cor1<-cor1[1:nrow(data_matrix),(nrow(data_matrix)+1):(nrow(data_matrix)+nrow(data_matrixB))]
  
  corpval1=apply(cor1,2,function(x){corPvalueStudent(x,n=ncol(goodfeats_temp[,-c(1:2)]))})
  
  #  ##save(cor1,goodfeats_temp,corpval1,file="cor2.Rda")
  
  fdr_adjust_pvalue<-try(fdrtool(as.vector(cor1[upper.tri(cor1)]),statistic="correlation",verbose=FALSE,plot=FALSE),silent=TRUE)
  
  cor1[(abs(cor1)<abs.cor.thresh)]<-0
  newnet <- cor1
  newnet[upper.tri(newnet)][fdr_adjust_pvalue$qval > cor.fdrthresh] <- 0
  newnet[upper.tri(newnet)][as.vector(corpval1[upper.tri(corpval1)]) > pvalue.thresh] <- 0
  
  newnet[lower.tri(newnet)] <- t(newnet)[lower.tri(newnet)]
  newnet <- as.matrix(newnet)
  
  corqval1=newnet
  diag(corqval1)<-0
  upperTriangle<-upper.tri(cor1, diag=F)
  lowerTriangle<-lower.tri(cor1, diag=F)
  corqval1[upperTriangle]<-fdr_adjust_pvalue$qval
  corqval1[lowerTriangle]<-corqval1[upperTriangle]
  
  cor1=newnet
  rm(newnet)
  
  # rownames(cor1)<-paste(goodfeats_temp[,c(1)],goodfeats_temp[,c(2)],sep="_")
  #  colnames(cor1)<-rownames(cor1)
  
  cor_range<-round(range(cor1[upperTriangle],na.rm=TRUE),2)
  
  mainlab1<-paste("Pairwise correlation of variables between two datasets\n correlation range: ",cor_range[1]," to ",cor_range[2],sep="")
  
  # print(mainlab1)
  #xlab="Samples",ylab="features"
  
  #dendrogram="none",
  h1<-heatmap.2(cor1,col=brewer.pal(11,"RdBu"),Rowv=TRUE,Colv=TRUE,scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none",main=mainlab1,cexRow = 0.5,cexCol = 0.5,cex.main=0.7,xlab="Dataset B",ylab="Dataset A")
  
  
  
  
  return(cor1)
  
}

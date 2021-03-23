plot.pairwise.correlation.bipartite <-
function(data_matrix,data_matrixB,newdevice=TRUE,abs.cor.thresh=0.4,pvalue.thresh=0.05,
                                              cor.fdrthresh=0.2,cex.plots=0.7,output.device.type="pdf",plots.width=8,
                                              plots.height=8,plots.res=600, plots.type="cairo",cor.method="spearman",
                                              do.clustering=TRUE){
  
 
  cnames<-colnames(data_matrixB)
  cnames<-tolower(cnames)
  
  check_names<-grep(cnames,pattern="^name$")
  
  names_with_mz_time<-NA
  
  X<-data_matrixB
  if(length(check_names)>0){
    
    if(check_names==1){
      
      rownames(data_matrix)<-data_matrix$Name
      rownames(data_matrixB)<-data_matrixB$Name
      check_names1<-grep(cnames,pattern="^mz$")
      check_names2<-grep(cnames,pattern="^time$")
      
      
      if(length(check_names1)<1 & length(check_names2)<1){
        
        check_ind<-gregexpr(cnames,pattern="^name$")
        check_ind<-which(check_ind>0)
        col_ind_rm<-c(check_ind)
        
      }else{
        
        if(length(check_names1)>0 & length(check_names2)>0){
          
        
          
          check_ind1<-gregexpr(cnames,pattern="^mz$")
          check_ind1<-which(check_ind1>0)
          
          check_ind2<-gregexpr(cnames,pattern="^time$")
          check_ind2<-which(check_ind2>0)
          
          col_ind_rm<-c(check_ind,check_ind1,check_ind2)
          
          #write.table(names_with_mz_time,file="Stage1/Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        }
      }
      
    }
  }else{
    
    
    
    check_names1<-grep(cnames[1],pattern="^mz$")
    check_names2<-grep(cnames[2],pattern="^time$")
    if(length(check_names1)<1 || length(check_names2)<1){
      stop("Invalid feature table format. The format should be either Name in column A or mz and time in columns A and B. Please check example files.")
    }
    if(length(check_names1)>0 & length(check_names2)>0){
      
      
      
      check_ind1<-gregexpr(cnames,pattern="^mz$")
      check_ind1<-which(check_ind1>0)
      
      check_ind2<-gregexpr(cnames,pattern="^time$")
      check_ind2<-which(check_ind2>0)
      
      col_ind_rm<-c(check_ind1,check_ind2)
      
      #write.table(names_with_mz_time,file="Stage1/Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
    }
  }
  
 # data_matrixB<-X
  
  
  goodfeats_temp=rbind(data_matrix,data_matrixB)
  print(newdevice)
  print(output.device.type)
  
  if(newdevice==TRUE){
    if(output.device.type=="pdf"){
      
      try(dev.off(),silent=TRUE)
      temp_filename_1<-"Pairwise.correlation.plots.pdf"
      
      #pdf(temp_filename_1)
      pdf(temp_filename_1,width=plots.width,height=plots.height)
      
    }else{
      
      temp_filename_1<-"Pairwise.correlation.plots.png"
      
      png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      
    }
  }
  
  par(mfrow=c(1,1),family="sans",cex=cex.plots,cex.main=0.7)
  
  #save(goodfeats_temp,data_matrix,file="pcormat.Rda")
  
  cor1<-WGCNA::cor(t(goodfeats_temp[,-c(col_ind_rm)]),method=cor.method)
  
  cor1<-cor1[1:nrow(data_matrix),(nrow(data_matrix)+1):(nrow(data_matrix)+nrow(data_matrixB))]
  
  corpval1=apply(cor1,2,function(x){corPvalueStudent(x,n=ncol(goodfeats_temp[,-c(col_ind_rm)]))})
  
  #save(cor1,goodfeats_temp,corpval1,file="cor2.Rda")
 
  cor1[(abs(cor1)<abs.cor.thresh)]<-0
  newnet <- cor1
  

  #upperTriangle<-upper.tri(cor1, diag=F)
  #lowerTriangle<-lower.tri(cor1, diag=F)
  
  if(is.na(cor.fdrthresh)==FALSE){
  fdr_adjust_pvalue<-try(suppressWarnings(fdrtool(as.vector(cor1),statistic="correlation",verbose=FALSE,plot=FALSE)),silent=TRUE)
  #corqval1[upperTriangle]<-fdr_adjust_pvalue$qval
  #
  newnet[fdr_adjust_pvalue$qval > cor.fdrthresh] <- 0
  newnet[as.vector(corpval1) > pvalue.thresh] <- 0
  
 # newnet[lower.tri(newnet)] <- t(newnet)[lower.tri(newnet)]
  newnet <- as.matrix(newnet)
  
  cor1=newnet
  }
  
  correlations_melted<-na.omit(melt(cor1, value.name ="correlationCoef"))
 # correlations_melted<-na.omit(melt(newnet[upper.tri(newnet)], value.name ="correlationCoef")) #use melt to reshape the matrix into triplets, na.omit to get rid of the NA rows
  
  colnames(correlations_melted)<-c("from", "to", "weight")

  #corqval1[lowerTriangle]<-corqval1[upperTriangle]
  
  #cor1=newnet
  rm(newnet)
  
  # rownames(cor1)<-paste(goodfeats_temp[,c(1)],goodfeats_temp[,c(2)],sep="_")
  #  colnames(cor1)<-rownames(cor1)
  
  cor_range<-round(range(cor1,na.rm=TRUE),2)
  
  mainlab1<-paste("Pairwise correlation of variables \n correlation range: ",cor_range[1]," to ",cor_range[2],sep="")
  
  # print(mainlab1)
  #xlab="Samples",ylab=""
  
  #dendrogram="none",
  if(do.clustering==TRUE){
  h1<-heatmap.2(cor1,col=rev(brewer.pal(11,"RdBu")),Rowv=TRUE,Colv=TRUE,scale="none",key=TRUE, symkey=FALSE, 
                density.info="none", trace="none",main=mainlab1,cexRow = 0.5,cexCol = 0.5,cex.main=0.7,
                xlab="Dataset B",ylab="Dataset A")
  }else{
    h1<-heatmap.2(cor1,col=rev(brewer.pal(11,"RdBu")),Rowv=FALSE,Colv=FALSE,scale="none",key=TRUE, symkey=FALSE, 
                  density.info="none", trace="none",main=mainlab1,cexRow = 0.5,cexCol = 0.5,cex.main=0.7,
                  xlab="Dataset B",ylab="Dataset A",dendrogram="none")
    
  }
  correlations_melted$from<-factor(correlations_melted$from,levels=unique(correlations_melted$from))
  correlations_melted$to<-factor(correlations_melted$to,levels=unique(correlations_melted$to))
  p<-get_bubbleplot(correlations_melted[,c(2,1,3)],statistic.type = "correlation",newdevice = FALSE,reverse.yaxis = TRUE)
  print(p)
  if(newdevice==TRUE){
    
      
      try(dev.off(),silent=TRUE)
  }
  return(list("cor.matrix"=cor1,"link.matrix"=correlations_melted))
  
}

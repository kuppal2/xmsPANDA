get_partial_cornet <-
function(data_m_fc_withfeats, sigfeats.index=NA,targeted.index=NA,networkscope="global",cor.method="spearman",abs.cor.thresh=0.4,
                             cor.fdrthresh=0.2,outloc,net_node_colors,net_legend=FALSE,netrandseed=555,pcor.method="pcor.shrink",newdevice=FALSE){
  
  suppressMessages(library(GeneNet))
  suppressMessages(library(igraph))
  
  setwd(outloc)
  l1<-list.files(".")
  
  data_m_fc_withfeats<-as.data.frame(data_m_fc_withfeats)
  
  metab_names<-paste(data_m_fc_withfeats$mz,data_m_fc_withfeats$time,sep="_")
  
  
  
  if(pcor.method=="cor2pcor"){
    
    cormat<-get_full_cormat(data_m_fc_withfeats,cor.method=cor.method)
    
    
    
    fname<-paste("complete_correlation_matrix.txt",sep="")
    
    cormat<-round(cormat,2)
    
    #write.table(cormat,file=fname,sep="\t",row.names=TRUE)
    
    
    
    p1<-cor2pcor(cormat)
  }else{
    
    p1<-suppressWarnings(pcor.shrink(data_mt))
    
    p1<-replace(p1,which(p1>1),1)
    p1<-replace(p1,which(p1<(-1)),-1)
  }
  colnames(p1)<-as.character(metab_names)
  rownames(p1)<-as.character(metab_names)
  
  suppressWarnings(dir.create("Tables",showWarnings = FALSE))
  p1<-round(p1,2)
  ####savep1,file="partial_cor.Rda")
  write.table(p1,file="Tables/partial_cor.txt",sep="\t",row.names=TRUE)
  
  
  
  colnames(p1)<-as.character(metab_names)
  rownames(p1)<-as.character(metab_names)
  
  
  if(is.na(targeted.index[1])==FALSE){
    
    networkscope="targeted"
  }
  
  
  
  
  #only raw p-values
  if(is.na(cor.fdrthresh)==TRUE){
    
    #p1sig<-cor0.test(p1, kappa=16)
    
    fname<-paste("net3_corpval",cor.fdrthresh,".Rda",sep="")
    
    edge.list<-network.test.edges(p1,fdr=FALSE,verbose=FALSE,plot=FALSE)
    
    
    net<-extract.network(edge.list,cutoff.ggm=0.80,verbose=FALSE)
    
    net2<-net[order(net$node1),]
    
    net3<-net2[,c(2,3,1,4)]
    
    
    ####savenet3,file=fname)
    
    
    
    
    
    
    
    netcormat<-matrix(data=0,nrow=dim(p1)[2],ncol=dim(p1)[2])
    colnames(netcormat)<-colnames(p1)
    rownames(netcormat)<-rownames(p1)
    
    for(i in 1:(dim(net3)[1])){
      
      r<-net3$node1[i]
      c<-net3$node2[i]
      p<-net3$pcor[i]
      
      netcormat[r,c]<-p
      netcormat[c,r]<-p
    }
    diag(netcormat)<-1
    #print(dim(netcormat))
    if(is.na(targeted.index[1])==TRUE){
      
      targeted.index=seq(1,dim(netcormat)[2])
    }
    
    if(is.na(sigfeats.index[1])==FALSE){
      netcormat<-netcormat[c(sigfeats.index),c(targeted.index)]
      
    }else{
      netcormat<-netcormat[,c(targeted.index)]
    }
    
  }else{
    
    
    if(cor.fdrthresh!=(-1) | (is.na(cor.fdrthresh)==TRUE)){
      
      #p1sig<-cor0.test(p1, kappa=16)
      fname<-paste("net3_corpval",cor.fdrthresh,".Rda",sep="")
      
      edge.list<-network.test.edges(p1,fdr=TRUE,verbose=FALSE,plot=FALSE)
      
      
      net<-extract.network(edge.list,method.ggm="qval",cutoff.ggm=(1-cor.fdrthresh),verbose=FALSE)
      
      net2<-net[order(net$node1),]
      
      net3<-net2[,c(2,3,1,4,5)]
      
      
      ####savenet3,file=fname)
      
      
      
      netcormat<-matrix(data=0,nrow=dim(p1)[2],ncol=dim(p1)[2])
      colnames(netcormat)<-colnames(p1)
      rownames(netcormat)<-rownames(p1)
      
      for(i in 1:(dim(net3)[1])){
        netcormat[net3$node1[i],net3$node2[i]]=net3$pcor[i]
        netcormat[net3$node2[i],net3$node1[i]]=net3$pcor[i]
      }
      diag(netcormat)<-1
      if(is.na(targeted.index)==TRUE){
        
        targeted.index=seq(1,dim(netcormat)[2])
      }
      if(is.na(sigfeats.index)==FALSE){
        netcormat<-netcormat[c(sigfeats.index),c(targeted.index)]
        
      }else{
        netcormat<-netcormat[,c(targeted.index)]
      }
    }else{
      if(is.na(sigfeats.index)==FALSE){
        netcormat<-p1[c(sigfeats.index),c(targeted.index)]
        
      }else{
        netcormat<-p1[,c(targeted.index)]
      }
    }
    
    
    
    
  }
  
  if(length(sigfeats.index)<dim(data_m_fc_withfeats)[1]){
    fname<-paste("Tables/significant_correlations_for_network.txt",sep="")
    
    # netcormat<-round(netcormat,2)
    
    write.table(netcormat,file=fname,sep="\t",row.names=TRUE)
    
    net_fname<-paste("Figures/partial_corsig",cor.fdrthresh,"network.tiff",sep="")
    
    max_cor_check<-max(abs(netcormat),na.rm=TRUE,warnings=FALSE)
    
    #if(max_cor_check>0)
    
    if(max_cor_check>=abs.cor.thresh){
      
      netcormat<-t(netcormat)
      colnames(netcormat)<-paste("Y",seq(1,dim(netcormat)[2]),sep="")
      rownames(netcormat)<-paste("X",seq(1,dim(netcormat)[1]),sep="")
      #tiff(net_fname, width=plots.width,height=plots.height,res=plots.res, compression="lzw")
      
      if(newdevice==TRUE){
      pdf(net_fname,width=8,height=10)
      }
      par_rows=1
      par(mfrow=c(par_rows,1))
      set.seed(netrandseed)
      net_result<-try(network(mat=as.matrix(nonsig_vs_fdr0.05_pearson_mat_bool_filt_1), cutoff=abs.cor.thresh,color.node = net_node_colors,
                              shape.node = c("rectangle", "circle"),
                              color.edge = c("red", "blue"), lwd.edge = 1,
                              show.edge.labels = FALSE, interactive = FALSE,cex.node.name=0.6),silent=TRUE)
      
      if(is(net_result,"try-error")){
        
        if(FALSE){
          set.seed(netrandseed)
          net_result<-try(network(mat=as.matrix(netcormat), threshold=abs.cor.thresh,color.node = net_node_colors,shape.node = c("rectangle", "circle"),
                                  color.edge = c("red", "blue"), lwd.edge = 1,show.edge.labels = FALSE, interactive = FALSE,cex.node.name=0.6),silent=TRUE)
        }
      }
      
      #legend("bottomright",c("Row #","Colum #"),pch=c(22,19),col=net_node_colors, cex=0.8)
      print(net_result)
      
      if(net_legend==TRUE){
        (legend("bottomright",c("Row #","Colum #"),pch=c(22,19),col=net_node_colors, cex=cex.plots,title="Network matrix values:"))
      }
      
      
      # ###savenet_result,"metabnet.Rda")
      
      
      write.graph(net_result$gR, file = "Figures/network_cytoscape_format.gml", format = "gml")
      
      if(newdevice==TRUE){
        
      try(dev.off(),silent=TRUE)
      }
    }else{
      print("No significant correlations found.")
    }
  }else{
    
    print("Network plot can not be generated. Rows and columns need to be exclusive.")
  }
  net4<-netcormat
  

  
  return(net4)
}

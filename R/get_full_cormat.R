get_full_cormat <-
function(data_m_fc_withfeats,targeted.index=NA,cor.method="spearman",write.files=TRUE){
  Sys.sleep(0.2)
  
  allsig_pcornetwork<-{}
  cormat<-{}
  
  allsig_pcornetwork_fdr0.05<-{}
  
  metab_names<-paste(data_m_fc_withfeats$mz,data_m_fc_withfeats$time,sep="_")
  rnames<-paste("mzid_",seq(1,dim(data_m_fc_withfeats)[1]),sep="")
  
  if(dim(data_m_fc_withfeats)[1]>1)
  {
    m1<-apply(data_m_fc_withfeats[,-c(1:2)],2,as.numeric)
    
    
    rownames(m1)=rnames
    
    data_mt<-t(m1)
  }else{
    
    m1<-as.numeric(data_m_fc_withfeats[,-c(1:2)])
    
    
    data_mt<-(m1)
  }
  data_mt<-as.matrix((data_mt))
  
  print(paste("Computing ",cor.method," correlation matrix",sep="")) 
  cormat<-WGCNA::cor(data_mt,use="pairwise.complete.obs",method=cor.method)
  colnames(cormat)<-as.character(metab_names)
  rownames(cormat)<-as.character(metab_names)
  
  num_samp<-dim(data_m_fc_withfeats[,-c(1:2)])[2]
  pearson_resvec<-as.vector(cormat)
  
  cl<-makeCluster(3)
  
  pearson_pvalue<-parLapply(cl,1:length(pearson_resvec),function(x){
    
    return(WGCNA::corPvalueStudent(pearson_resvec[x],num_samp))
    
    
  })
  
  stopCluster(cl)
  
  complete_pearsonpvalue_mat<-unlist(pearson_pvalue)
  dim(complete_pearsonpvalue_mat)<-dim(cormat)
  
  
  cormat<-as.data.frame(cormat)
  cormat<-as.matrix(cormat)
  pearson_Res_all<-{}
  
  cormat<-round(cormat,2)
  
  complete_pearsonpvalue_mat<-as.data.frame(complete_pearsonpvalue_mat)
  
  complete_pearsonpvalue_mat<-round(complete_pearsonpvalue_mat,2)
  
  if(write.files==TRUE)
  {
    write.table(cormat,file="correlationmatrix.txt",sep="\t",row.names=TRUE)
    
    write.table(complete_pearsonpvalue_mat,file="pvalue_matrix.txt",sep="\t",row.names=TRUE)
    
    ####savecormat,file="full_correlation.Rda")
    
    if(is.na(targeted.index)==FALSE){
      
      cormat<-cormat[targeted.index,]
      complete_pearsonpvalue_mat<-complete_pearsonpvalue_mat[targeted.index,]
      
      write.table(cormat,file="targetcorrelationmatrix.txt",sep="\t",row.names=TRUE)
      write.table(complete_pearsonpvalue_mat,file="targetpvalue_matrix.txt",sep="\t",row.names=TRUE)
      
    }
  }
  
  return(cormat)
}

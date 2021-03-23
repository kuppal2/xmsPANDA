get_fdr_adjusted_pvalue <-
function(data_matrix,filename=NA,fdrmethod="BH")
{
  
  toMatch=c("P.value","p.value","pvalue","Pvalue","p-value","P-value")
  
  
  if(is.na(filename)==FALSE){
    f1<-read.table(filename,sep="\t",header=TRUE)
    
    pval_columns <- unique (grep(paste(toMatch,collapse="|"),colnames(f1),value=FALSE))
    
    cnames1<-colnames(f1)
    
    pvalues_matrix<-as.data.frame(f1[,pval_columns])
    
    colnames(pvalues_matrix)<-cnames1[pval_columns]
  }else{
    
    pval_columns <- unique (grep(paste(toMatch,collapse="|"),colnames(data_matrix),value=FALSE))
    
    cnames1<-colnames(data_matrix)
    
    pvalues_matrix<-as.data.frame(data_matrix[,pval_columns])
    
    colnames(pvalues_matrix)<-cnames1[pval_columns]
    
    
  }
  
  
  if(is.vector(pvalues_matrix)==TRUE){
    
    pvalues_matrix<-as.data.frame(pvalues_matrix)
    
    fdr_res<-get_fdr_adjusted_pvalue_child(pvalues_matrix[,1],fdrmethod)
    
  }else{
    
    if(ncol(pvalues_matrix)>=1){
      
      fdr_res<-apply(pvalues_matrix,2,get_fdr_adjusted_pvalue_child,fdrmethod)
      
    }
    
    
  }
  cnames1<-colnames(fdr_res)
  
  cnames1<-gsub(cnames1,pattern=paste(toMatch,collapse="|"),replacement="FDR.adjusted.pvalue")
  
  colnames(fdr_res)<-cnames1
  
  return(fdr_res)
}

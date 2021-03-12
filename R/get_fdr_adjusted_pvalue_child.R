get_fdr_adjusted_pvalue_child <-
function(pvalues_vec,fdrmethod="BH")
{
  
  pvalues1=pvalues_vec
  if(fdrmethod=="none"){
    fdr_adjust_pvalue1<-pvalues1
    
  }
  
  if(fdrmethod=="BH"){
    fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BH")
    
  }else{
    if(fdrmethod=="ST"){
      
      fdr_adjust_pvalue1<-try(qvalue(pvalues1),silent=TRUE)
      
      if(is(fdr_adjust_pvalue1,"try-error")){
        
        fdr_adjust_pvalue1<-qvalue(pvalues1,lambda=max(pvalues1,na.rm=TRUE))
      }
      fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qvalues
      
    }else{
      if(fdrmethod=="Strimmer"){
        pdf("fdrtool.pdf")
        #par_rows=1
        #par(mfrow=c(par_rows,1))
        fdr_adjust_pvalue1<-suppressWarnings(fdrtool(as.vector(pvalues1),statistic="pvalue",verbose=FALSE))
        fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qval
        try(dev.off(),silent=TRUE)
      }else{
        
        if(fdrmethod=="BY"){
          fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BY")
          
        }else{
          if(fdrmethod=="bonferroni"){
            fdr_adjust_pvalue1<-p.adjust(pvalues1,method="bonferroni")
          }
        }
        
      }
    }
  }
  
  
  return(fdr_adjust_pvalue1)
  
}

runlmreg <-
function(X,Y,fdrmethod="BH",fdrthresh=0.05,pvalue.thresh=0.05,logistic_reg=FALSE,poisson_reg=FALSE,robust.estimate=FALSE){
  options(warn=-1)
  data_m_fc<-X[,-c(1:2)]
  
  classlabels_response_mat<-Y
  data_m_fc_withfeats<-X
  
  poisson_reg=FALSE
  rm(X)
  
  if(logistic_reg==FALSE){
    fileheader="lmreg"
  }else{
    fileheader="logitreg"
    
  }
  
  
  res1<-apply(data_m_fc,1,function(x){
    xvec<-x
    
    data_mat_anova<-cbind(xvec,classlabels_response_mat)
    
    cnames<-colnames(data_mat_anova)
    cnames[1]<-"Response"
    
    colnames(data_mat_anova)<-cnames
    
    
    
    anova_res<-diffexplmreg(dataA=data_mat_anova,logistic_reg=logistic_reg,poisson_reg=poisson_reg,robust.estimate=robust.estimate)
    
    return(anova_res)
  })
  
  
  
  main_pval_mat<-{}
  
  posthoc_pval_mat<-{}
  pvalues<-{}
  
  all_inf_mat<-{}
  
  #for(i in 1:length(res1))
  all_inf_mat<-lapply(1:length(res1),function(i)
  {
    
    main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
    pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
    
    cur_pvals<-t(res1[[i]]$mainpvalues)
    cur_est<-t(res1[[i]]$estimates)
    cur_stderr<-t(res1[[i]]$stderr)
    cur_tstat<-t(res1[[i]]$statistic)
    cur_res<-cbind(cur_pvals,cur_est,cur_stderr,cur_tstat)
    
    
    
    #all_inf_mat<-rbind(all_inf_mat,cur_res)
    
    return(cur_res)
    
  })
  
  cur_pvals<-t(res1[[1]]$mainpvalues)
  
  all_inf_mat<-do.call(rbind,all_inf_mat)
  
  pvalues=all_inf_mat[,1]
  
  pvalues<-replace(pvalues,which(is.na(pvalues)==TRUE),1)
  all_inf_mat[,1]=pvalues
  
  #cnames_1<-c(paste(colnames(cur_pvals),".pvalue",sep=""),paste(colnames(cur_pvals),".estimate",sep=""),paste(colnames(cur_pvals),".tstatistic",sep=""))
  
  cnames_1<-c(paste("P.value_",colnames(cur_pvals),sep=""),paste("Estimate_",colnames(cur_pvals),sep=""),paste("StdError_var_",colnames(cur_pvals),sep=""),paste("t-statistic_",colnames(cur_pvals),sep=""))
  
  # ##save(cnames_1,file="c1.Rda")
  ###save(data_allinf_withfeats,file="dA.Rda")
  ###save(all_inf_mat,file="dB.Rda")
  
  data_allinf_withfeats<-cbind(all_inf_mat,data_m_fc_withfeats)
  
  colnames(data_allinf_withfeats)<-c(cnames_1,colnames(data_m_fc_withfeats))
  
  if(fdrmethod=="BH"){
    fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
  }else{
    if(fdrmethod=="ST"){
      #fdr_adjust_pvalue<-qvalue(pvalues)
      #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
      
      
      fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
      
      if(is(fdr_adjust_pvalue,"try-error")){
        
        fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
      }
      
      fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
      
      
      
    }else{
      if(fdrmethod=="Strimmer"){
        pdf("fdrtool.pdf")
        
        fdr_adjust_pvalue<-suppressWarnings(fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE))
        fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
        try(dev.off(),silent=TRUE)
      }else{
        if(fdrmethod=="none"){
          #fdr_adjust_pvalue<-pvalues
          fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
        }else{
          if(fdrmethod=="BY"){
            fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
          }else{
            if(fdrmethod=="bonferroni"){
              fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
              
            }
          }
        }
      }
    }
    
    
    
  }
  
  
  
  cnames_tab<-colnames(data_m_fc_withfeats)
  #cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
  
  pvalues<-as.data.frame(pvalues)
  
  final.pvalues<-pvalues
  
  sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
  
  #data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
  
  #colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
  
  
  filename<-paste(fileheader,"results.allfeatures.txt",sep="")
  
  #data_allinf_withfeats<-cbind(fdr_adjust_pvalue,all_inf_mat,data_m_fc_withfeats)
  
  pval_columns<-grep(colnames(data_allinf_withfeats),pattern="P.value")
  
  fdr_adjusted_pvalue<-get_fdr_adjusted_pvalue(data_matrix=data_allinf_withfeats,fdrmethod=fdrmethod)
  
  #   data_allinf_withfeats1<-cbind(data_allinf_withfeats[,pval_columns],fdr_adjusted_pvalue,data_allinf_withfeats[,-c(pval_columns)])
  
  cnames_tab1<-c(cnames_tab[pval_columns],colnames(fdr_adjusted_pvalue),cnames_tab[-pval_columns])
  pval_columns<-grep(colnames(data_allinf_withfeats),pattern="P.value")
  
  fdr_adjusted_pvalue<-get_fdr_adjusted_pvalue(data_matrix=data_allinf_withfeats,fdrmethod=fdrmethod)
  
  data_allinf_withfeats<-cbind(data_allinf_withfeats[,pval_columns],fdr_adjusted_pvalue,data_allinf_withfeats[,-c(pval_columns)])
  
  cnames_tab1<-c(cnames_tab[pval_columns],colnames(fdr_adjusted_pvalue),cnames_tab[-pval_columns])
  
  
  
  write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)
  filename<-paste(fileheader,"results.selectedfeatures.txt",sep="")
  
  data_select_withfeats<-data_allinf_withfeats[sel.diffdrthresh,]
  
  write.table(data_select_withfeats, file=filename,sep="\t",row.names=FALSE)
  
  
  cnames_tab<-colnames(data_m_fc_withfeats)
  
  classlabels_response_mat<-as.data.frame(classlabels_response_mat)
  
  class_column_names<-colnames(classlabels_response_mat)
  
  # cnames_tab<-c(paste("P.value_var",1:dim(classlabels_response_mat)[2],sep=""),
  #paste("Estimate_var",1:dim(classlabels_response_mat)[2],sep=""), paste("StdError_var",1:dim(classlabels_response_mat)[2],sep=""),
  #paste("statistic_var",1:dim(classlabels_response_mat)[2],sep=""),cnames_tab)
  
  #cnames_1<-c(paste("P.value_",colnames(cur_pvals),sep=""),paste("Estimate_",colnames(cur_pvals),sep=""),paste("StdError_var_",colnames(cur_pvals),sep=""),paste("t-statistic_",colnames(cur_pvals),sep=""))
  
  
  if(ncol(classlabels_response_mat)>1){
    cnames_1<-c(paste("FDR.adjusted.P.value_",class_column_names,sep=""),paste("P.value_",class_column_names,sep=""),
                paste("Estimate_",class_column_names,sep=""), paste("StdError_",class_column_names,sep=""),
                paste("t-statistic_",class_column_names,sep=""))
  }else{
    cnames_1<-c(paste("FDR.adjusted.P.value",sep=""),paste("P.value",sep=""),
                paste("Estimate",sep=""), paste("StdError",sep=""),
                paste("t-statistic",sep=""))
    
  }
  
  
  # cnames_1<-c(paste("P.value_",colnames(cur_pvals),sep=""),paste("Estimate_",colnames(cur_pvals),sep=""),paste("StdError_",colnames(cur_pvals),sep=""),paste("t-statistic_",colnames(cur_pvals),sep=""))
  
  cnames_tab<-c(cnames_1,cnames_tab)
  
  
  colnames(data_allinf_withfeats)<-as.character(cnames_tab)
  
  write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)
  
  #data_allinf_withfeats<-data_allinf_withfeats[order(data_allinf_withfeats[,2]),]
  
  sel.diffdrthresh=which(sel.diffdrthresh==TRUE)
  sel.diffdrthresh<-sel.diffdrthresh[order(data_allinf_withfeats[sel.diffdrthresh,2])]
  
  options(warn=0)
  return(list("all"=data_allinf_withfeats,"selected"=data_select_withfeats,"selected.index"=sel.diffdrthresh))
  
}

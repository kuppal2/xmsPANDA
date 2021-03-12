do_eigenms_norm <-
function(data_m,classlabels,featselmethod=NA,feature_id_vector=NA,pairedanalysis=FALSE){
  
  cnames_1<-colnames(classlabels)
  
  if(dim(classlabels)[2]>2){
    
    if(is.na(featselmethod)==FALSE){
      if(featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="lm2wayanova" | featselmethod=="lm2wayanovarepeat"){
        
        
        
        if(pairedanalysis==FALSE){
          treatment_vec<-as.factor(paste(classlabels[,2],":",classlabels[,3],sep=""))
        }else{
          treatment_vec<-as.factor(paste(classlabels[,3],":",classlabels[,4],sep=""))
          
        }
        
        
        
      }else{
        
        if(featselmethod=="limma1way" | featselmethod=="limma1wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm1wayanovarepeat"){
          check_factor2<-grep(cnames_1,pattern="Factor1")
          
          
          
          if(pairedanalysis==FALSE){
            treatment_vec<-as.factor(classlabels[,2])
          }else{
            treatment_vec<-as.factor(classlabels[,3])
            
          }
          
        }else{
          
          treatment_vec<-as.factor(classlabels[,2])
        }
        
        
        #treatment_vec<-as.factor(classlabels[,2])
        
      }
      
    }else{
      
      
      treatment_vec<-as.factor(classlabels[,2])
      
      
      
    }
    
    
  }else{
    
    treatment_vec<-as.factor(classlabels[,2])
  }
  edata<-data_m #log2(data_m+1)
  edata<-as.matrix(edata)
  
  #save(edata,treatment_vec,feature_id_vector,file="testeig.Rda")
  set.seed(123)
  
  if(is.na(feature_id_vector)==FALSE){
    ints_eig1 = eig_norm1(m=edata, treatment=treatment_vec, prot.info=cbind(paste('ID_',seq(1:nrow(data_m)), sep=''),feature_id_vector))
  }else{
    ints_eig1 = eig_norm1(m=edata, treatment=treatment_vec, prot.info=cbind(paste('ID_',seq(1:nrow(data_m)), sep=''),seq(1:nrow(data_m))))
    
  }
  
  set.seed(123)
  data_m = eig_norm2(rv=ints_eig1)
  data_m = data_m$norm_m
  return(data_m)
  
}

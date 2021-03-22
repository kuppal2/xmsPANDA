do_sva_norm <-
function(data_m,classlabels,featselmethod=NA){
  
  
  edata<-(data_m)
  
  pdata<-classlabels
  
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
          #if(length(check_factor2)>0)
          
          
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
  
  classlabels<-as.matrix(classlabels)
  
  trainMod = model.matrix(~treatment_vec)
  trainMod0 = model.matrix(~1,data=pdata)
  
  edata<-as.matrix(edata)
  
  
  trainSv = sva(edata,trainMod,trainMod0)
  
  fsvaobj = fsva(dbdat=edata,mod=trainMod,sv=trainSv,newdat=edata)
  
  data_m<-(fsvaobj$db)
  
  return(data_m)
  
}

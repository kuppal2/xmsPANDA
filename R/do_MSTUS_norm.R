do_MSTUS_norm <-
function(data_m,missing.val=0){
  
  total_sigs<-apply(data_m,1,function(x){
    if(is.na(missing.val)==FALSE){return(length(which(x>missing.val)))
    }else{
      return(length(which(is.na(x)==FALSE)))
    }})
  
  
  useful_metabs_1<-which(total_sigs>=dim(data_m)[2])
  
  useful_metabs_2<-which(total_sigs>=dim(data_m)[2]*0.95)
  
  if(length(useful_metabs_1)>0){
    
    data_useful<-data_m[useful_metabs_1,]
    
    
  }else{
    
    if(length(useful_metabs_2)>0){
      data_useful<-data_m[useful_metabs_2,]
      
    }else{
      useful_metabs_3<-which(total_sigs>=dim(data_m)[2]*0.9)
      
      if(length(useful_metabs_3)>0){
        
        data_useful<-data_m[useful_metabs_3,]
      }else{
        
        useful_metabs_4<-which(total_sigs>=dim(data_m)[2]*0.8)
        
        if(length(useful_metabs_4)>0){
          
          data_useful<-data_m[useful_metabs_4,]
        }
        
      }
    }
  }
  median_int<-apply(data_m,2,function(x){median(x,na.rm=TRUE)})
  sum_int<-apply(data_useful,2,function(x){sum(x,na.rm=TRUE)})
  
  data_m<-sweep(data_m,2,(sum_int/median_int),'/')
  
  return(data_m)
  
  
}

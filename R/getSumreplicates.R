getSumreplicates <-
function(curdata,alignment.tool,numreplicates,numcluster,rep.max.missing.thresh,summary.method="mean",
                           summary.na.replacement="zeros",missing.val=0,col_end=NA)
{
  mean_replicate_difference<-{}
  sd_range_duplicate_pairs<-{}
  #print(alignment.tool)
  if(is.na(alignment.tool)==FALSE){
  if(alignment.tool=="apLCMS")
  {
    col_end=2
  }
  else
  {
    if(alignment.tool=="XCMS")
    {
      col_end=2
    }
    else
    {
      #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
    }
  }
  }
  
  if(is.na(col_end)==FALSE){
  curdata_mz_rt_info=curdata[,c(1:col_end)]
  curdata=curdata[,-c(1:col_end)]
  }else{
    curdata_mz_rt_info=curdata
    
  }
  
  
  cl<-parallel::makeCluster(numcluster)
  numfeats=dim(curdata)[1]
  numsamp=dim(curdata)[2]	
  
  if(numsamp%%numreplicates>0){
    
    stop("Not all samples have replicates or the numreplicates value is not correct.")
  }
  
  clusterEvalQ(cl, "getSumreplicateschild")
  sub_samp_list<-list()
  
  sampcount=1
  for(samp in seq(1,(numsamp),numreplicates))
  {
    i=samp
    j=i+numreplicates-1
    if(dim(curdata[,c(i:j)])[1]>0){
      sub_samp_list[[sampcount]]=curdata[,c(i:j)]
    }
    sampcount=sampcount+1
  }
  
  avg.res<-parSapply(cl,sub_samp_list,getSumreplicateschild,alignment.tool=alignment.tool,numreplicates=numreplicates,rep.max.missing.thresh=rep.max.missing.thresh,method=summary.method,missing.val=missing.val)
  #avg.res<-getAvgreplicateschild(sub_samp_list[[1]],alignment.tool,numreplicates)
  #print("done")
  
  
  stopCluster(cl)
  
  
  
  final_set<-as.data.frame(avg.res)
  colnames_data<-colnames(curdata)
  colnames_data<-colnames_data[seq(1,(numsamp),numreplicates)]
  colnames(final_set)<-colnames_data
  rownames(final_set)=NULL
  #final_set<-cbind(curdata_mz_rt_info,final_set)
  
  final_set<-apply(final_set,2,as.numeric)
  #	write.table(final_set,file="final_Set.txt",sep="\t",row.names=FALSE)
  
  if(summary.na.replacement=="zeros"){
    
    final_set<-replace(final_set,which(is.na(final_set)==TRUE),0)
  }else{
    if(summary.na.replacement=="halfsamplemin"){
      
      
      final_set<-apply(final_set,2,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
    }else{
      
      if(summary.na.replacement=="halfdatamin"){
        
        
        min_val<-min(final_set,na.rm=TRUE)*0.5
        final_set<-replace(final_set,which(is.na(final_set)==TRUE),min_val)
        
        #data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
      }else{
        if(summary.na.replacement=="halffeaturemin"){
          
          
          final_set<-apply(final_set,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
          final_set<-t(final_set)
        }else{
          
          if(summary.na.replacement=="knn"){
            
            #	##save(final_set,file="final_set1.Rda")
            # print("imputing")
            suppressMessages(library(impute))
            final_set<-impute.knn(final_set,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
            
            final_set<-final_set$data
            final_set<-as.data.frame(final_set)
            
          }else{
            
            if(summary.na.replacement=="randomforest"){
              
              final_set<-RFimpute(t(data_m))
              final_set<-t(final_set$ximp)
              final_set<-as.data.frame(final_set)
              
            }else{
              
              if(summary.na.replacement=="QRILC"){
                
                
                library(tmvtnorm)
                final_set<-QRILCimpute(data_m) #rows: features; cols: samples
                ##save(final_set,file="final_set.Rda")
                final_set<-ldply(final_set,rbind) 
                final_set<-t(final_set)               
                final_set<-as.data.frame(final_set)
                
                if(length(which(final_set<0))>0){
                  final_set<-replace(as.matrix(final_set),which(final_set<0),0)
                }
                
              }
              
              
            }
            
          }
          # ##save(final_set,file="final_set.Rda")
          
        }
      }
    }
    
    
  }
  final_set<-as.data.frame(final_set)
  return(final_set)
}

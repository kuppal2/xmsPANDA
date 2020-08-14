find.Overlapping <-
function (dataA, dataB, mz.thresh = 10, time.thresh = 30, best_one=FALSE)
{
  data_a <- as.data.frame(dataA)
  data_b <- as.data.frame(dataB)
  rm(dataA)
  rm(dataB)
  
  colnames(data_a)[1] = "mz"
  colnames(data_b)[1] = "mz"
  if (is.na(time.thresh) == FALSE) {
    colnames(data_a)[2] = "time"
    colnames(data_b)[2] = "time"
    mznames = c("index.A", "mz.data.A", "time.data.A", "index.B",
                "mz.data.B", "time.data.B", 'mz.difference(ppm)', "time.difference(s)")
    print("Using the 1st columns as 'mz' and 2nd columsn as 'retention time'")
  }else {
    mznames = c("index.A", "mz.data.A", "index.B", "mz.data.B", 'mz.difference(ppm)')
    print("Using the 1st columns as 'mz'")
  }
  
  mz_groups <- mclapply(1:dim(data_a)[1], function(j) {
    commat = {}
    
    ppmb = (mz.thresh) * (data_a$mz[j]/1e+06)
    getbind_same <- which(abs(data_b$mz - data_a$mz[j]) <= ppmb)
    if (is.na(time.thresh) == FALSE) {
      if (length(getbind_same) > 0) {
        
        tmp = suppressWarnings(cbind(cbind(j, data_a[j, c(1, 2)]),cbind(getbind_same, data_b[getbind_same, c(1, 2)])))
        tmp[,"mzdiff"]=(abs(tmp[,2]-tmp[,5])/tmp[,2])*10^6
        tmp[,"timediff"]=abs(tmp[,3]-tmp[,6])
        colnames(tmp)=mznames
        tmp = tmp[tmp$time.difference<=time.thresh,]
        if(best_one==TRUE){
          if(dim(tmp)[1]>0){
            commat=tmp[which(tmp$time.difference==min(tmp$time.difference)),]
          }
        }else{
          commat = tmp
        }
        
      }
      
    } else {
      if (length(getbind_same) > 0) {
        
        tmp = suppressWarnings(cbind(as.data.frame(cbind(j, data_a[j, c(1)])),as.data.frame(cbind(getbind_same, data_b[getbind_same, 1]))))
        tmp[,"mzdiff"]=(abs(tmp[,2]-tmp[,4])/tmp[,2])*10^6
        colnames(tmp)=mznames
        commat = tmp
        
      }
    }
    return(as.data.frame(commat))
  },mc.cores=detectCores())
  
  commat = data.frame()
  if (length(mz_groups) > 0) {
    
    for(j in 1:length(mz_groups)){
      if(dim(mz_groups[[j]])[1]>0){
        commat <- rbind(commat,mz_groups[[j]])
      }
    }
    if(dim(commat)[1]>0){
      rownames(commat)=1:dim(commat)[1]
    }else{
      
    }
    return(commat)
  }else{
    return(commat)
  }
  
}

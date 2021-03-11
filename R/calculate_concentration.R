calculate_concentration <-
function(targeted_table,ref,seq,outloc,missing= 0.3,batch_qstd_missing=0.4){
  
  count=max(grep("_error",colnames(targeted_table)))
  
  #summarize the number of batches
  num_batches<- length(unique(seq[,4]))
  
  #create list object with each batch in slot
  qstd_sample_mapping<- seq[grepl('qstd',seq[,3],ignore.case = TRUE),]
  qstd_int<- targeted_table[,as.character(qstd_sample_mapping[,1])]
  
  
  #show error if number of reference samples is not 6x number of batches
  #if(ncol(qstd_int)/6 != length(unique(sample.mapping[,4]))){
  #    stop("Number of reference samples is not multiple of number of batches", call.=TRUE)
  #}
  
  print(paste("Number of batches: ",length(unique(seq[,4])), sep= ""))
  write.table(paste("Number of batches: ",length(unique(seq[,4])), sep= ""), paste(outloc,"step2_summary.txt",sep="/"), sep="\t", row.names=FALSE, col.names =FALSE)
  print(paste("Number of reference standards(QSTD samples): ",ncol(qstd_int), paste= ""))
  write.table(paste("Number of reference standards: ",ncol(qstd_int), sep= ""), paste(outloc,"step2_summary.txt",sep="/"), sep="\t", append = TRUE, row.names=FALSE, col.names =FALSE)
  
  
  #calculating response factor
  #create list of reference sample response factors for targeted, q-standard table
  rf_ls<- c()
  for(ii in 1:num_batches){
    batch_qs<- qstd_sample_mapping[qstd_sample_mapping[,4]==ii,]
    qstd_int_batch <- qstd_int[,as.character(batch_qs[,1])]
    qstd_int_batch[qstd_int_batch==0]<-NA
    pass <- apply(is.na(qstd_int_batch),1,sum)/dim(qstd_int_batch)[2]<batch_qstd_missing
    tmp <- rowMeans(qstd_int_batch,na.rm=TRUE)
    tmp[!(pass)]<-0
    rf_ls<-cbind(rf_ls, tmp)
  }
  
  
  ######restart here
  
  if(num_batches>1){
    #count number of zeroes in each row
    prop_missing<-rowSums(rf_ls==0)/num_batches
    #samples with number of missing values above proportion theshold (Qstd) and with no reference concentration
    remove_i<- unique(c(which(prop_missing> missing), which(is.na(targeted_table[,4]))))
    
    #updated data after removing missing values
    if(length(remove_i)!=0){
      rf_ls<-rf_ls[-remove_i,]
      targeted_table<- targeted_table[-remove_i,]
    }
    
    #replace missing values with group average (not including zeroes)
    rf_ls[rf_ls==0]<-NA
    #calculate means across study
    row_means<- rowMeans(rf_ls, na.rm=TRUE)
    
    #replace NA with row mean
    rf_ls_final<- rf_ls
    
    for(ii in 1:nrow(rf_ls)){
      sub_rf<- rf_ls[ii,]
      jj<- which(is.na(sub_rf), arr.ind=TRUE)
      if(length(jj)>0){
        rf_ls_final[ii,jj]<- row_means[ii]
      }
    }
    
    #evaluate CV
    cv_res<- rowSds(rf_ls_final)/rowMeans(rf_ls_final)*100
    print(paste("Average QSTD CV for all targets: ",round(mean(cv_res),2),"%", sep=""))
    write.table(paste("Average QSTD CV for all targets: ",round(mean(cv_res),2),"%", sep=""), paste(outloc,"step2_summary.txt",sep="/"), sep="\t", append = TRUE, row.names=FALSE, col.names =FALSE)
    
    pdf(paste(outloc,"CV_distribution_Qstd_intensity.pdf",sep="/"))
    hist(cv_res, breaks=40, xlab="CV", main="Distribution of QSTD Intensity CV")
    dev.off()
    
    rf_ls_final<-as.data.frame(rf_ls_final)
    colnames(rf_ls_final)=paste("Batch",seq(1:dim(rf_ls_final)[2]),sep="")
    CV<-cv_res
    qstd_avg_table<- cbind(targeted_table[,c(3,1,2)], CV, rf_ls_final)
    write.table(qstd_avg_table, paste(outloc,"Batchwise_Qstd_average_intensity.txt",sep="/"), sep= "\t", row.names=FALSE)
    #rf_df<- read.table("batchwise_response_factor.txt", sep= "\t", header=TRUE)
    
  }else{
    #samples with number of missing values above proportion theshold (Qstd) and with no reference concentration
    remove_i<- unique(c(which(is.na(targeted_table[,4]))))
    
    #updated data after removing missing values
    if(length(remove_i)!=0){
      rf_ls<-rf_ls[-remove_i,]
      targeted_table<- targeted_table[-remove_i,]
    }
    
    
    #remove all the targets with missing values
    rf_ls_final<- rf_ls[!rf_ls==0]
    targeted_table<- targeted_table[!rf_ls==0,]
    
    rf_ls_final<-as.data.frame(rf_ls_final)
    colnames(rf_ls_final)=paste("Batch",seq(1:dim(rf_ls_final)[2]),sep="")
    qstd_avg_table<- cbind(targeted_table[,c(3,1,2)], rf_ls_final)
    write.table(qstd_avg_table, paste(outloc,"Batchwise_Qstd_average_intensity.txt",sep="/"), sep= "\t", row.names=FALSE)
    
  }
  
  #create rf matrix
  rf_final<-{}
  for (ii in 1:nrow(rf_ls_final)){
    l1<- targeted_table[ii,4]/rf_ls_final[ii,]
    rf_final<- rbind(rf_final,l1)
  }
  
  rf_table<- cbind.data.frame(targeted_table[,c(3,1,2)], rf_final)
  write.table(rf_table, paste(outloc,"Batchwise_response_factors.txt",sep="/"), sep= "\t", row.names=FALSE)
  
  #calculate concentrations in each batch
  
  
  #includes q-standards
  res_concentrations<- data.frame(matrix(0, ncol=1, nrow=nrow(targeted_table)))
  for (ii in 1:num_batches){
    batch_feature_table<- targeted_table[, as.character(seq[seq[,4]==ii,1])]
    res_1<- c()
    for(jj in 1:nrow(batch_feature_table)){
      res_1<- rbind(res_1, rf_final[jj,ii]* batch_feature_table[jj,])
    }
    res_concentrations<- cbind(res_concentrations, res_1)
  }
  res_concentrations<- res_concentrations[,-1]
  
  final_concentrations<- cbind(targeted_table[1:(count+2)],res_concentrations)
  
  #create sample mapping
  sample_only_concentrations<- cbind(final_concentrations[,c(1:(count+2))],final_concentrations[,as.character(seq[grepl('study',seq[,3],ignore.case = TRUE),1])])
  sample_mean_concentration_batchwise<-final_concentrations[,c(1:(count+2))]
  
  for (ii in 1:num_batches){
    res_tmp <- final_concentrations[,as.character(seq[seq[,4]==ii & grepl('study',seq[,3],ignore.case = TRUE),1])]
    res_tmp[res_tmp==0] <- NA
    sample_mean_concentration_batchwise <- cbind(sample_mean_concentration_batchwise,apply(res_tmp,1,mean,na.rm=TRUE))
  }
  colnames(sample_mean_concentration_batchwise)[(count+3):((count+2)+num_batches)]=paste("Batch",seq(1:num_batches),sep="")
  sample_mean_concentration_batchwise[,(count+3):((count+2)+num_batches)][is.na(sample_mean_concentration_batchwise[,(count+3):((count+2)+num_batches)])] <- 0
  
  write.table(final_concentrations, paste(outloc,"Final_concentrations_with_Qstd_and_study.txt",sep="/"), sep= "\t", row.names=FALSE)
  write.table(sample_only_concentrations, paste(outloc,"Final_concentrations_with_study_only.txt",sep="/"), sep= "\t", row.names=FALSE)
  write.table(sample_mean_concentration_batchwise, paste(outloc,"Final_average_concentrations_with_study_only.txt",sep="/"), sep= "\t", row.names=FALSE)
  
  return(final_concentrations)
}

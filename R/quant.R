quant <-
function(Xmat=NA,Ymat=NA,Wmat=NA,Zmat=NA,feature_table,class_file,ref_list,foldchange_list,outloc,
                 num_replicates=1,
                 summarize_replicates=FALSE,
                 rep.max.missing.thresh=0.6,
                 summary.method="mean",
                 mass_error= 10,
                 time_error= 30,
                 percent_node=0.6,
                 foldchange_thresh=2,
                 steps="123",
                 min_num_nonmissing=3,
                 targetID=NA,
                 minhit=3,
                 groupcheck=TRUE,
                 highcolor='red',
                 lowcolor='blue',summary.na.replacement="zeros",missing.val=0,normalization.method="none",alphabetical.order=FALSE
) {
  theme_set(theme_classic())
  
  #read in data tables
  if(!is.na(feature_table)){
    feature.table<- read.table(feature_table, sep= "\t", header=TRUE)
  }else{
    if(is.data.frame(Xmat)){
      feature.table<-Xmat
    }else{
      stop("There is no feature table file. Please provide the path of feature table file or feature table R object.", call.=TRUE)
    }
  }
  if(is.numeric(feature.table[,1])==FALSE){
    stop('Please double check the mz column in your feature table, it is not numeric.')
  }
  if(is.numeric(feature.table[,2])==FALSE){
    stop('Please double check the time column in your feature table, it is not numeric.')
  }
  
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  
  #read in sample id mapping
  if(!is.na(class_file)){
    sample.mapping<- read.table(class_file, sep= "\t", header= TRUE)
  }else{
    if(is.data.frame(Ymat)){
      sample.mapping<-Ymat
    }else{
      stop("There is no class label file. Please provide the path of class label file or class label R object.", call.=FALSE)
    }
  }
  colnames(sample.mapping)[1:4]=c("FileName","SampleID","Sample_type","Batch")
  
  
  #read in reference standardization list
  if(!is.na(ref_list)){
    ref_table<- read.table(ref_list, sep= "\t", header= TRUE, quote= "\"")
  }else{
    if(is.data.frame(Wmat)){
      ref_table<-Wmat
    }else{
      stop("There is no reference standardization list. Please provide the path of reference standardization list or list R object.", call.=FALSE)
    }
  }
  colnames(ref_table)[c(1,2,3,4)]=c("ref_mz","ref_time","Chemical_name","Qstd3")
  if(is.numeric(ref_table$ref_mz)==FALSE){
    stop('Please double check the mz column in your reference standardization list, it is not numeric.')
  }
  if(is.numeric(ref_table$ref_time)==FALSE){
    stop('Please double check the time column in your reference standardization list, it is not numeric.')
  }
  if(is.numeric(ref_table$Qstd3)==FALSE){
    stop('Please double check the Qstd column in your reference standardization list, it is not numeric.')
  }
  
  if(length(grep("3",steps))!=0){
    #read in fold change list
    if(!is.na(foldchange_list)){
      foldchange_table<- read.table(foldchange_list, sep= "\t", header= TRUE, quote= "\"")
    }else{
      if(is.data.frame(Zmat)){
        foldchange_table<-Zmat
      }else{
        stop("There is no fold change list. Please provide the path of fold change list or list R object.", call.=FALSE)
      }
    }
  }
  
  #check if the file name can match with feature table
  match1 <- setdiff(colnames(feature.table[,-c(1,2)]),as.character(sample.mapping[,1]))
  match2 <- setdiff(as.character(sample.mapping[,1]),colnames(feature.table[,-c(1,2)]))
  if(length(match1)!=0 && length(match2)!=0){
    stop("The file name in class file can't be matched with the name in feature table. Please check if the fileName has '.' symbol. If so, please use '_' to replace it.", call.=FALSE)
  }else{
    feature.table=feature.table[,c("mz","time",as.character(sample.mapping[,1]))]
  }
  
  suppressWarnings(dir.create(outloc,showWarnings = FALSE))
  setwd(outloc)
  
  if(summarize_replicates==TRUE){
    capture.output(avg_feature_table<-data_preprocess(X=feature.table,feature_table_file=NA,parentoutput_dir=outloc,
                                                      
                                                      class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh=NA,summarize.replicates=TRUE,
                                                      summary.method=summary.method,all.missing.thresh=NA,group.missing.thresh=NA,
                                                      log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,
                                                      lowess_norm=FALSE,madscaling=FALSE,missing.val=missing.val,samplermindex=NA, 
                                                      rep.max.missing.thresh=rep.max.missing.thresh,summary.na.replacement=summary.na.replacement,
                                                      normalization.method = normalization.method,alphabetical.order=alphabetical.order), file='/dev/null')
    avg_feature_table<-avg_feature_table$data_matrix_afternorm_scaling
    unlink(paste(outloc,"Stage1",sep="/"), recursive = TRUE)
  }else{
    avg_feature_table<-feature.table
  }
  #create average sample mapfile
  sample.mapping<- sample.mapping[as.character(sample.mapping[,1])%in%colnames(avg_feature_table)[-c(1,2)],]
  
  ######################## pre-processing ############################
  
  #create feature table with reference masses
  #matches ref list massess with feature table
  capture.output(overlap<- find.Overlapping(ref_table[,c(1,2)], avg_feature_table, mz.thresh = mass_error, time.thresh = time_error), file='/dev/null')
  
  #show error if there is no match between feature table and reference list
  if(length(overlap)==0){
    stop("There is no match between feature table and reference list.",call.=FALSE)
  }
  
  #calculate ppm error between reference mass and feature table
  mz_error<- overlap$`mz.difference(ppm)`
  #merges detected metabolites in reference list with feature table; includes all samples
  if(is.null(overlap$time.difference)){
    if(length(grep("3",steps))!=0){
      foldchange=foldchange_table[overlap$index.B,3]
      r_targeted_table<- cbind(ref_table[overlap$index.A,], foldchange , mz_error, avg_feature_table[overlap$index.B,])
    }else{
      r_targeted_table<- cbind(ref_table[overlap$index.A,], mz_error, avg_feature_table[overlap$index.B,])
    }
  }else{
    if(length(grep("3",steps))!=0){
      foldchange=foldchange_table[overlap$index.B,3]
      r_targeted_table<- cbind(ref_table[overlap$index.A,], foldchange, mz_error, overlap$time.difference, avg_feature_table[overlap$index.B,])
      colnames(r_targeted_table)[match("overlap$time.difference",colnames(r_targeted_table))]='time_error'
    }else{
      r_targeted_table<- cbind(ref_table[overlap$index.A,], mz_error, overlap$time.difference, avg_feature_table[overlap$index.B,])
      colnames(r_targeted_table)[match("overlap$time.difference",colnames(r_targeted_table))]='time_error'
    }
  }
  
  ###################################################################
  
  #check distribution
  if(length(grep("1",steps))!=0){
    print("1. Plotting the intensity distribution of samples and qstd.")
    outloc1=paste(outloc,"/step1",sep="")
    suppressWarnings(dir.create(outloc1,showWarnings = FALSE))
    generate_distribution(targeted_table=r_targeted_table,seq=sample.mapping,outloc=outloc1,groupcheck=groupcheck,targetID=targetID,min_num_nonmissing=min_num_nonmissing)
    print("Step 1 complete")
  }
  
  #calculate the concentration
  if(length(grep("2",steps))!=0){
    print("2. Calculating the concentration for all samples.")
    outloc2=paste(outloc,"/step2",sep="")
    suppressWarnings(dir.create(outloc2,showWarnings = FALSE))
    final_study_qstd_conc=calculate_concentration(targeted_table=r_targeted_table,ref=ref_table,seq=sample.mapping,outloc=outloc2)
    
    generate_distribution_concentration(targeted_table=final_study_qstd_conc,seq=sample.mapping,outloc=outloc2,groupcheck=groupcheck,targetID=targetID,min_num_nonmissing=min_num_nonmissing)
    
    print("Step 2 complete")
  }
  
  #draw kegg map
  if(length(grep("3",steps))!=0){
    print("3. Pulling the KEGG maps out from KEGG website.")
    outloc3=paste(outloc,"/step3",sep="")
    suppressWarnings(dir.create(outloc3,showWarnings = FALSE))
    draw_keggmap(targeted_table=r_targeted_table,outloc=outloc3,foldchange.thresh=foldchange_thresh,minhit=minhit,highcolor=highcolor,lowcolor=lowcolor,percent_node=percent_node)
    print("Step 3 complete")
  }
  
  #suppressWarnings(dir.create(paste(outloc,"KEGGmaps",sep="/")))
  #try(file.rename(paste(getwd(),list.files(pattern=".png"),sep="/"),paste(paste(outloc,"KEGGmaps",sep="/"),list.files(pattern=".png"),sep="/")),silent=TRUE)
  
}

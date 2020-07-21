data_preprocess <-
function(Xmat=NA,Ymat=NA,feature_table_file=NA,parentoutput_dir,class_labels_file=NA,num_replicates=3,feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method="mean",
all.missing.thresh=0.5,group.missing.thresh=0.7,
log2transform=TRUE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=TRUE,lowess_norm=FALSE,rangescaling=FALSE,paretoscaling=FALSE,mstus=FALSE,sva_norm=FALSE,TIC_norm=FALSE,eigenms_norm=FALSE,madscaling=FALSE,vsn_norm=FALSE,
cubicspline_norm=FALSE,
missing.val=0,samplermindex=NA, rep.max.missing.thresh=0.5,summary.na.replacement="zeros",featselmethod=NA,pairedanalysis=FALSE,normalization.method="none",input.intensity.scale="raw",create.new.folder=TRUE){
    
    
    options(warn=-1)
    
    #read file; First row is column headers
    if(is.na(Xmat==TRUE)){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
    }else{
        data_matrix<-Xmat
        #rm(Xmat)
    }
    
    setwd(parentoutput_dir)
    
  if(create.new.folder==TRUE){
    dir.create("Stage1")
    setwd("Stage1")
  }else{
    
    dir.create("Tables")
    setwd("Tables")
    
    
  }
   
   normalization.method=tolower(normalization.method)
   
    if(length(normalization.method)==1){
        
        if(normalization.method=="log2quantilenorm" || normalization.method=="log2quantnorm"){
            print("Performing log2 transformation and quantile normalization")
            log2transform=TRUE
            quantile_norm=TRUE
            
        }else{
            if(normalization.method=="log2transform"){
                print("Performing log2 transformation")
                log2transform=TRUE
            }else{
                if(normalization.method=="znormtransform"){
                        print("Performing autoscaling")
                        znormtransform=TRUE
                
                }else{
                    if(normalization.method=="quantile_norm"){
                         print("Performing quantile normalization")
                        quantile_norm=TRUE
                    }else{
                        if(normalization.method=="lowess_norm"){
                            print("Performing Cyclic Lowess normalization")
                            lowess_norm=TRUE
                        }else{
                            
                            if(normalization.method=="rangescaling"){
                                print("Performing Range scaling")
                                rangescaling=TRUE
                            }else{
                                if(normalization.method=="paretoscaling"){
                                    print("Performing Pareto scaling")
                                    paretoscaling=TRUE
                                }else{
                                    
                                    if(normalization.method=="mstus"){
                                        
                                        print("Performing MS Total Useful Signal (MSTUS) normalization")
                                        mstus=TRUE
                                    }else{
                                        
                                        if(normalization.method=="sva_norm"){
                                            
                                            print("Performing Surrogate Variable Analysis (SVA) normalization")
                                            sva_norm=TRUE
                                            log2transform=TRUE
                                        }else{
                                            if(normalization.method=="eigenms_norm"){
                                                print("Performing EigenMS normalization")
                                                eigenms_norm=TRUE
                                                if(input.intensity.scale=="raw"){
                                                    log2transform=TRUE
                                                }
                                                
                                            }else{
                                                if(normalization.method=="vsn_norm"){
                                                    print("Performing variance stabilizing normalization")
                                                    vsn_norm=TRUE
                                                    
                                                }else{
						
								if(normalization.method=="tic_norm"){
								
									print("Performing totial ion intensity normalization")
									
									TIC_norm=TRUE
								}else{
								
									if(normalization.method=="cubicspline_norm"){
									
									
										print("Cubic spline normalization")
										
										cubicspline_norm=TRUE
									}
								
								}
						
						}
						
						
                                            }
                                        }
                                        
                                    }
                                    
                                    
                                }
                                
                            }
                            
                        }
                        
                    }
                    
                    
                }
            }
            
            
        }
        
    
    }
    cnames<-colnames(data_matrix)
    cnames<-tolower(cnames)
    
  
  if(input.intensity.scale=="log2"){
   
   log2transform=FALSE
  }
    check_names<-grep(cnames,pattern="^name$")
  
    names_with_mz_time<-NA
    
    X<-data_matrix
    if(length(check_names)>0){
        
        if(check_names==1){
            
            check_names1<-grep(cnames,pattern="^mz$")
            check_names2<-grep(cnames,pattern="^time$")
            
            
            if(length(check_names1)<1 & length(check_names2)<1){
                mz<-seq(1.00001,nrow(X)+1,1)
                time<-seq(1.01,nrow(X)+1,1.00)
                check_ind<-gregexpr(cnames,pattern="^name$")
                check_ind<-which(check_ind>0)
                X<-as.data.frame(X)
                
                
                Name<-as.character(X[,check_ind])
               
                X<-cbind(mz,time,X[,-check_ind])
                names_with_mz_time=cbind(Name,mz,time)
               
                names_with_mz_time<-as.data.frame(names_with_mz_time)
                X<-as.data.frame(X)
                
               
                write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
                
            }else{
                
                if(length(check_names1)>0 & length(check_names2)>0){
                    
                    check_ind<-gregexpr(cnames,pattern="^name$")
                    check_ind<-which(check_ind>0)
                    Name<-as.character(X[,check_ind])
                    X<-X[,-check_ind]
                    names_with_mz_time=cbind(Name,X$mz,X$time)
                    colnames(names_with_mz_time)<-c("Name","mz","time")
                    names_with_mz_time<-as.data.frame(names_with_mz_time)
                    X<-as.data.frame(X)
                    write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
                }
            }
            
        }
    }else{
        
       
        
        check_names1<-grep(cnames[1],pattern="^mz$")
        check_names2<-grep(cnames[2],pattern="^time$")
        if(length(check_names1)<1 || length(check_names2)<1){
            stop("Invalid feature table format. First two columns should be mz and time. Please check example files.")
        }
    }
    
     data_matrix<-X
    #print(head(data_matrix))
    
    if(is.na(all.missing.thresh)==TRUE){
        
        all.missing.thresh=(-1)
    }
    
    if(is.na(samplermindex)==FALSE){
        data_matrix<-data_matrix[,-c(samplermindex)]
    }
    
    #use only unique records
    data_matrix<-unique(data_matrix)
    
    if(is.na(missing.val)==FALSE){
        
        print("Replacing missing values with NAs.")
        data_matrix<-replace(as.matrix(data_matrix),which(data_matrix==missing.val),NA)
    }
    
    # print(data_matrix[1:10,1:5])

    #print("dim of original data matrix")
    #print(dim(data_matrix))
    data_matrix_orig<-data_matrix
    
    
    snames<-colnames(data_matrix)
    
   
    
    
    dir.create(parentoutput_dir,showWarnings=FALSE)
    #parentoutput_dir<-paste(parentoutput_dir,"/Stage1/",sep="")
    
    #dir.create(parentoutput_dir,showWarnings=FALSE)
    fheader="transformed_log2fc_threshold_"
    #setwd(parentoutput_dir)
    
    data_m<-as.matrix(data_matrix[,-c(1:2)])
    
    if(is.na(Xmat)==FALSE){
        
        #   write.table(Xmat,file="organized_featuretableA.txt",sep="\t",row.names=TRUE)
       
        
    }
    
    if(is.na(Ymat)==FALSE){
        # write.table(Ymat,file="organized_classlabelsA.txt",sep="\t",row.names=FALSE)
        
    }
    
    #Step 2) Average replicates
    if(summarize.replicates==TRUE)
    {
        if(num_replicates>1)
        {
            
            data_m<-getSumreplicates(data_matrix,alignment.tool="apLCMS",numreplicates=num_replicates,numcluster=10,rep.max.missing.thresh=rep.max.missing.thresh,summary.method=summary.method,summary.na.replacement, missing.val=missing.val)
            
            #data_m<-round(data_m,3)
            
            data_m<-replace(data_m,which(is.na(data_m)==TRUE),missing.val)
            
            if(summary.method=="mean"){
                print("Replicate averaging done")
                filename<-paste("Rawdata_averaged.txt",sep="")
            }else{
                if(summary.method=="median"){
                    print("Replicate median summarization done")
                    filename<-paste("Rawdata_median_summarized.txt",sep="")
                }
                
            }
            
            data_m_prenorm<-cbind(data_matrix[,c(1:2)],data_m)
            
            write.table(data_m_prenorm, file=filename,sep="\t",row.names=FALSE)
            
            data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
            #num_samps_group[[1]]=(1/num_replicates)*num_samps_group[[1]]
            #num_samps_group[[2]]=(1/num_replicates)*num_samps_group[[2]]
        }
    }
    
    
    data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
    
    data_matrix_orig<-data_matrix
    data_subjects<-data_m
    
    ordered_labels={}
    
    num_samps_group<-new("list")
    
    
    
    if(is.na(class_labels_file)==FALSE)
    {
        
        print("Class labels file:")
        print(class_labels_file)
        
        data_matrix={}
        
        if(is.na(Ymat)==TRUE){
            classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
        }else{
            classlabels<-Ymat
        }
        
        #class_labels_sampnames<-classlabels[,1]
        #data_matrix_sampnames<-colnames(data_m)
        
        #classlabels<-classlabels[match(class_labels_sampnames,data_matrix_sampnames),]
        
        classlabels<-as.data.frame(classlabels)
        if(pairedanalysis==TRUE){
             classlabels<-classlabels[,-c(2)]
            
        }
       
       
       
        
        cnames1<-colnames(classlabels)
        cnames1[1]<-c("SampleID")
        cnames1[2]<-c("Class")
        
  
        colnames(classlabels)<-cnames1 #c("SampleID","Class")
        
  
        f1<-table(classlabels$SampleID)
        
       
        
        classlabels<-as.data.frame(classlabels)
        classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
        #print(classlabels)
        class_labels_levels<-levels(as.factor(classlabels[,2]))
        bad_rows<-which(class_labels_levels=="")
        if(length(bad_rows)>0){
            class_labels_levels<-class_labels_levels[-bad_rows]
        }
        
        for(c in 1:length(class_labels_levels))
        {
            
            if(c>1){
                data_matrix<-cbind(data_matrix,data_subjects[,which(classlabels[,2]==class_labels_levels[c])])
            }else{
                data_matrix<-data_subjects[,which(classlabels[,2]==class_labels_levels[c])]
            }
            classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
            ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
            num_samps_group[[c]]<-length(classlabels_index)
            
        }
        
        #colnames(data_matrix)<-as.character(ordered_labels)
        data_matrix<-cbind(data_matrix_orig[,c(1:2)],data_matrix)
        data_m<-as.matrix(data_matrix[,-c(1:2)])
        
        
    }else
    {
        if(is.na(Ymat)==TRUE)
        {
            classlabels<-rep("A",dim(data_m)[2])
            classlabels<-as.data.frame(classlabels)
            ordered_labels<-classlabels
            num_samps_group[[1]]<-dim(data_m)[2]
            class_labels_levels<-c("A")
            data_m<-as.matrix(data_matrix[,-c(1:2)])
            
        }else{
            classlabels<-Ymat
            classlabels<-as.data.frame(classlabels)
            
            if(pairedanalysis==TRUE){
                classlabels<-classlabels[,-c(2)]
                
            }
            
            
            
            
            cnames1<-colnames(classlabels)
            cnames1[1]<-c("SampleID")
            cnames1[2]<-c("Class")
            
            colnames(classlabels)<-cnames1
            #colnames(classlabels)<-c("SampleID","Class")
            f1<-table(classlabels$SampleID)
            
            
            classlabels<-as.data.frame(classlabels)
            classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
            #print(classlabels)
            class_labels_levels<-levels(as.factor(classlabels[,2]))
            bad_rows<-which(class_labels_levels=="")
            if(length(bad_rows)>0){
                class_labels_levels<-class_labels_levels[-bad_rows]
            }
            
            for(c in 1:length(class_labels_levels))
            {
                #if(c>1){
                #data_matrix<-cbind(data_matrix,data_subjects[,which(classlabels[,2]==class_labels_levels[c])])
                #}else{
                #	data_matrix<-data_subjects[,which(classlabels[,2]==class_labels_levels[c])]
                #}
                classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
                ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
                num_samps_group[[c]]<-length(classlabels_index)
                
            }
            
            #colnames(data_matrix)<-as.character(ordered_labels)
            #data_matrix<-cbind(data_matrix_orig[,c(1:2)],data_matrix)
            #data_m<-as.matrix(data_matrix[,-c(1:2)])
        }
        
        
        
    }
    
   # #save(class_labels_levels,file="class_labels_levels.Rda")
  #      #save(classlabels,file="classlabels1.Rda")
  #  #save(num_samps_group,file="num_samps_group.Rda")
    ##save(ordered_labels,file="ordered_labels1.Rda")
    
    rnames_xmat<-colnames(data_matrix[,-c(1:2)])
    
    if(is.na(classlabels)==FALSE){
    
    
     rnames_ymat<-as.character(classlabels[,1])
    }else{
        rnames_ymat<-rnames_xmat
        
    }
 
if(FALSE){   
    #for(ind1 in 1:length(rnames_xmat))
    {
    check_ylabel<-regexpr(rnames_ymat[1],pattern="^[0-9]*",perl=TRUE)
    check_xlabel<-regexpr(rnames_xmat[1],pattern="^X[0-9]*",perl=TRUE)
    
    if(length(check_ylabel)>0 && length(check_xlabel)>0){
        if(attr(check_ylabel,"match.length")>0 && attr(check_xlabel,"match.length")>0){
            
            rnames_ymat<-paste("X",rnames_ymat,sep="") #gsub(rnames_ymat,pattern="\\.[0-9]*",replacement="")
            
            
        }
    }
    
    }

    rnames_xmat<-gsub(rnames_xmat,pattern=" |-",replacement=".")
    rnames_ymat<-gsub(rnames_ymat,pattern=" |-",replacement=".")	    

    match_names<-match(rnames_xmat,rnames_ymat)
    
    bad_colnames<-length(which(is.na(match_names)==TRUE))
  

 
    if(length(bad_colnames)>0){

	stop("Sample names do not match between feature table and classlabels. Please try replacing spaces and - with .")
    }else{ 
    classlabels<-classlabels[match_names,]

    }
}

    #Step 3a) Remove features if signal is not detected in at least x% of all samples
    ##################################################################################
    metab_zeros={}
    data_clean<-{}
    clean_metabs<-{}
    total_good_metabs<-{}

    total_sigs<-apply(data_m,1,function(x){
        if(is.na(missing.val)==FALSE){return(length(which(x>missing.val)))
        }else{
            return(length(which(is.na(x)==FALSE)))
        }})
    
    
    useful_metabs_1<-which(total_sigs>=dim(data_m)[2])
    
    useful_metabs_2<-which(total_sigs>=dim(data_m)[2]*0.95)
    
    
    if(is.na(all.missing.thresh)==FALSE)
    {
        
        total_sig_thresh<-dim(data_m)[2]*all.missing.thresh
        
        total_good_metabs<-which(total_sigs>total_sig_thresh)
        
    }
    
    #remove bad features based on all missing values criteria
    if(length(total_good_metabs)>0){
        data_m<-data_m[total_good_metabs,]
        data_matrix<-data_matrix[total_good_metabs,]
        #print(paste("Dimension of data matrix after overall ",all.missing.thresh,"% signal threshold filtering",sep=""))
        print(paste("Dimension of data matrix after using overall ",100*all.missing.thresh, "% signal criteria for filtering:"),sep="")
        print(dim(data_matrix))
    }else{
        stop(paste("None of the metabolites have signal in ",all.missing.thresh*100, "% of samples",sep=""))
    }
    
    
    #Step 3b) Find features for which the signal is not detected in at least x% of samples in either of the groups
    
    
    data_m<-data_matrix[,-c(1:2)]
    


    
    if(is.na(group.missing.thresh)==FALSE)
    {
        
    
                         if(length(class_labels_levels)>1){
                    

                   

                           clean_metabs<-lapply(1:dim(data_matrix)[1],function(metab_num)
                           {
                             clean_metabs<-NA
                             for(c in 1:length(class_labels_levels)){

		                  		classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
			                  	templabels<-classlabels[,2]
                            if(is.na(missing.val)==FALSE){
                                num_cursig<-length(which(data_m[metab_num,which(templabels==class_labels_levels[c])]>missing.val))
                            }else{
                                num_cursig<-length(which(is.na(data_m[metab_num,which(templabels==class_labels_levels[c])])==FALSE))
                            }
                            sig_thresh_cur<-length(which(templabels==class_labels_levels[c]))*group.missing.thresh
                            if(num_cursig>=sig_thresh_cur)
                            {
                                clean_metabs<-metab_num
                                break   #for(i in 1:4){if(i==3){break}else{print(i)}}
                                
                            }
                            
                             }
                             return(clean_metabs)
                        })
                }
                else{
                    
                    
                    
                    if(length(class_labels_levels)==1){
                        num_samps_group[[1]]<-num_samps_group[[1]]
                        
                        
                        sig_thresh_groupA<-group.missing.thresh*num_samps_group[[1]]
                        
                        
                        clean_metabs<-lapply(1:dim(data_matrix)[1],function(metab_num)
                        {
                            if(is.na(missing.val)==FALSE){
                                num_sigsA<-length(which(data_m[metab_num,1:num_samps_group[[1]]]>missing.val))
                                
                            }else{
                                
                                num_sigsA<-length(which(is.na(data_m[metab_num,1:num_samps_group[[1]]])==FALSE))
                            }
                            
                            if((num_sigsA>=sig_thresh_groupA) )
                            {
                                #clean_metabs<-c(clean_metabs,metab_num)
                              
                                return(metab_num)
                            }else{
                              
                               return(NA)
                            }
                            
                        })
                    }
                
            
                }
      
      clean_metabs<-unlist(clean_metabs)
      clean_metabs<-na.omit(clean_metabs)
        
        
        
    }else{
    

      
      clean_metabs<-seq(1,dim(data_matrix)[1])
    }
    ####################################################################################
    
    #Step 4) Replace missing values
    if(summarize.replicates==TRUE)
    {
        
        {
            
            if(is.na(missing.val)==FALSE){
                
                print("Replacing missing values with NAs.")
                data_m<-replace(as.matrix(data_m),which(data_m==missing.val),NA)
            }

            
            if(summary.na.replacement=="zeros"){
                data_m<-replace(data_m,which(is.na(data_m)==TRUE),0)
            }else{
                if(summary.na.replacement=="halfsamplemin"){
                    data_m<-apply(data_m,2,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
                }else{
                    
                    if(summary.na.replacement=="halfdatamin"){
                        
                        
                        min_val<-min(data_m,na.rm=TRUE)*0.5
                        data_m<-replace(data_m,which(is.na(data_m)==TRUE),min_val)
                        
                        #data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
                    }else{
                        if(summary.na.replacement=="halffeaturemin"){
                            data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
                            data_m<-t(data_m)
                        }else{
                            
                            
                            if(summary.na.replacement=="bpca"){
                                library(pcaMethods)
                                
                              
                                pc1 <- pcaMethods::pca(t(data_m), method="bpca", nPcs=3,scale="uv")
                               
                                data_m<-pcaMethods::completeObs(pc1)
                              
                                try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)
                                
                                data_m<-t(data_m)
                                
                            }else{

                                 if(summary.na.replacement=="knn"){
                                   suppressMessages(library(impute))
				
                                    data_m<-apply(data_m,2,as.numeric)
                                    data_m<-impute.knn(data.matrix(data_m),k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
                                    data_m<-data_m$data

                                 }else{

                                    if(summary.na.replacement=="featuremean"){
                            			data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(mean(x,na.rm=TRUE))}; return(x)})
                            			data_m<-t(data_m)
                                    }else{
                                        
                                        if(summary.na.replacement=="randomforest"){
                                            
                                                final_set<-RFimpute(t(data_m))
                                        
                                        }else{
                                            
                                            if(summary.na.replacement=="QRILC"){
                                                                   
                                                                       final_set<-QRILCimpute(data_m)
                                                               
                                            }else{
                                                data_m<-data_m
                                                
                                            }
                                           
                                            
                                        }
                                        
                                        
                                    }
                                    
                                }
                            }
                            
                            
                            
                        }
                    }
                }
                
                
            }
        }
    }else
    {
        data_m<-data_matrix[,-c(1:2)]
        
        if(is.na(missing.val)==FALSE){
            
            print("Replacing missing values with NAs.")
            data_m<-replace(as.matrix(data_m),which(data_m==missing.val),NA)
        }
        
        if(summary.na.replacement=="zeros"){
            data_m<-replace(data_m,which(is.na(data_m)==TRUE),0)
        }else{
            if(summary.na.replacement=="halfsamplemin"){
                data_m<-apply(data_m,2,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
            }else{
                
                if(summary.na.replacement=="halfdatamin"){
                    
                    
                    min_val<-min(data_m,na.rm=TRUE)*0.5
                    data_m<-replace(data_m,which(is.na(data_m)==TRUE),min_val)
                    
                    #data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
                }else{
                    if(summary.na.replacement=="halffeaturemin"){
                        data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
                        data_m<-t(data_m)
                    }else{
                        
                        if(summary.na.replacement=="bpca"){
                            
                            suppressMessages(library(pcaMethods))
                            pc1 <- pcaMethods::pca(t(data_m), method="bpca", nPcs=3, scale="uv")
                            
                            data_m<-completeObs(pc1)
                            data_m<-t(data_m)
                        }else{
                            if(summary.na.replacement=="knn"){
                                suppressMessages(library(impute))
                                data_m<-impute.knn(data.matrix(data_m),k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
                                data_m<-data_m$data
                            }else{
                                
                                
                                if(summary.na.replacement=="randomforest"){
                                                       
                                                           final_set<-RFimpute(t(data_m))
                                                   
                                                   }else{
                                                       
                                                       if(summary.na.replacement=="QRILC"){
                                                                              
                                                                                  final_set<-QRILCimpute(data_m)
                                                                          
                                                                          }
                                                      
                                                       
                                                   }
                                
                                
                                
                            }
                        }
                        
                    }
                }
            }
            
            
        }
        
    }
    
    
   
    
    #group-wise missing values
    if(length(clean_metabs)>0)
    {
        data_m<-data_m[clean_metabs,]
        data_matrix<-data_matrix[clean_metabs,]
        
        print(paste("Dimension of data matrix after using group-wise (Factor 1) ",100*group.missing.thresh, "% signal criteria for filtering:"),sep="")
        print(dim(data_matrix))
        
    }
    
    data_m<-as.data.frame(data_m)
    data_matrix<-as.data.frame(data_matrix)
    
    #save(data_matrix,file="data_matrix.Rda")
    #save(data_m,file="data_m.Rda")
    
    
     data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
    write.table(data_matrix,file="pretransformation.txt",sep="\t",row.names=FALSE)
    ####################################################################
    #Step 4) Data transformation and normalization
    
    data_m_prescaling<-data_m
    
    
    

    
    if(TIC_norm==TRUE){
        
        ####savedata_m,file="data_m_raw.Rda")
       
        
        data_m<-do_TIC_norm(data_m)
        
    }else{
        
        if(mstus==TRUE){
            
           
           data_m<-do_MSTUS_norm(data_m,missing.val=missing.val)
            
            
            
        }
        
    }
    
    if(cubicspline_norm==TRUE){
    
        data_m<-do_cubicspline_norm(data_m)
    
    }
      if(log2transform==TRUE)
    {
        data_m<-do_log2transform_norm(data_m)
        
       
    }
  #  #save(data_m,classlabels,data_matrix,file="debugnorm.Rda")
    
    if(eigenms_norm==TRUE){
        
        feature_id_vector<-paste(data_matrix[,c(1)],"_",data_matrix[,c(2)],sep="")
        
       
        
        if(featselmethod=="limma2wayrepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="spls2wayrepeat"){
            analysistype="twowayrepeat"
            pairedanalysis=TRUE
        }else{
            
            if(featselmethod=="limma1wayrepeat" | featselmethod=="lm1wayanovarepeat" | featselmethod=="spls1wayrepeat"){
                                   analysistype="onewayrepeat"
                                   pairedanalysis=TRUE
            }else{
                pairedanalysis=FALSE
                
            }
            
        }
        
        data_m<-do_eigenms_norm(data_m,classlabels,featselmethod,feature_id_vector,pairedanalysis=pairedanalysis)
       
    }
    
    if(sva_norm==TRUE){
        
       data_m<-do_sva_norm(data_m,classlabels,featselmethod)
        
    }
   
    
    
	
   if(quantile_norm==TRUE)
    {
        data_m<-do_quantile_norm(data_m)
        
        
    }
    
    if(vsn_norm==TRUE)
    {
        data_m<-do_vsn_norm(data_m)
       
        
        
    }
    
    if(lowess_norm==TRUE)
    {
        data_m<-do_loess_norm(data_m)
        #print("lowess")
    }
    
  
    
    if(medcenter==TRUE)
    {
       data_m<-do_medcenter_norm(data_m)
        
        
    }
    if(znormtransform==TRUE)
    {
        data_m<-do_znormtransform_norm(data_m)
    }
    
    
   
    if(paretoscaling==TRUE){
        #pareto scaling
        data_m<-do_paretoscaling_norm(data_m)
        
    }
    
    if(rangescaling==TRUE){
        #range scaling
        data_m<-do_rangescaling_norm(data_m)
    }
    
    
   

    
    data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
    
   
    data_m<-round(data_m,5)
    data_m<-as.data.frame(data_m)
    
    
    num_rows<-dim(data_m)[1]
    num_columns<-dim(data_m)[2]
    
    #print("num rows is ")
    #print(num_rows)
    #for apLCMS:
    rnames<-paste("mzid_",seq(1,num_rows),sep="")
    rownames(data_m)=rnames
    
    rnames_xmat<-colnames(data_m)
    rnames_ymat<-classlabels[,1]
    
    check_ylabel<-regexpr(rnames_ymat[1],pattern="^[0-9]*",perl=TRUE)
     check_xlabel<-regexpr(rnames_xmat[1],pattern="^X[0-9]*",perl=TRUE)
    
    
     if(length(check_ylabel)>0 && length(check_xlabel)>0){
         if(attr(check_ylabel,"match.length")>0 && attr(check_xlabel,"match.length")>0){
             
             rnames_ymat<-paste("X",rnames_ymat,sep="")
         }
     }
    classlabels<-classlabels[match(rnames_xmat,rnames_ymat),]
    
    
    filename<-paste("ordered_classlabels_file.txt",sep="")
    write.table(classlabels, file=filename,sep="\t",row.names=FALSE)
    
    filename<-paste("Normalized_sigthreshfilt_averaged_data.txt",sep="")
    data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
    write.table(data_matrix, file=filename,sep="\t",row.names=FALSE)
    data_matrix_prescaling<-cbind(data_matrix[,c(1:2)],data_m_prescaling)
    setwd("../")
    return(list(data_matrix_afternorm_scaling=data_matrix,data_matrix_prescaling=data_matrix_prescaling,classlabels=classlabels,names_with_mz_time=names_with_mz_time))
    #return(data_matrix)
}

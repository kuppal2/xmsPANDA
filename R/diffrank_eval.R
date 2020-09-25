diffrank_eval <-
function(feature_table_file=NA,class_labels_file=NA,X=NA,Y=NA,sigfeats=NA,sigfeatsind=NA,cor.method="pearson",num_nodes=2,abs.cor.thresh=0.4,pvalue.thresh=0.005,cor.fdrthresh=0.2,fdrmethod="Strimmer",degree.centrality.method="eigenvector"){
    
    suppressMessages(library(WGCNA))
    
    #print("degree eval")
    
    if(is.na(X[1])==TRUE){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
        
    }else{
        
        data_matrix<-X
    }
    if(is.na(Y[1])==TRUE){
        classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
        
    }else{
        
        classlabels<-Y
    }
    
    
    classlabels<-as.data.frame(classlabels)
    #print(dim(classlabels))
    # print(length(classlabels))
    
    
    if(dim(classlabels)[2]>2){
        classgroup<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
    }else{
        
        classgroup<-classlabels[,2]
    }
    classlabels<-as.data.frame(classlabels)
    
    class_labels_levels<-levels(as.factor(classgroup))
    
    rnames<-paste(sprintf("%.4f",data_matrix$mz),data_matrix$time,sep="_")
    
    data_matrix_orig<-data_matrix
    data_matrix<-data_matrix[,-c(1:2)]
    #data_matrix<-na.omit(data_matrix)
    
    rnamesAB<-gsub(pattern="NA_NA",replacement=NA,x=rnames)
    rnamesAB<-na.omit(rnamesAB)
    
    nSets = length(class_labels_levels);
    multiExpr = vector(mode = "list", length = nSets)
    data_matrix_list<-new("list")
    num_samps_groups<-new("list")
    degree_list<-new("list")
    
    data_matrix_all<-t(data_matrix)
    
    
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = pickSoftThreshold(data=data_matrix_all, dataIsExpr=TRUE,powerVector = powers, verbose = 0)
    power_val=sft$powerEstimate
    
    if(is.na(power_val)==TRUE){
        power_val=6
    }
    
 
    
    #degree_overall<-softConnectivity(datExpr=data_matrix_all,power=power_val,minNSamples=10)
    #degree_overall<-replace(degree_overall,which(degree_overall<0),0)
    
    
    sig_status<-rep(0,dim(data_matrix_orig)[1])
    
    colors_met<-cbind(data_matrix_orig[,c(1:2)],sig_status)
    
    colors_met<-as.data.frame(colors_met)
    
    colors_met$sig_status[sigfeatsind]=1
    
    colnames(colors_met)<-c("mz","time","sig_status")
    
    colnames_vec<-c("mz","time","sig_status")
    
    sig_ind<-which(colors_met$sig_status==1)
    
    
    names=paste(round(colors_met$mz,4),round(colors_met$time,0),sep="_")
    
    #names=paste(round(colors_met$mz,4),sep="_")
    

     cl<-makeCluster(max(2,detectCores()*0.5))
     clusterEvalQ(cl,library(WGCNA))
     #save(data_matrix,abs.cor.thresh,cor.method,pvalue.thresh,power_val,class_labels_levels,classgroup,file="diffrankdebug.Rda")
    adj_mat_list<-parLapply(cl,1:length(class_labels_levels),function(i,data_matrix,abs.cor.thresh,cor.method,pvalue.thresh,power_val,class_labels_levels,classgroup){
        
        data_matrix_list[[i]]<-t(data_matrix[,which(classgroup==class_labels_levels[i])])
        num_samps_groups[[i]]<-dim(data_matrix_list[[i]])[1]
        #print(dim(data_matrix_list[[i]]))
        multiExpr[[i]]<-list(data = as.data.frame(data_matrix_list[[i]]));
        rownames(multiExpr[[i]]$data)=c(paste(rep(class_labels_levels[i],num_samps_groups[[i]]),seq(1,num_samps_groups[[i]]),sep=""))
        
        cor.method="pearson"
        simmat<-WGCNA::cor(data_matrix_list[[i]],nThreads=num_nodes,method=cor.method,use = 'p')
        
        simmat<-round(simmat,3)
        pearson_resvec<-as.vector(simmat)
        
        #complete_pearsonqvalue_mat<-unlist(fdr_adjust_pvalue)
        
        cor_vec=seq(0,1,0.01)
        pvalues_vec<-corPvalueStudent(cor=cor_vec,nSamples=num_samps_groups[[i]])
        abs.cor.thresh=max(abs.cor.thresh,min(cor_vec[which(pvalues_vec<pvalue.thresh)],na.rm=TRUE),na.rm=TRUE)
        
        #dim(complete_pearsonqvalue_mat)<-dim(simmat)
        
        simmat[which(abs(simmat)<abs.cor.thresh)]<-0
        #simmat[which(complete_pearsonqvalue_mat>fdrthresh)]<-0
      #  save(simmat,sigfeatsind,rnames,rnamesAB,names,file="c1.Rda")
        rownames(simmat)<-rnamesAB
        colnames(simmat)<-rnamesAB
        simmat_sig<-simmat[sigfeatsind,sigfeatsind]
        
        fname=paste("Tables/cor.matrix.selectedfeatures.class",class_labels_levels[i],".corthresh",abs.cor.thresh,".txt",sep="")
        write.table(simmat_sig,file=fname,sep="\t")
       
        ADJdataOne<-adjacency.fromSimilarity(similarity=simmat,power=power_val)
        
        #degree_list[[i]]<-softConnectivity(datExpr=data_matrix_list[[i]],power=power_val,minNSamples=2)
        #
    
        return(ADJdataOne)
    },data_matrix=data_matrix,abs.cor.thresh=abs.cor.thresh,cor.method=cor.method,pvalue.thresh=pvalue.thresh,power_val=power_val,class_labels_levels=class_labels_levels,classgroup=classgroup)
stopCluster(cl)
    #the first class is used as reference
    diffrank_list<-lapply(2:length(class_labels_levels),function(i){
        
        res<-diffRank(adj_mat_list[[1]],adj_mat_list[[i]],degree.centrality.method=degree.centrality.method)
        return(res)
    })
    
    #savediffrank_list,file="diffrank_list.Rda")
    
    diffrank_mat<-do.call("cbind",diffrank_list)
    
    
    
    if(length(class_labels_levels)>2){
        
        diffrank_res<-apply(diffrank_mat,1,function(x){max(x,na.rm=TRUE)})
        
        diffrank_res_rank<-rank(-1*diffrank_res)
        colnames_vec<-c("mz","time","sigstatus","Max.Delta.Centrality","DiffRank",paste("DiffRank.",class_labels_levels[1],"vs",class_labels_levels[-1],sep=""))
        colors_met_all<-cbind(colors_met,diffrank_res,diffrank_res_rank,diffrank_mat)
        
        
    }else{
        diffrank_res<-diffrank_mat
        
        diffrank_res_rank<-rank(-1*diffrank_res)
        
        colnames_vec<-c("mz","time","sigstatus","Delta.Centrality","DiffRank")
        colors_met_all<-cbind(colors_met,diffrank_res,diffrank_res_rank)
        
        
    }
    
    #  print(head(colors_met_all))
   
    colors_met_all<-as.data.frame(colors_met_all)
    colnames(colors_met_all)<-colnames_vec
    
    #  print(head(colors_met_all))
    #print(table(MET[[1]]$validColors))
    #print(table(classAmoduleColors))
    write.table(colors_met_all,file="Tables/DiNA_eval_allfeats.txt",sep="\t",row.names=FALSE)
    
    sub_colors_met<-{}
    if(is.na(sigfeats)==FALSE){
        sub_colors_met<-colors_met_all[sigfeatsind,]
        write.table(sub_colors_met,file="Tables/DiNA_eval_selectfeats.txt",sep="\t",row.names=FALSE)
    }
    
    ####saveMET,file="MET.Rda")
    
    
    
    return(list(all=colors_met_all,sigfeats=sub_colors_met))
}

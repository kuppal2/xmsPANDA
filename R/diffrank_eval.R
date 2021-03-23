diffrank_eval <-
function(feature_table_file=NA,class_labels_file=NA,X=NA,Y=NA,sigfeats=NA,
                        sigfeatsind=NA,cor.method="pearson",num_nodes=2,abs.cor.thresh=0.4,pvalue.thresh=0.005,
                        cor.fdrthresh=0.2,fdrmethod="Strimmer",degree.centrality.method="eigenvector",alphabetical.order=FALSE,
                        node_names=NA,plot_graph_bool=FALSE,networktype="complete"){
  
  suppressMessages(library(WGCNA))
  
  #print("degree eval")
  
  if(typeof(X)=="logical"){
    data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
    
  }else{
    
    data_matrix<-X
  }
  if(typeof(Y)=="logical"){
    classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
    
  }else{
    
    classlabels<-Y
  }
  
  
  classlabels<-as.data.frame(classlabels)
  #print(dim(classlabels))
  # print(length(classlabels))
  
  
  if(dim(classlabels)[2]>2){
    
    if(alphabetical.order==FALSE){
      
      classlabels[,2]<-factor(classlabels[,2],levels=unique(classlabels[,2]))
      classlabels[,3]<-factor(classlabels[,3],levels=unique(classlabels[,3]))
      
    }
    
    classgroup<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
  }else{
    
    if(alphabetical.order==FALSE){
      
      classlabels[,2]<-factor(classlabels[,2],levels=unique(classlabels[,2]))
     
    }
    
    classgroup<-classlabels[,2]
  }
  
  classlabels<-as.data.frame(classlabels)
  
  if(alphabetical.order==FALSE){
    
    classgroup<-factor(classgroup,levels=unique(classgroup))
    class_labels_levels<-levels(as.factor(classgroup))
    
  }else{
    class_labels_levels<-levels(as.factor(classgroup))
  }
#  rnames<-paste(sprintf("%.4f",data_matrix$mz),data_matrix$time,sep="_")
  
  rnames<-paste(data_matrix$mz,data_matrix$time,sep="_")
  
  data_matrix_orig<-data_matrix
  
  #addition start
  cnames<-colnames(data_matrix)
  cnames<-tolower(cnames)
  
  check_names<-grep(cnames,pattern="^name$")
  
  
  if(length(check_names)>0){
    
    if(check_names==1){
      
      check_names1<-grep(cnames,pattern="^mz$")
      check_names2<-grep(cnames,pattern="^time$")
      
      
      if(length(check_names1)<1 & length(check_names2)<1){
        mz<-seq(1.00001,nrow(data_matrix)+1,1)
        time<-seq(1.01,nrow(data_matrix)+1,1.00)
        check_ind<-gregexpr(cnames,pattern="^name$")
        check_ind<-which(check_ind>0)
        data_matrix<-as.data.frame(data_matrix)
        
        Name<-as.character(data_matrix[,check_ind])
        name=Name
        data_matrix<-cbind(mz,time,data_matrix[,-check_ind])
        names_with_mz_time=cbind(Name,mz,time)
        
        names_with_mz_time<-as.data.frame(names_with_mz_time)
        data_matrix<-as.data.frame(data_matrix)
        
        write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        
        
      }else{
        
        if(length(check_names1)>0 & length(check_names2)>0){
          
          check_ind<-gregexpr(cnames,pattern="^name$")
          check_ind<-which(check_ind>0)
          Name<-as.character(data_matrix[,check_ind])
          name=Name
          data_matrix<-data_matrix[,-check_ind]
          names_with_mz_time=cbind(Name,data_matrix$mz,data_matrix$time)
          colnames(names_with_mz_time)<-c("Name","mz","time")
          names_with_mz_time<-as.data.frame(names_with_mz_time)
          data_matrix<-as.data.frame(data_matrix)
          write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        }
      }
      
    }
  }else{
    
    
    check_names1<-grep(cnames[1],pattern="^mz$")
    check_names2<-grep(cnames[2],pattern="^time$")
    if(length(check_names1)<1 || length(check_names2)<1){
      stop("Invalid feature table format. The format should be either Name in column A or mz and time in columns A and B. Please check example files.")
    }
    
    data_matrix<-data_matrix[,-c(1:2)]
    
    #names_with_mz_time=cbind(paste(data_matrix$mz,"_",data_matrix$time,sep=""),data_matrix$mz,data_matrix$time)
    #colnames(names_with_mz_time)<-c("Name","mz","time")
    #names_with_mz_time<-as.data.frame(names_with_mz_time)
    
  }
  
  
  

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
  clusterEvalQ(cl,library(parallel))
  clusterEvalQ(cl,library(corpcor))
  # save(data_matrix,abs.cor.thresh,cor.method,pvalue.thresh,power_val,class_labels_levels,node_names,classgroup,file="diffrankdebug.Rda")
  
  adj_mat_list<-parLapply(cl,1:length(class_labels_levels),function(i,data_matrix,abs.cor.thresh,cor.method,pvalue.thresh,power_val,class_labels_levels,classgroup){
  # adj_mat_list<-lapply(1:length(class_labels_levels),function(i,data_matrix,abs.cor.thresh,cor.method,pvalue.thresh,power_val,class_labels_levels,classgroup){
    
     print(i)
    data_matrix_list<-t(data_matrix[,which(classgroup==class_labels_levels[i])])
    num_samps_groups<-dim(data_matrix_list)[1]
    
    num_nodes=max(2,detectCores()*0.5)
    
    if(networktype=="complete"){
    #cor.method="pearson"
    simmat<-suppressWarnings(WGCNA::cor(data_matrix_list,nThreads=num_nodes,method=cor.method,use = 'p'))
    }else{
      
      simmat<-suppressWarnings(pcor.shrink(data_matrix_list))
      
      simmat<-replace(simmat,which(simmat>1),1)
      simmat<-replace(simmat,which(simmat<(-1)),-1)
      
    }
    if(length(which(is.na(simmat)==TRUE))>0){
      
      simmat<-replace(simmat,which(is.na(simmat)==TRUE),0)
      #ADJdataOne<-NA
      
      cor_vec=seq(0,1,0.01)
      #pvalues_vec<-corPvalueStudent(cor=cor_vec,nSamples=num_samps_groups[i])
     # abs.cor.thresh=max(abs.cor.thresh,min(cor_vec[which(pvalues_vec<pvalue.thresh)],na.rm=TRUE),na.rm=TRUE)
      
      simmat[which(abs(simmat)<abs.cor.thresh)]<-0
      
      rownames(simmat)<-rnamesAB
      colnames(simmat)<-rnamesAB
      
      
      #".corthresh",abs.cor.thresh,
      fname=paste("Tables/cor.matrix.allfeatures.class",class_labels_levels[i],".txt",sep="")
      #write.table(simmat,file=fname,sep="\t")
      
      simmat_sig<-simmat[sigfeatsind,sigfeatsind]
      
      fname=paste("Tables/cor.matrix.selectedfeatures.class",class_labels_levels[i],".txt",sep="")
      write.table(simmat_sig,file=fname,sep="\t")
      
      power_val=1
      #adjacency matrix for all features
      ADJdataOne<-simmat
    }else{
    simmat<-round(simmat,3)
    
    #simmat<-replace(simmat,which(is.na(simmat)==TRUE),0)
    pearson_resvec<-as.vector(simmat)
    
    
    cor_vec=seq(0,1,0.01)
   # pvalues_vec<-corPvalueStudent(cor=cor_vec,nSamples=num_samps_groups[i])
    #abs.cor.thresh=max(abs.cor.thresh,min(cor_vec[which(pvalues_vec<pvalue.thresh)],na.rm=TRUE),na.rm=TRUE)
    
    simmat[which(abs(simmat)<abs.cor.thresh)]<-0
   
    rownames(simmat)<-rnamesAB
    colnames(simmat)<-rnamesAB
    
    
    
    fname=paste("Tables/cor.matrix.allfeatures.class",class_labels_levels[i],".corthresh",abs.cor.thresh,".txt",sep="")
    #write.table(simmat,file=fname,sep="\t")
    
    simmat_sig<-simmat[sigfeatsind,sigfeatsind]
    
    fname=paste("Tables/cor.matrix.selectedfeatures.class",class_labels_levels[i],".corthresh",abs.cor.thresh,".txt",sep="")
    write.table(simmat_sig,file=fname,sep="\t")
    
    power_val=1
    #adjacency matrix for all features
    ADJdataOne<-simmat #adjacency.fromSimilarity(similarity=simmat,power=power_val,type="signed")
    }
   
    return(ADJdataOne)
  },data_matrix=data_matrix,abs.cor.thresh=abs.cor.thresh,cor.method=cor.method,pvalue.thresh=pvalue.thresh,power_val=power_val,
  class_labels_levels=class_labels_levels,classgroup=classgroup)
   
  stopCluster(cl)
  
  
  if(plot_graph_bool==TRUE){
    
    pdf("Figures/DiNa_graphs.pdf")
    
  }
  
  
  count_connections_perclass<-lapply(1:length(class_labels_levels),function(i){
    print(i)
    
    diag(adj_mat_list[[i]])<-0
    
    length(which(abs(adj_mat_list[[i]])>abs.cor.thresh))
                  
  })
  
  count_connections_perclass<-unlist(count_connections_perclass)
  
 # save(count_connections_perclass,adj_mat_list,class_labels_levels,file="debugnames.Rda")
  
  names(count_connections_perclass)<-class_labels_levels #rownames(adj_mat_list[[1]])
  
  if(plot_graph_bool==TRUE){
    
    barplot(count_connections_perclass,main=paste("total # of connections per class\n |r| threshold: ", abs.cor.thresh,sep=""),
            horiz=TRUE,las=2,cex.text=0.7,col="#69b3a2",xlab="Number of connections",cex.lab=0.7,cex.axis=0.7,cex.names=0.5)
    par(mfrow=c(2,2))
  }
  
  if(min(count_connections_perclass,na.rm=TRUE)>0){
    
    ignore.edge.weights=FALSE
  }else{
    
    ignore.edge.weights=TRUE
  }
  
  ignore.edge.weights=FALSE
  
  
  #the first class is used as reference
  diffrank_list<-lapply(2:length(class_labels_levels),function(i){
   # print(i)
   # if(typeof(adj_mat_list[[i]])!="logical"){
    #save(i,adj_mat_list,degree.centrality.method,class_labels_levels,file="m1.Rda")
    
    if(count_connections_perclass[i]>0){
      
      #Returns delta centrality between groups
    res<-diffRank(adj_mat_list[[1]],adj_mat_list[[i]],degree.centrality.method=degree.centrality.method,
                  class_labels=c(class_labels_levels[1],class_labels_levels[i]),node_names=node_names,plot_graph_bool=plot_graph_bool,ignore.edge.weights=ignore.edge.weights)
    
    
    
    return(res)
    }else{
      
      res<-rep(NA,nrow(adj_mat_list[[1]]))
    }
    
  })
  if(plot_graph_bool==TRUE){
    
  dev.off()
    
  }
#  save(diffrank_list,class_labels_levels,file="diffrank_list.Rda")
  
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
  
  
  
  colors_met_all<-as.data.frame(colors_met_all)
  colnames(colors_met_all)<-colnames_vec
  
  
  #write.table(colors_met_all,file="Tables/DiNA_eval_allfeats.txt",sep="\t",row.names=FALSE)
  
  sub_colors_met<-{}
  if(is.na(sigfeats)==FALSE){
    sub_colors_met<-colors_met_all[sigfeatsind,]
    #write.table(sub_colors_met,file="Tables/DiNA_eval_selectfeats.txt",sep="\t",row.names=FALSE)
  }
  
  ####saveMET,file="MET.Rda")
  
  colnames(diffrank_mat)<-paste("DiffRank.",class_labels_levels[1],"vs",class_labels_levels[-1],sep="")
  
  rownames(diffrank_mat)<-rownames(adj_mat_list[[1]])
  write.table(diffrank_mat,file="Tables/DiNA_pairwisedeltacentrality.txt",sep="\t",row.names=TRUE)
  
  if(plot_graph_bool==TRUE){
    pdf("Figures/DiNA_delta_centrality.pdf") 
    
    for(i in 1:ncol(diffrank_mat)){
      
      try(barplot(diffrank_mat[,i],main=paste("", colnames(diffrank_mat)[i],sep=""),
              horiz=TRUE,las=2,cex.text=0.7,col="#69b3a2",xlab=paste("Delta centrality (",degree.centrality.method,")",sep=""),
              cex.lab=0.7,cex.axis=0.7,cex.names=0.5,xlim=c(0,1)),silent=TRUE)
    }
    
    dev.off()
  }
  
  return(list(all=colors_met_all,sigfeats=sub_colors_met))
}

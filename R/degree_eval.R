degree_eval <-
function(feature_table_file=NA,class_labels_file=NA,X=NA,Y=NA,sigfeats=NA,sigfeatsind=NA){
  
  
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
    classgroup<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
  }else{
    
    classgroup<-classlabels[,2]
  }
  classlabels<-as.data.frame(classlabels)
  
  class_labels_levels<-levels(as.factor(classgroup))
  
  rnames<-paste(data_matrix$mz,data_matrix$time,sep="_") #sprintf("%.4f",
  
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
  
  
  degree_overall<-softConnectivity(datExpr=data_matrix_all,power=power_val,minNSamples=10)
  degree_overall<-replace(degree_overall,which(degree_overall<0),0)
  
  
  sig_status<-rep(0,dim(data_matrix_orig)[1])
  
  colors_met<-cbind(data_matrix_orig[,c(1:2)],degree_overall,sig_status)
  
  colors_met<-as.data.frame(colors_met)
  
  colors_met$sig_status[sigfeatsind]=1
  
  colnames(colors_met)<-c("mz","time","degree_overall","sig_status")
  
  colnames_vec<-c("mz","time","degree_overall","sig_status")
  
  # print(sigfeats$mz)
  
  
  colors_met_all<-colors_met
  
  colors_met<-colors_met[order(colors_met$degree_overall,decreasing=TRUE),]
  
  sig_ind<-which(colors_met$sig_status==1)
  
  # print(summary(colors_met$degree_overall))
  
  #pdf("DICE_plots.pdf")
  
  names=paste(round(colors_met$mz,3),round(colors_met$time,0),sep="_")
  
  names=paste(round(colors_met$mz,4),sep="_")
  
  # print(names[sig_ind])
  
  if(FALSE)
  {
    plot(colors_met$degree_overall,cex=0.2,main="Overall degree distribution",ylab="Degree",xlab="Feature Index",col="orange",type="b",lwd=0.5)
    
    #for(i in 1:length(sig_ind)){
    lapply(1:length(sig_ind),function(i){
      points(y=colors_met$degree_overall[sig_ind[i]],x=sig_ind[i],col="darkgreen",cex=0.8,lwd=2)
      if(i%%2>0){
        text(y=colors_met$degree_overall[sig_ind[i]],x=sig_ind[i],names[sig_ind[i]],cex=0.31,adj=c(1,2))
      }else{
        text(y=colors_met$degree_overall[sig_ind[i]],x=sig_ind[i],names[sig_ind[i]],cex=0.31,adj=c(0,-1))
      }
      
    })
  }
  
  
  for(i in 1:length(class_labels_levels)){
    
    data_matrix_list[[i]]<-t(data_matrix[,which(classgroup==class_labels_levels[i])])
    num_samps_groups[[i]]<-dim(data_matrix_list[[i]])[1]
    #print(dim(data_matrix_list[[i]]))
    multiExpr[[i]]<-list(data = as.data.frame(data_matrix_list[[i]]));
    rownames(multiExpr[[i]]$data)=c(paste(rep(class_labels_levels[i],num_samps_groups[[i]]),seq(1,num_samps_groups[[i]]),sep=""))
    
    degree_list[[i]]<-softConnectivity(datExpr=data_matrix_list[[i]],power=power_val,minNSamples=2)
    #
    colnames_vec<-c(colnames_vec,paste("DegreeClass",i,sep=""))
    
    degree_list[[i]]<-replace(degree_list[[i]],which(degree_list[[i]]<0),0)
    
    #degree_list[[i]][which(is.na(degree_list[[i]])==TRUE)]=1
    
    colors_met<-cbind(data_matrix_orig[,c(1:2)],degree_list[[i]],sig_status)
    
    colors_met<-as.data.frame(colors_met)
    
    colors_met_all<-cbind(colors_met_all,degree_list[[i]])
    
    colors_met$sig_status[sigfeatsind]=1
    
    colnames(colors_met)<-c("mz","time","DegreeClass","sig_status")
    
    #colnames_vec<-c("mz","time","DegreeClass","sig_status")
    colors_met<-as.data.frame(colors_met)
    
    colors_met<-colors_met[order(colors_met$DegreeClass,decreasing=TRUE),]
    
    sig_ind<-which(colors_met$sig_status==1)
    
    names=paste(round(colors_met$mz,4),sep="_")
    
    # print(names[sig_ind])
    mainlab=paste("Class ",i," degree distribution",sep="")
    
    if(FALSE)
    {
      plot(colors_met$DegreeClass,cex=0.2,main=mainlab,ylab="Degree",xlab="Feature Index",col="orange",type="b",lwd=0.5)
      #for(i in 1:length(sig_ind))
      lapply(1:length(sig_ind),function(i)
      {
        points(y=colors_met$DegreeClass[sig_ind[i]],x=sig_ind[i],col="darkgreen",cex=0.8,lwd=2)
        
        if(i%%2>0){
          text(y=colors_met$DegreeClass[sig_ind[i]],x=sig_ind[i],names[sig_ind[i]],cex=0.31,adj=c(1,2))
        }else{
          text(y=colors_met$DegreeClass[sig_ind[i]],x=sig_ind[i],names[sig_ind[i]],cex=0.31,adj=c(0,-1))
        }
      })
      
    }
    
  }
  
  # dev.off()
  
  colors_met_all<-as.data.frame(colors_met_all)
  colnames(colors_met_all)<-colnames_vec
  
  
  #print(table(MET[[1]]$validColors))
  #print(table(classAmoduleColors))
  write.table(colors_met_all,file="Tables/Degree_eval_allfeats.txt",sep="\t",row.names=FALSE)
  
  sub_colors_met<-{}
  #if(is.na(sigfeats)==FALSE){
  
  if(typeof(sigfeats)!="logical"){
    sub_colors_met<-colors_met_all[sigfeatsind,]
    write.table(sub_colors_met,file="Tables/Degree_eval_selectfeats.txt",sep="\t",row.names=FALSE)
  }
  
  ####saveMET,file="MET.Rda")
  
  
  
  return(list(all=colors_met_all,sigfeats=sub_colors_met))
}

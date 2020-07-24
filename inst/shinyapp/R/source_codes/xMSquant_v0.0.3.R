library(tidyr)
library(parallel)
library(ggplot2)
library(genefilter)
library(raster)
library(RColorBrewer)
#library(xmsPANDA)
source("R/source_codes/xmsPANDA_v1.0.8.43.R")
library(gplots)
theme_set(theme_classic())

find.Overlapping <- function (dataA, dataB, mz.thresh = 10, time.thresh = 30, best_one=FALSE) 
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

generate_distribution<-function(targeted_table,seq,outloc,groupcheck=FALSE,targetID=NA,min_num_nonmissing=3){
  
  #summarize the number of batches
  num_batches<- length(unique(seq[,4]))
  
  #check batchwise distribution
  pdf(paste(outloc,"batchwise_intensity_distribution.pdf",sep="/"))
  for(ii in 1:dim(targeted_table)[1]){
    
    inten_study<-c()
    batch_study<-c()
    inten_qstd<-c()
    batch_qstd<-c()
    for(jj in 1:num_batches){
      study_list=as.character(seq[seq[,4]==jj & grepl('study',seq[,3],ignore.case = TRUE),1])
      qstd_list=as.character(seq[seq[,4]==jj & grepl('qstd',seq[,3],ignore.case = TRUE),1])
      tmp_intensity_study<-as.numeric(targeted_table[ii,study_list])
      tmp_intensity_study<-tmp_intensity_study[!tmp_intensity_study==0]
      tmp_intensity_qstd<-as.numeric(targeted_table[ii,qstd_list])
      tmp_intensity_qstd<-tmp_intensity_qstd[!tmp_intensity_qstd==0]
      if(length(tmp_intensity_study)>=min_num_nonmissing){
        inten_study<-c(inten_study,tmp_intensity_study)
        batch_study<-c(batch_study,rep(paste("Batch",jj,sep=""),times=length(tmp_intensity_study)))
        inten_qstd<-c(inten_qstd,tmp_intensity_qstd)
        batch_qstd<-c(batch_qstd,rep(paste("Batch",jj,sep=""),times=length(tmp_intensity_qstd)))
      }
    }
    if(length(inten_study)>0){
      plotdata_study<-data.frame(inten_study,batch_study)
      plotdata_study$batch_study<-factor(plotdata_study$batch_study,levels=as.character(unique(plotdata_study$batch_study)))
      plotdata1<-cbind(plotdata_study,"black")
      colnames(plotdata1)<-c("inten","batch","col")
      if(length(inten_qstd)>0){
        plotdata_qstd<-data.frame(inten_qstd,batch_qstd)
        plotdata_qstd$batch_qstd<-factor(plotdata_qstd$batch_qstd,levels=as.character(unique(plotdata_qstd$batch_qstd)))
        plotdata2<-cbind(plotdata_qstd,"red")
        colnames(plotdata2)<-c("inten","batch","col")
        plotdata<-rbind(plotdata1,plotdata2)
      }else{
        plotdata<-plotdata1
      }
      
      if(length(grep("time_error",colnames(targeted_table)))>0){
        subtitle=paste(
          paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
          paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
          paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm","  ","time error:",round(targeted_table[ii,"time_error"],1),"s",sep=""),
          sep="\n")
      }else{
        subtitle=paste(
          paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
          paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
          paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm",sep=""),
          sep="\n")
      }
      
      #draw the plot
      print(ggplot() +
              geom_boxplot(data=plotdata_study, aes(x=batch_study, y=inten_study,fill=batch_study)) + 
              geom_point(data=plotdata, aes(x=batch, y=inten,fill=batch,color=col)) +
              scale_y_continuous(name="intensity") +
              scale_x_discrete(name="",limits = rev(levels(plotdata_study$batch_study)))+
              coord_flip() +
              labs(subtitle=subtitle, 
                   title=as.character(targeted_table[ii,"Chemical_name"])) +
              theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                    axis.line.y = element_line(size = 0.5, colour = "black"),
                    axis.line = element_line(size=1, colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.title=element_text(size = 20),
                    legend.title=element_text(size=13),
                    legend.text=element_text(size=10),
                    text=element_text(size = 16),
                    #plot.margin=margin(10,5,10,1),
                    axis.text.x=element_text(colour="black", size = 10),
                    axis.text.y=element_text(colour="black", size = 10))+
              scale_colour_manual(name ='Category', 
                                  values =c('black'='black','red'='red'), labels = c('Sample','QSTD'))+
              guides(fill="none"))
    }
    
  }
  dev.off()
  
  
  
  #check overall distribution 
  pdf(paste(outloc,"overall_intensity_distribution.pdf",sep="/"))
  for(ii in 1:dim(targeted_table)[1]){
    study_list=as.character(seq[grepl('study',seq[,3],ignore.case = TRUE),1])
    qstd_list=as.character(seq[grepl('qstd',seq[,3],ignore.case = TRUE),1])
    tmp_intensity_study<-as.numeric(targeted_table[ii,study_list])
    tmp_intensity_study<-tmp_intensity_study[!tmp_intensity_study==0]
    tmp_intensity_qstd<-as.numeric(targeted_table[ii,qstd_list])
    tmp_intensity_qstd<-tmp_intensity_qstd[!tmp_intensity_qstd==0]
    if(length(tmp_intensity_study)>=min_num_nonmissing){
      plotdata_study<-cbind(as.data.frame(tmp_intensity_study),"black")
      colnames(plotdata_study)<-c("intensity","col")
      if(length(tmp_intensity_qstd)>0){
        plotdata_qstd<- cbind(as.data.frame(tmp_intensity_qstd),"red")
        colnames(plotdata_qstd)<-c("intensity","col")
        plotdata<-rbind(plotdata_study,plotdata_qstd)
      }else{
        plotdata<-plotdata_study
      }
      
      if(any(is.na(targetID))){
        
        if(length(grep("time_error",colnames(targeted_table)))>0){
          subtitle=paste(
            paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
            paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
            paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm","  ","time error:",round(targeted_table[ii,"time_error"],1),"s",sep=""),
            sep="\n")
        }else{
          subtitle=paste(
            paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
            paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
            paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm",sep=""),
            sep="\n")
        }
        
        
        print(ggplot() +
                geom_boxplot(data=plotdata_study, aes(x="", y=intensity)) + 
                geom_point(data=plotdata, aes(x="", y=intensity,color=col)) +
                scale_y_continuous(name="intensity") +
                scale_x_discrete(name="")+
                coord_flip() +
                labs(subtitle=subtitle, 
                     title=as.character(targeted_table[ii,"Chemical_name"])) +
                theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                      axis.line.y = element_line(size = 0.5, colour = "black"),
                      axis.line = element_line(size=1, colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      plot.title=element_text(size = 20),
                      legend.title=element_text(size=13),
                      legend.text=element_text(size=10),
                      text=element_text(size = 16),
                      #plot.margin=margin(10,5,10,1),
                      axis.text.x=element_text(colour="black", size = 10),
                      axis.text.y=element_text(colour="black", size = 10))+
                scale_colour_manual(name ='Category', 
                                    values =c('black'='black','red'='red'), labels = c('Sample','QSTD'))+
                guides(fill="none"))
        
      }else{
        if(!any(as.character(seq[,2])%in%targetID)){
          stop("There is none of targetIDs that can find a match in class label file.", call.=TRUE)
        }
        if(length(targetID)>10){
          stop("The maximum number of allowed targetID is 10.", call.=TRUE)
        }
        plotdata_target<-seq[as.character(seq[,2])%in%targetID,c(1,2)]
        plotdata_target<-cbind(plotdata_target,as.numeric(targeted_table[ii,as.character(plotdata_target[,1])]))
        colnames(plotdata_target)<-c("FileName","SampleID","intensity")
        plotdata_target<-plotdata_target[!plotdata_target$intensity==0,]
        if(length(plotdata_target$intensity)>0){
          col=c(palette()[-c(1,2)],c("darkorange1","brown","gold1","deeppink1"))
          col<-col[1:length(plotdata_target$intensity)]
          plotdata_target<-cbind(plotdata_target,col)
          plotdata_target$col<-factor(plotdata_target$col,levels=as.character(unique(plotdata_target$col)))
          col_len<-c('black','red',as.character(plotdata_target$col))
          names(col_len)<-c('black','red',as.character(plotdata_target$col))
          plotdata[,"size"]=1
          plotdata=rbind(plotdata,cbind(plotdata_target[,c("intensity","col")],size=8))
          
          if(length(grep("time_error",colnames(targeted_table)))>0){
            subtitle=paste(
              paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
              paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
              paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm","  ","time error:",round(targeted_table[ii,"time_error"],1),"s",sep=""),
              sep="\n")
          }else{
            subtitle=paste(
              paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
              paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
              paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm",sep=""),
              sep="\n")
          }
          
          print(ggplot() +
                  geom_boxplot(data=plotdata_study, aes(x="", y=intensity)) + 
                  geom_point(data=plotdata, aes(x="", y=intensity,color=col,size=size)) +
                  scale_y_continuous(name="intensity") +
                  scale_x_discrete(name="")+
                  coord_flip() +
                  labs(subtitle=subtitle, 
                       title=as.character(targeted_table[ii,"Chemical_name"])) +
                  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                        axis.line.y = element_line(size = 0.5, colour = "black"),
                        axis.line = element_line(size=1, colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        plot.title=element_text(size = 20),
                        legend.title=element_text(size=13),
                        legend.text=element_text(size=10),
                        text=element_text(size = 16),
                        #plot.margin=margin(10,5,10,1),
                        axis.text.x=element_text(colour="black", size = 10),
                        axis.text.y=element_text(colour="black", size = 10))+
                  scale_colour_manual(name ='Category', 
                                      values =col_len, labels = c('Sample','QSTD',as.character(plotdata_target$SampleID)))+
                  guides(fill="none",size="none"))
        }else{
          
          if(length(grep("time_error",colnames(targeted_table)))>0){
            subtitle=paste(
              paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
              paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
              paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm","  ","time error:",round(targeted_table[ii,"time_error"],1),"s",sep=""),
              sep="\n")
          }else{
            subtitle=paste(
              paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
              paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
              paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm",sep=""),
              sep="\n")
          }
          
          print(ggplot() +
                  geom_boxplot(data=plotdata_study, aes(x="", y=intensity)) + 
                  geom_point(data=plotdata, aes(x="", y=intensity,color=col)) +
                  scale_y_continuous(name="intensity") +
                  scale_x_discrete(name="")+
                  coord_flip() +
                  labs(subtitle=subtitle, 
                       title=as.character(targeted_table[ii,"Chemical_name"])) +
                  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                        axis.line.y = element_line(size = 0.5, colour = "black"),
                        axis.line = element_line(size=1, colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        plot.title=element_text(size = 20),
                        legend.title=element_text(size=13),
                        legend.text=element_text(size=10),
                        text=element_text(size = 16),
                        #plot.margin=margin(10,5,10,1),
                        axis.text.x=element_text(colour="black", size = 10),
                        axis.text.y=element_text(colour="black", size = 10))+
                  scale_colour_manual(name ='Category', 
                                      values =c('black'='black','red'='red'), labels = c('Sample','QSTD'))+
                  guides(fill="none"))
        }
      }
    }
  }
  dev.off()
  
  
  #check group intensity distribution
  if(groupcheck==TRUE){
    if(length(grep("group",colnames(seq),ignore.case =TRUE))==0){
      stop("There is no column called 'group' in class label file. Please assign 'Group' to the name of the column with group label in class label file", call.=TRUE)
    }
    colnames(seq)[grep("group",colnames(seq),ignore.case =TRUE)]='Group'
    seq_group=seq[!is.na(seq$Group),]
    group=unique(as.character(seq_group$Group))
    pdf(paste(outloc,"groupwise_intensity_distribution.pdf",sep="/"))
    for(ii in 1:dim(targeted_table)[1]){
      inten_group<-c()
      label_group<-c()
      name_group<-c()
      for(jj in 1:length(group)){
        group_list=as.character(seq_group[seq_group$Group==group[jj],1])
        tmp_intensity_group<-as.numeric(targeted_table[ii,group_list])
        group_list<-group_list[!tmp_intensity_group==0]
        tmp_intensity_group<-tmp_intensity_group[!tmp_intensity_group==0]
        if(length(tmp_intensity_study)>=min_num_nonmissing){
          inten_group<-c(inten_group,tmp_intensity_group)
          label_group<-c(label_group,rep(group[jj],length(tmp_intensity_group)))
          name_group<-c(name_group,group_list)
        }
      }
      if(length(inten_group)>0){
        plotdata_group<-data.frame(name_group,inten_group,label_group)
        plotdata_group$label_group<-factor(plotdata_group$label_group,levels=as.character(unique(plotdata_group$label_group)))
        plotdata<-cbind(plotdata_group,"black")
        colnames(plotdata)<-c("name","inten","group","col")
        plotdata$col<-as.character(plotdata$col)
        plotdata<-merge(plotdata,seq[,c(1,2)],by.x='name',by.y='FileName',all.x=TRUE)
        
        colnames(plotdata)[ncol(plotdata)]="sampleID"
        
        if(!any(is.na(targetID))){
          if(!any(as.character(seq_group[,2])%in%targetID)){
            stop("There is none of targetIDs that can find a match in class label file.", call.=TRUE)
          }
          if(length(targetID)>10){
            stop("The maximum number of allowed targetID is 10.", call.=TRUE)
          }
          if(any(plotdata$name%in%as.character(seq_group[as.character(seq_group[,2])%in%targetID,1]))){
            plotdata=cbind(plotdata,size=1)
            plotdata[plotdata$name%in%as.character(seq_group[as.character(seq_group[,2])%in%targetID,1]),"size"]=8
            col=c(palette()[-c(1,2)],c("darkorange1","brown","gold1","deeppink1"))
            col<-col[1:dim(plotdata[plotdata$name%in%as.character(seq_group[as.character(seq_group[,2])%in%targetID,1]),])[1]]
            plotdata[plotdata$name%in%as.character(seq_group[as.character(seq_group[,2])%in%targetID,1]),"col"]=col
            plotdata$col<-factor(plotdata$col,levels=as.character(unique(plotdata$col)))
            col_len<-c('black',as.character(col))
            names(col_len)<-c('black',as.character(col))
            
            #draw the plot
            if(length(grep("time_error",colnames(targeted_table)))>0){
              subtitle=paste(
                paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
                paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
                paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm","  ","time error:",round(targeted_table[ii,"time_error"],1),"s",sep=""),
                sep="\n")
            }else{
              subtitle=paste(
                paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
                paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
                paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm",sep=""),
                sep="\n")
            }
            
            print(ggplot() +
                    geom_boxplot(data=plotdata, aes(x=group, y=inten,fill=group)) + 
                    geom_point(data=plotdata, aes(x=group, y=inten,fill=group,color=col,size=size)) +
                    scale_y_continuous(name="intensity") +
                    scale_x_discrete(name="",limits = rev(levels(plotdata_study$batch_study)))+
                    coord_flip() +
                    labs(subtitle=subtitle, 
                         title=as.character(targeted_table[ii,"Chemical_name"])) +
                    theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                          axis.line.y = element_line(size = 0.5, colour = "black"),
                          axis.line = element_line(size=1, colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          plot.title=element_text(size = 20),
                          legend.title=element_text(size=13),
                          legend.text=element_text(size=10),
                          text=element_text(size = 16),
                          #plot.margin=margin(10,5,10,1),
                          axis.text.x=element_text(colour="black", size = 10),
                          axis.text.y=element_text(colour="black", size = 10))+
                    scale_colour_manual(name ='Category', 
                                        values =col_len, labels = c('Sample',as.character(plotdata[plotdata$name%in%as.character(seq_group[as.character(seq_group[,2])%in%targetID,1]),"sampleID"])))+
                    guides(fill="none",size="none"))
          }else{
            if(length(grep("time_error",colnames(targeted_table)))>0){
              subtitle=paste(
                paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
                paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
                paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm","  ","time error:",round(targeted_table[ii,"time_error"],1),"s",sep=""),
                sep="\n")
            }else{
              subtitle=paste(
                paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
                paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
                paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm",sep=""),
                sep="\n")
            }
            
            print(ggplot() +
                    geom_boxplot(data=plotdata, aes(x=group, y=inten,fill=group)) + 
                    geom_point(data=plotdata, aes(x=group, y=inten,fill=group,color=col)) +
                    scale_y_continuous(name="intensity") +
                    scale_x_discrete(name="",limits = rev(levels(plotdata_study$batch_study)))+
                    coord_flip() +
                    labs(subtitle=subtitle, 
                         title=as.character(targeted_table[ii,"Chemical_name"])) +
                    theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                          axis.line.y = element_line(size = 0.5, colour = "black"),
                          axis.line = element_line(size=1, colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          plot.title=element_text(size = 20),
                          legend.title=element_text(size=13),
                          legend.text=element_text(size=10),
                          text=element_text(size = 16),
                          #plot.margin=margin(10,5,10,1),
                          axis.text.x=element_text(colour="black", size = 10),
                          axis.text.y=element_text(colour="black", size = 10))+
                    scale_colour_manual(name ='Category', 
                                        values =c('black'='black'), labels = c('Sample'))+
                    guides(fill="none"))
          }
          
          
        }else{
          #draw the plot
          if(length(grep("time_error",colnames(targeted_table)))>0){
            subtitle=paste(
              paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
              paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
              paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm","  ","time error:",round(targeted_table[ii,"time_error"],1),"s",sep=""),
              sep="\n")
          }else{
            subtitle=paste(
              paste("mz:",round(targeted_table[ii,"mz"],5),"  ","time:",round(targeted_table[ii,"time"],1),sep=""),
              paste("reference mz:",round(targeted_table[ii,"ref_mz"],5),"  ","reference time:",round(targeted_table[ii,"ref_time"],1),sep=""),
              paste("mz error:",round(targeted_table[ii,"mz_error"]),"ppm",sep=""),
              sep="\n")
          }
          
          print(ggplot() +
                  geom_boxplot(data=plotdata, aes(x=group, y=inten,fill=group)) + 
                  geom_point(data=plotdata, aes(x=group, y=inten,fill=group,color=col)) +
                  scale_y_continuous(name="intensity") +
                  scale_x_discrete(name="",limits = rev(levels(plotdata_study$batch_study)))+
                  coord_flip() +
                  labs(subtitle=subtitle, 
                       title=as.character(targeted_table[ii,"Chemical_name"])) +
                  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                        axis.line.y = element_line(size = 0.5, colour = "black"),
                        axis.line = element_line(size=1, colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        plot.title=element_text(size = 20),
                        legend.title=element_text(size=13),
                        legend.text=element_text(size=10),
                        text=element_text(size = 16),
                        #plot.margin=margin(10,5,10,1),
                        axis.text.x=element_text(colour="black", size = 10),
                        axis.text.y=element_text(colour="black", size = 10))+
                  scale_colour_manual(name ='Category', 
                                      values =c('black'='black'), labels = c('Sample'))+
                  guides(fill="none"))
        }
        
      }
      
    }
    dev.off()
  }
  
}


calculate_concentration<-function(targeted_table,ref,seq,outloc,missing= 0.3,batch_qstd_missing=0.4){
  
  count=max(grep("_error",colnames(targeted_table)))
  
  #summarize the number of batches
  num_batches<- length(unique(seq[,4]))
  
  #create list object with each batch in slot
  qstd_sample_mapping<- seq[grepl('qstd',seq[,3],ignore.case = TRUE),]
  qstd_int<- targeted_table[,as.character(qstd_sample_mapping[,1])]
  
  
  #show error if number of reference samples is not 6x number of batches
  #if(ncol(qstd_int)/6 != length(unique(sample.mapping[,4]))){
  #	stop("Number of reference samples is not multiple of number of batches", call.=TRUE)
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
  
}


draw_keggmap<-function(targeted_table,outloc,foldchange.thresh=2,minhit=3,highcolor='red',lowcolor='blue',percent_node=0.2){
  
  if(length(grep('keggid',colnames(targeted_table),ignore.case=TRUE))==0){
    stop("No KEGG ID was detected in reference list. To draw kegg map, please add one more column named 'KEGGID' to reference list for KEGG ID.")
  }
  
  iddata=targeted_table[,"KEGGID"]
  
  targeted_table[,'map'] <- unlist(mclapply(1:length(iddata), function(j) {
    API=paste("http://rest.kegg.jp/get/",as.character(iddata[j]),sep="")
    content<-try(readLines(API),silent=TRUE)
    if (class(content) == "try-error") {
      return(NA)
    }else{
      return(paste(gsub('^.*(map[0-9]*).*$','\\1',content[grep('map[0-9]',content)],perl=TRUE),collapse=","))
    }
  },mc.cores=detectCores()))
  
  totalmap=unique(unlist(strsplit(paste(targeted_table[,"map"][grep("map",targeted_table[,"map"])],collapse=","), ",")))
  mapdata=data.frame(mapID=totalmap)
  
  for(ii in 1:length(totalmap)){
    mapdata[ii,"ref_mz_time"]=paste(paste(round(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"ref_mz"],5),round(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"ref_time"],1),sep="_"),collapse="/")
    mapdata[ii,"mz_time"]=paste(paste(round(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"mz"],5),round(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"time"],1),sep="_"),collapse="/")
    mapdata[ii,"keggid"]=paste(as.character(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"KEGGID"]),collapse="/")
    mapdata[ii,"foldchange"]=paste(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"foldchange"],collapse="/")
  }

  write.table(mapdata, paste(outloc,"keggmapdata.txt",sep="/"), sep= "\t", row.names=FALSE)
  
  suppressWarnings(dir.create(paste(outloc,"KEGGmaps",sep="/")))
  
  print_png <- function(i,mapdata,foldchange.thresh,lowcolor,highcolor,minhit,outloc){
    
    #load the required packages
    library(pathview)
    
    if(as.character(mapdata[i,1])%in%c('map01100','map01110','map01130','map01120','map00121','map01060')){
      
    }else{
      cpd.data=as.numeric(unlist(strsplit(as.character(mapdata[i,"foldchange"]), "/")))
      names(cpd.data)=unlist(strsplit(as.character(mapdata[i,"keggid"]), "/"))
      cpd.data=cpd.data[!duplicated(cpd.data)]
      if(length(cpd.data)<minhit){
        
      }else{
        log.cpd.data <- log(cpd.data,base = foldchange.thresh)
        pathwayid=gsub("map","",as.character(mapdata[i,"mapID"]))
        suffix=paste("highlighted","_",length(log.cpd.data),'hits',sep="")
        setwd(paste(outloc,"KEGGmaps",sep="/"))
        try(pathview(gene.data =NULL, cpd.data=log.cpd.data, pathway.id=pathwayid, cpd.idtype='kegg', same.layer=TRUE, plot.col.key =TRUE,new.signature=FALSE,low = list(cpd = lowcolor), mid =list(cpd = "grey"),high=list(cpd = highcolor),out.suffix=suffix,species = "ko"),silent=TRUE)
        unlink(paste(paste("ko",pathwayid,sep=""),"png",sep="."), recursive = FALSE)
        unlink(paste(paste("ko",pathwayid,sep=""),"xml",sep="."), recursive = FALSE)
        print(getwd())
      }
    }
  }
  
  nu_cores <- detectCores()
  cl <- makeCluster(ceiling(nu_cores*as.numeric(percent_node)))
  parLapply(cl,1:nrow(mapdata),print_png,mapdata=mapdata,foldchange.thresh=foldchange.thresh,lowcolor=lowcolor,highcolor=highcolor,minhit=minhit,outloc=outloc)
  stopCluster(cl)
  
}



quant<- function(Xmat=NA,Ymat=NA,Wmat=NA,Zmat=NA,feature_table,class_file,ref_list,foldchange_list,outloc,
                    num_replicates=1,
                    summarize_replicates=FALSE,
                    rep.max.missing.thresh=0.5,
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
                    lowcolor='blue'
) {
  
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
  
  suppressWarnings(dir.create(outloc))
  setwd(outloc)
  
  if(summarize_replicates==TRUE){
    capture.output(avg_feature_table<-data_preprocess(X=feature.table,feature_table_file=NA,parentoutput_dir=outloc,class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method=summary.method,all.missing.thresh=NA,group.missing.thresh=NA,log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,samplermindex=NA, rep.max.missing.thresh=rep.max.missing.thresh,summary.na.replacement="zeros"), file='/dev/null')
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
    suppressWarnings(dir.create(outloc1))
    generate_distribution(targeted_table=r_targeted_table,seq=sample.mapping,outloc=outloc1,groupcheck=groupcheck,targetID=targetID,min_num_nonmissing=min_num_nonmissing)
    print("Step 1 complete")
  }
  
  #calculate the concentration
  if(length(grep("2",steps))!=0){
    print("2. Calculating the concentration for all samples.")
    outloc2=paste(outloc,"/step2",sep="")
    suppressWarnings(dir.create(outloc2))
    calculate_concentration(targeted_table=r_targeted_table,ref=ref_table,seq=sample.mapping,outloc=outloc2)
    print("Step 2 complete")
  }
  
  #draw kegg map
  if(length(grep("3",steps))!=0){
    print("3. Pulling the KEGG maps out from KEGG website.")
    outloc3=paste(outloc,"/step3",sep="")
    suppressWarnings(dir.create(outloc3))
    draw_keggmap(targeted_table=r_targeted_table,outloc=outloc3,foldchange.thresh=foldchange_thresh,minhit=minhit,highcolor=highcolor,lowcolor=lowcolor,percent_node=percent_node)
    print("Step 3 complete")
  }
  
  #suppressWarnings(dir.create(paste(outloc,"KEGGmaps",sep="/")))
  #try(file.rename(paste(getwd(),list.files(pattern=".png"),sep="/"),paste(paste(outloc,"KEGGmaps",sep="/"),list.files(pattern=".png"),sep="/")),silent=TRUE)
  
}


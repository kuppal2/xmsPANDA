generate_distribution_concentration <-
function(targeted_table,seq,outloc,groupcheck=FALSE,targetID=NA,min_num_nonmissing=3){
  
  #summarize the number of batches
  num_batches<- length(unique(seq[,4]))
  
  
  #check batchwise distribution
  pdf(paste(outloc,"batchwise_concentration_distribution.pdf",sep="/"))
  for(ii in 1:dim(targeted_table)[1]){
    
    inten_study<-c()
    batch_study<-c()
    inten_qstd<-c()
    batch_qstd<-c()
    for(jj in 1:num_batches){
      study_list=as.character(seq[seq[,4]==jj & grepl('study',seq[,3],ignore.case = TRUE),1])
      qstd_list=as.character(seq[seq[,4]==jj & grepl('qstd',seq[,3],ignore.case = TRUE),1])
      study_list<-which(colnames(targeted_table)%in%study_list)
      qstd_list<-which(colnames(targeted_table)%in%qstd_list)
      
      
      
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
              stat_boxplot(geom='errorbar',width=0.2) +
              geom_boxplot(data=plotdata_study, aes(x=batch_study, y=inten_study,fill=batch_study)) +
              geom_point(data=plotdata, aes(x=batch, y=inten,fill=batch,color=col)) +
              scale_y_continuous(name="concentration") +
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
  pdf(paste(outloc,"overall_concentration_distribution.pdf",sep="/"))
  for(ii in 1:dim(targeted_table)[1]){
    study_list=as.character(seq[grepl('study',seq[,3],ignore.case = TRUE),1])
    qstd_list=as.character(seq[grepl('qstd',seq[,3],ignore.case = TRUE),1])
    
    study_list<-which(colnames(targeted_table)%in%study_list)
    qstd_list<-which(colnames(targeted_table)%in%qstd_list)
    
    
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
                stat_boxplot(geom='errorbar',width=0.2) +
                geom_boxplot(data=plotdata_study, aes(x="", y=intensity)) +
                geom_point(data=plotdata, aes(x="", y=intensity,color=col)) +
                scale_y_continuous(name="concentration") +
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
                  stat_boxplot(geom='errorbar',width=0.2) +
                  geom_boxplot(data=plotdata_study, aes(x="", y=intensity)) +
                  geom_point(data=plotdata, aes(x="", y=intensity,color=col,size=size)) +
                  scale_y_continuous(name="concentration") +
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
                  stat_boxplot(geom='errorbar',width=0.2) +
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
      #stop("There is no column called 'group' in class label file. Please assign 'Group' to the name of the column with group label in class label file", call.=TRUE)
      Group<-rep(1,nrow(seq))
      seq<-cbind(seq,Group)
    }
    print(head(seq))
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
        
        group_list<-which(colnames(targeted_table)%in%group_list)
        
        
        
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
                    stat_boxplot(geom='errorbar',width=0.2) +
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
                    stat_boxplot(geom='errorbar',width=0.2) +
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
                  stat_boxplot(geom='errorbar',width=0.2) +
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

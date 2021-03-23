get_spaghettiplots <-
function(X,Y,feature_table_file,parentoutput_dir,class_labels_file,boxplot.col.opt="journal",alphacol=0.3,newdevice=TRUE,cex.plots=0.8,
                             replace.by.NA=FALSE,pairedanalysis=FALSE,filename="",ylabel="Intensity",
                             alphabetical.order=FALSE,name=NA,add.jitter=TRUE,add.pvalues=TRUE,class.levels=NA,fill.plots=FALSE,connectpairedsamples=FALSE,
                             boxplot.type="ggplot",study.design=c("multiclass","onewayanova","twowayanova","onewayanovarepeat",
                                                                  "twowayanovarepeat"),multiple.figures.perpanel=TRUE,
                             ggplot.type1=TRUE,replace.outliers=FALSE,plot.height=8,plot.width=8,...)
{
  
  analysistype=study.design[1]
  
  if(boxplot.type=="ggplot"){
    suppressMessages(library(ggplot2))
  }else{
    suppressMessages(library(ggpubr))
  }
  
  # multiple.figures.perpanel=TRUE
  
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
  rm(X)
  rm(Y)
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
        
        data_matrix<-cbind(mz,time,data_matrix[,-check_ind])
        names_with_mz_time=cbind(Name,mz,time)
        
        names_with_mz_time<-as.data.frame(names_with_mz_time)
        data_matrix<-as.data.frame(data_matrix)
        
        #write.table(names_with_mz_time,file="Stage1/Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        
        write.table(names_with_mz_time,file="Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        
        
      }else{
        
        if(length(check_names1)>0 & length(check_names2)>0){
          
          check_ind<-gregexpr(cnames,pattern="^name$")
          check_ind<-which(check_ind>0)
          Name<-as.character(data_matrix[,check_ind])
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
  }
  
  
  
  dir.create(parentoutput_dir,showWarnings = FALSE)
  setwd(parentoutput_dir)
  
  data_m<-data_matrix[,-c(1:2)]
  
  data_m<-as.matrix(data_m)
  
  mzvec<-data_matrix[,1]
  timevec<-data_matrix[,2]
  goodfeats<-data_m			
  rm(data_m)
  
  #library(extrafont)
  #loadfonts()
  
  multiple.groups=FALSE
  
  #c("multiclass","onewayanova","twowayanova","onewayanovarepeat","twowayanovarepeat")
  
  if(study.design=="onewayanovarepeat" | study.design=="twowayanovarepeat" | study.design=="twowayrepeat" | study.design=="onewayrepeat"){
    
    pairedanalysis=TRUE
  }
  
  if(pairedanalysis==TRUE){
    
    print("Using column 2 as subject identifiers")
    paireddesign=classlabels[,2]
    
    classlabels<-classlabels[,-c(2)]
  }
  
  
  
  if(dim(classlabels)[2]>2){
    
    
    if(study.design=="twowayanova" | study.design=="twowayanovarepeat" | study.design=="twoway" | study.design=="twowayrepeat"){
      # print("More than two columns found in the class labels file. ")
      
      
      Class<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
      
      multiple.groups=TRUE
      
    }else{
      
      if(study.design=="onewayanova" | study.design=="onewayanovarepeat" | study.design=="oneway"){
        #     print("More than two columns found in the class labels file. ")
        
        
        Class<-classlabels[,2] #classlabels[,2]:classlabels[,3]
        
        
      }else{
        Class<-classlabels[,2]
        
      }
      
    }
    
  }else{
    
    Class<-classlabels[,2]
  }
  
  
  #par(mfrow=c(2,2),family="sans",cex=cex.plots)
  
  if(alphabetical.order==FALSE){
    Class <- factor(Class, levels=unique(Class))
  }
  
  class_levels<-levels(as.factor(Class))
  
  class_labels_levels<-levels(as.factor(Class))
  ordered_labels<-Class
  
  class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
  
  if(is.na(boxplot.col.opt)==TRUE){
    
    col_vec<-rep(c("white"),length(class_labels_levels))
    boxplot.col.opt<-col_vec
  }
  
  if(boxplot.col.opt=="default"){
    
    col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
               "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
               "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
               "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
    
  }else{ 
    if(boxplot.col.opt=="topo"){
      #col_vec<-topo.colors(256) #length(class_labels_levels)) 
      
      #col_vec<-col_vec[seq(1,length(col_vec),)]
      
      col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
    }else{
      if(boxplot.col.opt=="heat"){
        #col_vec<-heat.colors(256) #length(class_labels_levels))
        
        col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
      }else{
        if(boxplot.col.opt=="rainbow"){
          #col_vec<-heat.colors(256) #length(class_labels_levels))
          col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
          
          #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
        }else{
          
          if(boxplot.col.opt=="terrain"){
            #col_vec<-heat.colors(256) #length(class_labels_levels))
            #col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
            
            col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
          }else{
            
            if(boxplot.col.opt=="colorblind"){
              #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
              # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
              
              if(length(class_labels_levels)<9){
                
                col_vec <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "#E64B35FF", "grey57")
                
              }else{
                
                #   col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                           "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                
              }
              
              
            }else{
              
              check_brewer<-grep(pattern="brewer",x=boxplot.col.opt)
              
              if(length(check_brewer)>0){
                
                boxplot.col.opt=gsub(x=boxplot.col.opt,pattern="brewer.",replacement="")
                col_vec <- colorRampPalette(brewer.pal(10, boxplot.col.opt))(length(class_labels_levels))
                
              }else{
                
                if(boxplot.col.opt=="journal"){
                  
                  col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                             "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                             "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                             
                             "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                             "#E64B3519","#4DBBD519","#631879E5","grey75")
                  if(length(class_labels_levels)<8){
                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75")
                    
                    #col_vec2<-brewer.pal(n = 8, name = "Dark2")
                    
                  }else{
                    if(length(class_labels_levels)<=28){
                      # col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "grey75","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666","#1B9E77", "#7570B3", "#E7298A", "#A6761D", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
                      
                      col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                 "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                 "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                 
                                 "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF", "#8BD76BFF",
                                 "#E64B3519","#9DBBD0FF","#631879E5","#666666","grey75")
                      
                    }else{
                      
                      
                      
                      
                      colfunc <-colorRampPalette(c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75"));col_vec<-colfunc(length(class_labels_levels))
                      
                      col_vec<-col_vec[sample(col_vec)]
                      
                      
                    }
                  }
                  
                  
                  
                }else{
                  col_vec <-boxplot.col.opt
                  col_vec <- rep(col_vec,length(class_labels_levels))
                }
                
              }
              
            }
          }
          
          
        }
        
      }
      
    }	
  }
  
  
  
  ordered_labels={}
  num_samps_group<-new("list")
  num_samps_group[[1]]<-0
  groupwiseindex<-new("list")
  groupwiseindex[[1]]<-0
  
  for(c in 1:length(class_labels_levels))
  {
    
    classlabels_index<-which(Class==class_labels_levels[c])
    #ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
    num_samps_group[[c]]<-length(classlabels_index)
    groupwiseindex[[c]]<-classlabels_index
  }
  
  sampleclass<-{}
  patientcolors<-{}
  
  if(length(mzvec)>4){
    max_per_row<-3
    
    
    par_rows<-ceiling(9/max_per_row)
    
  }else{
    max_per_row<-length(mzvec)
    par_rows<-1
  }
  
  
  #name=goodfeats_name,
  #class_labels_levels<-paste("x",seq(1,length(class_labels_levels)),sep="")
  
  file_ind<-0
  boxplots_fname<-paste(filename,".pdf",sep="")
  #tiff(boxplots_fname, width=plots.width,height=plots.height,res=plots.res, compression="lzw")
  
  #tiff(boxplots_fname, width=2000,height=3000,res=plots.res, compression="lzw")
  
  if(newdevice==TRUE){
    pdf(boxplots_fname,width=plot.width,height=plot.height)
  }
  #par(mfrow=c(par_rows,max_per_row))
  #save(Class,goodfeats,class_labels_levels,file="debug1.Rda")
  
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  
  
  # text(5,9,"Description:",font=2,col="blue")
  
  text(5, 8, "This PDF includes spaghetti plots for each variable representing\n the longitudinal variation in
  individual subjects over time.
       \n\n\n Note: The panels are grouped by factor 1 (e.g. group).",cex=1.5,font=2)
  
  
  temp_dm<-cbind(as.character(Class),(t(goodfeats)))
  temp_dm<-as.data.frame(temp_dm)
  colnames(temp_dm)<-c("Class",rownames(goodfeats))
  class_levels<-levels(as.factor(temp_dm$Class))
  
  #save(temp_dm,file="temp_dm.Rda")
  #keeps the class order same as in the input file; avoids arrangement by alphabetical order
  if(alphabetical.order==FALSE){
    temp_dm$Class <- factor(temp_dm$Class, levels=unique(temp_dm$Class))
  }
  
  #save(temp_dm,file="temp_dm1.Rda")
  
  if(is.na(class.levels)==FALSE){
    
    match_test<-match(class.levels,levels(as.factor(temp_dm$Class)))
    if(length(which(is.na(match_test)==TRUE))<1){
      
      if(alphabetical.order==FALSE){
        temp_dm$Class <- factor(temp_dm$Class, levels=unique(class.levels))
      }
    }else{
      stop(paste("User defined class.levels ", paste(class.levels,sep=" ",collapse=""), " do not match the levels in the class labels matrix, ",paste(class_levels,sep=" ",collapse=""),sep=""))    
    }
  }
  
  #par(mfrow=c(1,1),family="sans",cex=cex.plots)
  

    {
    
    plot_res<-lapply(1:dim(goodfeats)[1],function(m)
    {
      
      if(m%%9==0){
        
        file_ind<-file_ind+1
        boxplots_fname<-paste("boxplots_file",file_ind,".tiff",sep="")
        
      }
      
      round_mzval<-mzvec[m] #sprintf("%.4f",mzvec[m])
      
      round_timeval<-timevec[m] #sprintf("%.1f",timevec[m])
      
      if(is.na(name[1])==TRUE){            
        
        if(length(check_names)>0){	
          if(check_names==1){
            
            mzname<-as.character(names_with_mz_time[m,1])
          }else{  
            
            mzname<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
          }
          
        }else{
          
          mzname<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
        }
      }else{
        
        mzname=as.character(name[m])
      }
      if(length(class_labels_levels)>=2)
      {
        if(length(class_labels_levels)>=1)
        {
          t1<-table(sampleclass)
          cur_d<-new("list")
          feat_vec<-{}
          class_vec<-{}
          id_vec<-{}
          
          
          for(c in 1:length(class_labels_levels))
          {
            num_samps_group[[1]]<-t1[1]
            cvec<-as.vector(t(goodfeats[m,c(groupwiseindex[[c]])]))
            if(replace.outliers==TRUE){
              cvec<-replace_outliers(cvec,replace.by.NA)
            }
            cur_d[[c]]<-cvec
            feat_vec<-c(feat_vec,cvec)
            class_vec<-c(class_vec,rep(class_labels_levels[c],length(which(Class==class_labels_levels[c]))))
            id_vec<-c(id_vec,names(goodfeats[m,c(groupwiseindex[[c]])]))
          }
          
          
          temp_dm<-cbind(as.character(class_vec),as.vector(feat_vec),as.vector(id_vec))
          
          temp_dm<-as.data.frame(temp_dm)
          #save(temp_dm,class_labels_levels,goodfeats,Class,groupwiseindex,file="t1.Rda")
          temp_dm2<-temp_dm #[,c(1,(m+1))]
          temp_dm2<-as.data.frame(temp_dm2)
          colnames(temp_dm2)<-c("Class","Feature","SubjectID")
          temp_dm2$Feature<-as.numeric(as.character(temp_dm2$Feature))
          
          sum_yval1=summary(temp_dm2$Feature,na.rm=TRUE)
          max_yval1=max(temp_dm2$Feature,na.rm=TRUE)+(sum_yval1[5]-sum_yval1[2])
          
         
          
          fname1<-paste("temp_dm2",mzname,"A1.Rda")
          #save(Class,class.levels,temp_dm2,alphabetical.order,file=fname1)
          Class<-temp_dm2
          
          if(is.na(class.levels)==FALSE){
            
            match_test<-match(class.levels,levels(as.factor(temp_dm2$Class)))
            if(length(which(is.na(match_test)==TRUE))<1){
              
              if(alphabetical.order==FALSE){
                temp_dm2$Class <- factor(temp_dm2$Class, levels=unique(class.levels))
              }
            }else{
              stop(paste("User defined class.levels ", paste(class.levels,sep=" ",collapse=""), " do not match the levels in the class labels matrix, ",paste(class_levels,sep=" ",collapse=""),sep=""))
            }
          }
          
          
          if(multiple.groups==TRUE){
            Factor1<-gsub(temp_dm2$Class,pattern=":([\\w|\\W])*",replacement="",perl=TRUE)
            Factor2<-gsub(temp_dm2$Class,pattern="([\\w|\\W])*:",replacement="",perl=TRUE)
            
            
            fname1<-paste("temp_dm2",mzname,"A.Rda")
            #save(Class,Factor1,Factor2,temp_dm2,file=fname1)
            temp_dm2<-cbind(temp_dm2,Factor1,Factor2)
          }
          
          
          
          
          #  par(mfrow=c(2,2),family="sans",cex=cex.plots)
          
          w <- 0.1
          par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
          
          #    print(cex)
          ###save(temp_dm2,cur_d,mzname,boxplot.col.opt,col_vec,ylabel,Class,goodfeats,class_labels_levels,cex,boxplot.type,multiple.groups,add.jitter,add.pvalues,fill.plots,file="temp_dm2.Rda")
          
          #ID<-paste("SID",seq(1,nrow(temp_dm2)),sep="")
          temp_dm2$SubjectID<-gsub(temp_dm2$SubjectID,pattern="[a-z|0-9|A-Z]*_",replacement="")
          
         # temp_dm2<-cbind(ID,temp_dm2)
          temp_dm2<-as.data.frame(temp_dm2)
            
            
           
              
              suppressMessages(library(ggplot2))
              if(multiple.groups==FALSE){
                p <- ggplot(temp_dm2, aes(y=as.numeric(Feature),x=Class,group=SubjectID)) + geom_line() + labs(title=mzname)
                
              }else{
                
                fname1<-paste("temp_dm2",mzname,".Rda")
                
                #save(temp_dm2,mzname,class_labels_levels,col_vec,file=fname1)
                
                 #p <- ggplot(temp_dm2, aes(y=as.numeric(Feature),x=Factor1,fill=Factor2)) + labs(title=mzname) + facet_wrap(~Factor2, scale="free")  
                if(ggplot.type1==TRUE){
                  p <- ggplot(temp_dm2, aes(y=as.numeric(Feature),x=Factor1,group=SubjectID,color=SubjectID)) + geom_point()+geom_line(size=1.25)+labs(title=mzname) + facet_wrap(~Factor2, scale="free_x")  
                }else{
                  p <- ggplot(temp_dm2, aes(y=as.numeric(Feature),x=Factor2,group=SubjectID,color=SubjectID))+ geom_point()+geom_line(size=1.25)+labs(title=mzname) + facet_wrap(~Factor1, scale="free_x")  
                  
                  
                }
                
              }
              
             if(FALSE){
               
              p<-p+ labs(y=ylabel) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(),
                                                        panel.spacing=unit(1,"lines"),
                                                        axis.line = element_line(colour = "black",size=1),
                                                        axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=18*cex.plots,face="bold"),
                                                        plot.title = element_text(hjust = 0.5,size=18*cex.plots),
                                                        axis.ticks.length = unit(-0.05, "in"),
                                                        axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
                                                        axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
                                                        axis.ticks.x = element_blank(),
                                                        aspect.ratio = 1,
                                                        legend.background = element_rect(color = "black", fill = "white"),
                                                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                        strip.text.x = element_text(size = 14*cex.plots, colour = "black"),
                                                        strip.text = element_text(face="bold")) +scale_fill_manual(values=(col_vec[1:length(unique(temp_dm2$ID))]))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
             }
              
              p<-p+ labs(y=ylabel) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(),
                                                        panel.spacing=unit(1,"lines"),
                                                        axis.line = element_line(colour = "black",size=1),
                                                        axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=18*cex.plots,face="bold"),
                                                        plot.title = element_text(hjust = 0.5,size=18*cex.plots),
                                                       
                                                        legend.background = element_rect(color = "black", fill = "white"),
                                                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                        strip.text.x = element_text(size = 14*cex.plots, colour = "black"),
                                                       strip.text = element_text(face="bold")) +scale_fill_manual(values=(col_vec[1:length(unique(temp_dm2$ID))]))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
              p=p + theme(legend.title = element_text(size = 13*cex.plots),
                          legend.text = element_text(size=13*cex.plots))
              
              #save(p,file="p1.Rda")
            
             
            
            return(p)
            #print(p)
        }}
          
    })
    
    #save(plot_res,multiple.figures.perpanel,file="boxplot_plot.Rda")
    
    if(boxplot.type=="ggplot"){
      
      suppressMessages(library(ggpubr))
      
      
      
      if(length(plot_res)>0){
        if(multiple.figures.perpanel==FALSE){
          res<-lapply(seq(1,length(plot_res),1),function(i){
            p1=plot_res[[i]]
            figure<-ggpubr::ggarrange(p1,ncol = 1, nrow = 1,heights=c(6),width=c(10),legend=TRUE,align = c("hv"))
            return(figure)
          })
        }else{
          res<-lapply(seq(1,length(plot_res),4),function(i){
            p1=plot_res[[i]]
            p1=plot_res[[i]]
            p2={}
            p3={}
            p4={}
            
            if((i+1)<length(plot_res)){
              p2=plot_res[[i+1]]
            }
            if((i+2)<length(plot_res)){
              p3=plot_res[[i+2]]
            }
            if((i+3)<length(plot_res)){
              
              p4=plot_res[[i+3]]
            }
            figure<-ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2,heights=c(4,4),width=c(6,6),legend=FALSE,align = c("hv"))
            # gg##save(res,file="t.pdf")
            #figure<-ggarrange(p1,ncol = 2, nrow = 2,heights=c(4,4),width=c(6,6),legend=FALSE,align = c("hv"))
            return(figure)
          })
          
        }
        
      }
      
      
      
      #save(res,file="res.Rda")
      
      for(i in 1:length(res)){
        p_fixed <- set_panel_size(res[[i]][[1]],
                                  width  = unit(plot.width-0.5, "in"),
                                  height = unit(plot.height-0.5, "in"))
        
        
        grid.newpage()
        grid.draw(p_fixed)
        
        #print(p_fixed) #res[[i]][[1]])
      }
      
      
      legend <- ggpubr::get_legend(plot_res[[i]])
      
      suppressMessages(library(grid))
      
      grid.newpage()
      grid.draw(legend)
      
      
    }
  }
  
  if(newdevice==TRUE){
    try(dev.off(boxplots_fname),silent=TRUE)
  }
  
  par(mfrow=c(1,1))
}

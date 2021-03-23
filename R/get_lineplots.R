get_lineplots <-
function(X=NA,Y=NA,feature_table_file=NA,parentoutput_dir=NA,
                        class_labels_file=NA,alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,
                        point.cex.val=4,legendlocation="topright",pca.ellipse=TRUE,
                        ellipse.conf.level=0.95,filename="all",newdevice=FALSE,
                        lineplot.col.opt=c("journal"),ylabel="Intensity",error.bar=TRUE,
                        cex.plots=0.8,
                        lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                        timeseries.lineplots=FALSE,name=NA,study.design="oneway",
                        alphabetical.order=TRUE,output.format="pdf",
                        multiple.figures.perpanel=TRUE,sizeval=1.5,plot.height=8,plot.width=8)
{
  
  
  cex.val=cex.plots
  
  suppressMessages(library(ggpubr))
  suppressMessages(library(ggplot2))
  
  analysistype=study.design
  
  #lineplot.col.opt=tolower(lineplot.col.opt)
  if(length(lineplot.col.opt)==1){
  lineplot.col.opt=get_hexcolors_for_palettes(color.palette=lineplot.col.opt,alpha.col=alphacol[1])
  }
  
  if(is.na(parentoutput_dir)==TRUE){
    
    parentoutput_dir=getwd()
  }
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
        
        #write.table(names_with_mz_time,file="Stage1/Name_mz_time_mapping.txt",sep="\t",row.names=FALSE)
        
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
  }
  
  #addition end
  
  sample.col.opt=lineplot.col.opt
  
  par(mfrow=c(1,1),family="sans",cex=cex.val)
  
  #print("here")
  dir.create(parentoutput_dir,showWarnings = FALSE)
  setwd(parentoutput_dir)
  
  data_m<-data_matrix[,-c(1:2)]
  
  data_m<-as.matrix(data_m)
  
  mzvec<-data_matrix[,1]
  rtvec<-data_matrix[,2]
  
  rnames<-paste(mzvec,rtvec,sep="_")
  
  
  if(is.na(name[1])==TRUE){
    
    #mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
    rownames(data_m)<-as.character(rnames)
  }else{
    
    rownames(data_m)<-name
    
    rnames<-name
  }
  
  #  rownames(data_m)<-as.character(rnames)
  
  
  classlabels_orig<-classlabels
  
  if(analysistype=="onewayrepeat" | analysistype=="twowayrepeat" | analysistype=="onewayanovarepeat" | analysistype=="twowayanovarepeat"){
    
    pairedanalysis=TRUE
  }
  
  if(pairedanalysis==TRUE){
    
    # paireddesign=classlabels_orig[,2]
    
    #classlabels_orig<-classlabels_orig[,-c(2)]
  }
  
  if(dim(classlabels_orig)[2]>2){
    
    #classlabels_orig[,2]<-as.factor(paste("A",as.character(classlabels_orig[,2]),sep=""))
    #classlabels_orig[,3]<-as.factor(paste("B",as.character(classlabels_orig[,3]),sep=""))
    # print(head(classlabels_orig))
    
    if(analysistype=="twoway" | analysistype=="twowayrepeat" | analysistype=="twowayanova"){
      
      if(alphabetical.order==FALSE){
        
        classlabels_orig[,2]<-factor(classlabels_orig[,2],levels=unique(classlabels_orig[,2]))
        
        classlabels_orig[,3]<-factor(classlabels_orig[,3],levels=unique(classlabels_orig[,3]))
      }
      
      classgroup<-paste(classlabels_orig[,2],":",classlabels_orig[,3],sep="") #classlabels_orig[,2]:classlabels_orig[,3]
      do_pca_anova=FALSE
    }else{
      
      if(alphabetical.order==FALSE){
        
        classlabels_orig[,2]<-factor(classlabels_orig[,2],levels=unique(classlabels_orig[,2]))
      }
      
      
      classgroup<-classlabels_orig[,2]
      do_pca_anova=TRUE
    }
    
    col_class_levels<-levels(as.factor(classlabels_orig[,2]))
  }else{
    
    if(alphabetical.order==FALSE){
      
      classlabels_orig[,2]<-factor(classlabels_orig[,2],levels=unique(classlabels_orig[,2]))
    }
    
    classgroup<-classlabels_orig[,2]
    col_class_levels<-levels(as.factor(classlabels_orig[,2]))
    
    do_pca_anova=TRUE
  }
  
  
  
  
  
  
  
  
  
  class_labels_levels<-levels(as.factor(classgroup))
  ordered_labels<-classgroup
  
  class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
  
  if(newdevice==TRUE){
    
    fname<-paste("timeseriesplots",filename,".pdf",sep="")
    pdf(fname)
  }
  
  if(is.na(sample.col.opt)==FALSE)
  {
    if(sample.col.opt=="default"){
      
      col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
                 "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
                 "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
                 "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
      
    }else{
      if(sample.col.opt=="topo"){
        #col_vec<-topo.colors(256) #length(class_labels_levels))
        
        #col_vec<-col_vec[seq(1,length(col_vec),)]
        
        col_vec <- topo.colors(length(col_class_levels), alpha=alphacol)
      }else{
        if(sample.col.opt=="heat"){
          #col_vec<-heat.colors(256) #length(class_labels_levels))
          
          col_vec <- heat.colors(length(col_class_levels), alpha=alphacol)
        }else{
          if(sample.col.opt=="rainbow"){
            #col_vec<-heat.colors(256) #length(class_labels_levels))
            col_vec<-rainbow(length(col_class_levels), start = 0, end = alphacol)
            
            #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
          }else{
            
            if(sample.col.opt=="terrain"){
              #col_vec<-heat.colors(256) #length(class_labels_levels))
              #col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
              
              col_vec <- cm.colors(length(col_class_levels), alpha=alphacol)
            }else{
              if(is.na(sample.col.opt)==TRUE){
                col_vec<-c("black")
              }else{
                
                if(sample.col.opt=="colorblind"){
                  #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                  # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                  
                  if(length(col_class_levels)<9){
                    
                    col_vec <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "#E64B35FF", "grey57")
                    
                  }else{
                    
                    #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                               "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                    
                  }
                  
                  
                }else{
                  
                  check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                  
                  if(length(check_brewer)>0){
                    
                    sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                    
                    col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(col_class_levels))
                    
                  }else{
                    
                    #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                    if(sample.col.opt=="journal"){
                      
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
                      #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                      
                    #  if(length(sample.col.opt)==1){
                       # col_vec <-rep(sample.col.opt,length(col_class_levels))
                      #}else{
                        
                       # colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(col_class_levels))
                        
                      #}
                      
                      
                      if(length(sample.col.opt)==1){
                        col_vec <-rep(sample.col.opt,length(class_labels_levels))
                      }else{
                        
                        #if(length(sample.col.opt)>=length(class_labels_levels)){
                          
                          col_vec <-sample.col.opt
                          col_vec <- rep(col_vec,length(class_labels_levels))
                          
                          
                        #}else{
                         # colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                        #}
                        
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
  }else{
    
    col_vec<-c("black")
    
  }
  
  #    print(class_labels_levels)
  
  ordered_labels={}
  num_samps_group<-new("list")
  num_samps_group[[1]]<-0
  groupwiseindex<-new("list")
  groupwiseindex[[1]]<-0
  
  
  col_vec=alpha(colour=col_vec,alpha=alphacol)
  
  class_col_vec=col_vec
  S<-new("list")
  for(c in 1:length(class_labels_levels))
  {
    
    classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
    
    num_samps_group[[c]]<-length(classlabels_index)
    groupwiseindex[[c]]<-classlabels_index
    
    
    
    
  }
  
  
  
  sampleclass<-{}
  patientcolors<-{}
  
  classlabels<-as.data.frame(classlabels)
  
  
  for(c in 1:length(class_labels_levels)){
    
    num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
    
    
    sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
    
    patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group[[c]]))
  }
  
  pca.cex.val=point.cex.val
  
  if(length(mzvec)>4){
    max_per_row<-3
    
    
    par_rows<-ceiling(9/max_per_row)
    
  }else{
    max_per_row<-length(mzvec)
    par_rows<-1
  }
  
  # The palette with black:
  #col_vec <- c("#E69F00","#000000","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # To use for fills, add
  # scale_fill_manual(values=col_vec)
  
  # To use for line and point colors, add
  #scale_colour_manual(values=col_vec)
  
  file_ind<-0
  boxplots_fname<-paste("xyplots.pdf",sep="")
  
  if(dim(data_m)[1]>=1){
    
    
    t1<-table(classgroup)
    l1<-levels(classgroup)
    
    l1<-levels(as.factor(classgroup))
    
    
    patientcolors <- rep(col_vec[1:length(t1)], t1)
    
    
    get_pooled_sp<-function(n1,n2,S1,S2){a<-(n1-1)*S1;b<-(n2-1)*S2;c<-(n1+n2-2);return((a+b)/c)}
    
    
    class_levels<-levels(as.factor(classlabels_orig[,2]))
    
    
    df_matrix<-{}
    
    class_col_vec=col_vec
    
    if(dim(classlabels_orig)[2]>2){
      class_levels<-levels(as.factor(classlabels_orig[,2]))
      t1<-table(classlabels_orig[,3])
      
      
      
      
  #  save(classlabels_orig,classgroup,data_m,file="myData1.Rda")
      
      
      myData<-cbind(classgroup,as.data.frame(classlabels_orig[,2]),as.data.frame(classlabels_orig[,3]),t(data_m))
      
      myData<-as.data.frame(myData)
      
      #
      if(alphabetical.order==FALSE){
        myData[,2] <- factor(myData[,2], levels=unique(myData[,2]))
        myData[,3] <- factor(myData[,3], levels=unique(myData[,3]))
        t1=table(myData[,3])
      }
      
      
     # save(myData,file="myData_A.Rda")
     if(FALSE){ 
      c1<-WGCNA::cor(t(myData[,-c(1:3)]),use = 'pairwise.complete.obs')
      diag(c1)<-NA
      c2=cbind(myData[,c(1:3)],c1)
      
      myData_sum_cor<-{}
      levels_groups<-levels(as.factor(c2[,1]))
      
      myData_sum_cor<-lapply(levels_groups,function(g){
        
        x<-c1[which(c2[,1]==g),which(c2[,1]==g)]
        
        x<-as.numeric(as.character(x))
        se_val<-round(sd(x,na.rm=TRUE)/sqrt(length(x)),2)
        mean_val<-round(mean(x,na.rm=TRUE),2)
        upper.conf.limit=mean_val+1.96*se_val
        lower.conf.limit=mean_val-1.96*se_val
        
        y=as.data.frame(c(group=g,mean.Pearson.cor = mean_val, 
                          sd.Pearson.cor = round(sd(x,na.rm=TRUE),2),n = length(which(c2[,1]==g)),se.Pearson.cor=se_val,
          lower.conf.limit=lower.conf.limit, 
          upper.conf.limit=upper.conf.limit,
          width.conf.interval=round(upper.conf.limit-lower.conf.limit,2)))
        colnames(y)<-NULL
        return(t(y))
        
      })
      
     myData_sum_cor<-ldply(myData_sum_cor,rbind)    
     myData_sum_cor<-as.data.frame(myData_sum_cor)
     
     #write.table(myData_sum_cor,file="myData_sum_cor.txt",sep="\t")
     myData_sum_cor$mean.Pearson.cor<-as.numeric(as.character(myData_sum_cor$mean.Pearson.cor))
     myData_sum_cor$lower.conf.limit<-as.numeric(as.character(myData_sum_cor$lower.conf.limit))
     myData_sum_cor$upper.conf.limit<-as.numeric(as.character(myData_sum_cor$upper.conf.limit))
     
   #  save(myData_sum_cor,sizeval,cex.val,class_col_vec,file="mysamp.Rda")
     group_num<-gsub(myData_sum_cor$group,pattern=":[a-z|A-Z|0-9|_|-]*",replacement="")
     myData_sum_cor<-cbind(myData_sum_cor,group_num)
     Factor2<-gsub(myData_sum_cor$group,pattern="[a-z|A-Z|0-9|-|_]*:",replacement="")
     myData_sum_cor<-cbind(myData_sum_cor,Factor2)
     
     plot_sample_cor<-ggplot(myData_sum_cor,aes(y=mean.Pearson.cor,x=Factor2,color=group_num2))+
       geom_point()+facet_wrap(~group_num,nrow=1)+geom_errorbar(aes(ymin=lower.conf.limit,ymax=upper.conf.limit))+ theme_bw() + theme(panel.border = element_blank(), 
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       axis.line = element_line(colour = "black",size=sizeval),
                                       axis.text= element_text(size=14*cex.val), 
                                       axis.title=element_text(size=16*cex.val,face="bold"),
                                       strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
strip.text = element_text(face="bold"))+scale_y_continuous(name="mean Pearson correlation\nbetween biological replicates in each group",breaks = scales::pretty_breaks(n = 10))
     #plot_sample_cor<-plot_sample_cor+scale_color_manual(values=unique(class_col_vec))
     #plot_sample_cor<-plot_sample_cor+scale_x_discrete(name ="Class")
    plot_sample_cor<-plot_sample_cor+theme(axis.text.x = element_text(angle = 45, hjust = 1))
      #guides(fill = guide_legend(title = FALSE),colour=FALSE)+
    #  guides(colour = guide_legend(override.aes = list(shape = unique(shape_vec1))))+
      
   # print(plot_sample_cor)
    }
      
      
      
      myData_sum <- do.call(data.frame,aggregate(list(myData[,-c(1:3)]),
                                                 by = list(myData[,2],myData[,3]),
                                                 FUN = function(x){
                                                   
                                                   x<-as.numeric(as.character(x))
                                                   
                                                   c(mean = mean(x), sd = sd(x),n = length(x),se=sd(x)/sqrt(length(x)))
                                                   
                                                 }))
      
      
      
      label_inc_list<-seq(3,dim(myData_sum)[2],4)
      
   #   save(myData_sum,file="myDatasum.Rda")
      
      pch_vec<-c(19,17,23,22,13,0:12)
      #save(label_inc_list,file="label_inc_list.Rda")
      #save(rnames,file="rnames.Rda")
      #save(data_m,file="data_m.Rda")
      
      if(output.format=="pdf"){
        
        #par(mfrow=c(2,4),family="sans",cex=cex.val)
      }
      

      plot_res<-lapply(seq(3,dim(myData_sum)[2],4),function(pc)
      {
        
        df.summary<-myData_sum[,c(pc:(pc+3),1,2)] #cbind(xvec,classgroup,classlabels_orig[,2],classlabels_orig[,3])
        
        get_label_ind<-which(label_inc_list==pc)
        if(pairedanalysis==FALSE){
          
          
          if(error.bar==FALSE){
            #mzname<-paste(rnames[get_label_ind]," distribution in each group using ",filename," feats vs factors",sep="")
            
            
          }else{
            #mzname<-paste(rnames[get_label_ind]," distribution with 95% confidence interval in each group using ",filename," feats vs factors",sep="")
          }
        }else{
          
          if(error.bar==FALSE){
            mzname<-paste(rnames[get_label_ind]," distribution in each group using ",filename," feats vs time",sep="")
          }else{
            
            mzname<-paste(rnames[get_label_ind]," distribution 95% confidence interval in each group using ",filename," feats vs time",sep="")
          }
          
        }
        
        mzname<-rnames[get_label_ind]
        
        colnames(df.summary)<-c("Intensity","sd","number","se","Class","time")
        ymax = df.summary$Intensity + 1.96*df.summary$se
        
        ymin = df.summary$Intensity - 1.96*df.summary$se
        
        df_write_temp<-cbind(mzname,df.summary[,c(5,6,3,1,2,4)],ymin,ymax)
        colnames(df_write_temp)<-c("Name","Class","time","Number of subjects","mean","Std.deviation","Std.error","lower.limit.95%CI","upper.limit.95%CI")
        
       #write.table(df_write_temp,file=) 
        
        
        df.summary<-as.data.frame(df.summary)
        df.summary$Class<-as.numeric(as.factor(df.summary$Class))
        
        df_time<-levels(df.summary$time)
        #df.summary$time<-as.numeric(as.factor(df.summary$time))
        
       
        #   max_yval<-ceiling(max(df.summary$Intensity+(2.5* df.summary$se),na.rm=TRUE)) #max(ymax,na.rm=TRUE)
        
        #min_yval<-max(0,floor(min(df.summary$Intensity-(2.5* df.summary$se),na.rm=TRUE)))
        
        max_yval<-ceiling(max(df.summary$Intensity+(2.5* df.summary$se),na.rm=TRUE)) #round(max( df.summary$Intensity+(4* df.summary$se),na.rm=TRUE))
        
        
        min_yval<-floor(min(df.summary$Intensity-(2.5* df.summary$se),na.rm=TRUE))
        
        
        
        
        
        if(alphabetical.order==FALSE){
          classlabels_orig[,2] <- factor(classlabels_orig[,2], levels=unique(classlabels_orig[,2]))
          classlabels_orig[,3] <- factor(classlabels_orig[,3], levels=unique(classlabels_orig[,3]))
          class_levels_time<-levels(as.factor(classlabels_orig[,3]))
          class_levels<-levels(as.factor(classlabels_orig[,2]))
        }else{
          
          class_levels_time<-levels(as.factor(classlabels_orig[,3]))
        }
        
        Class<-{}
        for(cnum in 1:length(class_levels)){
          
          Class<-c(Class,rep(class_levels[cnum],length(t1)))
          
          df.summary$Class[which(df.summary$Class==cnum)]<-class_levels[cnum]
          
        }
        
        if(alphabetical.order==FALSE){
        df.summary$Class<-factor(df.summary$Class,levels=unique(df.summary$Class))
        df.summary$time<-factor(df.summary$time,levels=unique(df.summary$time))
        
        }
        Class<-unique(Class)
        
        t1<-table(df.summary$Class)
        
        
        
        #print(df.summary)
        
        df.summary$x<-factor(df.summary$time) # time.hour
       # save(df.summary,file="df.summary.Rda")
        #savecol_vec,file="col_vec.Rda")
        #saveclass_levels,file="class_levels.Rda")
        #savelist=ls(),file="debuglinep.Rda")
        
        shape_vec=c(0:2,4:25)
        t1shape=levels(factor(df.summary$Class))
        
        shape_vec1<-(df.summary$Class)
        
        for(i in 1:length(t1shape)){
          shape_vec1[which(shape_vec1==t1shape[i])]<-shape_vec[i]
        }
        
        #  shape_vec1<-lapply(1:length(t1shape),function(i){
        
        #   return(rep(shape_vec[i],t1shape[i]))
        #})
        
        #shape_vec1<-unlist(shape_vec1)
        
        
        #col_vec<-colfunc(length(t1))
        #ylim(0,ymax) +
        sizeval=1.5
        
        unique_class_col_vec<-unique(class_col_vec)
        
        shape_vec1<-as.numeric(as.character(shape_vec1))
        
        
         fname1<-paste("lp",pc,".Rda",sep="")
      #  save(shape_vec1,df.summary,class_labels_levels,unique_class_col_vec,
       #      cex.val,sizeval,class_col_vec,ymin,ymax,pairedanalysis,ylabel,
        #     df_time,file=fname1) #"lp1.Rda")
        
        
        if(pairedanalysis==TRUE){
          
          options(repr.plot.width = 3, repr.plot.height = 2)
          #ylim(0,max_yval) +
          
          plot_res<-suppressWarnings(ggplot(df.summary, 
                                            aes(x = as.factor(x), y = Intensity,color = Class,linetype=Class)) +
                                       geom_point(size = cex.val,shape=shape_vec1) + 
                                       geom_line(aes(group =Class),size=sizeval)  + 
                                       geom_errorbar(aes(ymin = ymin, ymax = ymax),size=0.3,width=0.1) + 
                                       ylab(ylabel) + theme_bw() + theme(panel.border = element_blank(), 
                                                                         panel.grid.major = element_blank(),
                                                                         panel.grid.minor = element_blank(),
                                                                         axis.line = element_line(colour = "black",size=sizeval),
                                                                         axis.text= element_text(size=14*cex.val), axis.title=element_text(size=16*cex.val,face="bold"),
                                                                         strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                         strip.text = element_text(face="bold")) + 
                                       scale_y_continuous(name=ylabel,breaks = scales::pretty_breaks(n = 10)))
          
                                                                                                                                                                                                                                                                                                                 
          if(length(unique_class_col_vec)==1){
            plot_res<-plot_res+scale_color_manual(values=rep(unique(class_col_vec),length(class_labels_levels)))
            
          }else{
            plot_res<-plot_res+scale_color_manual(values=unique(class_col_vec))
          }
          
          plot_res<-plot_res+scale_x_discrete(name ="Time",limits=unique(df_time))
          
          plot_res<-plot_res+guides(fill = guide_legend(title = FALSE),colour=FALSE)+guides(colour = guide_legend(override.aes = list(shape = unique(shape_vec1))))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
          
          
          
          #+ scale_colour_manual(values=col_vec)
        }else{
          
          if(timeseries.lineplots==TRUE){
            
            
            #savelist=ls(),file="r.Rda")
            #print("here this ABC")
            #ylim(0,max_yval) +
            
            
        #    save(shape_vec1,df.summary,class_labels_levels,sizeval,ymin, ymax,ylabel,cex.val,
         #        unique_class_col_vec,class_col_vec,ymin,ymax,pairedanalysis,
          #      timeseries.lineplots,df_time,file="lp2.Rda")
            
            
            #if(FALSE)
            {
              shape_vec1<-as.numeric(as.character(shape_vec1))
              options(repr.plot.width = 3, repr.plot.height = 2)
              
              if(is.na(lineplot.lty.option)==TRUE){
              plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(x), y = Intensity,colour = Class,
                                                                linetype=Class)) + 
                                           geom_point(size = 4,shape=shape_vec1) + 
                                           geom_line(aes(group =Class),size=sizeval) + 
                                           geom_errorbar(aes(ymin = ymin, ymax = ymax),size=0.3,width=0.1)+ ylab(ylabel) + theme_bw() +
                                           theme(panel.border = element_blank(),
                                                 panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                 axis.text= element_text(size=14*cex.val), axis.title=element_text(size=16,face="bold"),
                                                 strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                 strip.text = element_text(face="bold")) + 
                                           scale_y_continuous(name=ylabel,
                                                              breaks = scales::pretty_breaks(n = 10)))
                                                                                                                                                                                                                                                                                                               
              }else{
                #linetype=lineplot.lty.option
                plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(x), y = Intensity,colour = Class
                                                                  )) + 
                                             geom_point(size = 4,shape=shape_vec1) + 
                                             geom_line(aes(group =Class),size=sizeval) + 
                                             geom_errorbar(aes(ymin = ymin, ymax = ymax),size=0.3,width=0.1)+ ylab(ylabel) + theme_bw() +
                                             theme(panel.border = element_blank(),
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                   axis.text= element_text(size=14*cex.val), axis.title=element_text(size=16,face="bold"),
                                                   strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                   strip.text = element_text(face="bold")) + 
                                             scale_y_continuous(name=ylabel,
                                                                breaks = scales::pretty_breaks(n = 10)))
                
                
                
              }
              if(length(unique_class_col_vec)==1){
                plot_res<-plot_res+scale_color_manual(values=rep(unique(class_col_vec),length(class_labels_levels)))
                
              }else{
                plot_res<-plot_res+scale_color_manual(values=unique(class_col_vec))
              }
              #save(plot_res,file="p1.Rda")
              plot_res<-plot_res+scale_x_discrete(name ="Time",limits=unique(df_time))
              
              #save(plot_res,file="p2.Rda")
              plot_res<-plot_res+guides(fill = guide_legend(title = FALSE),colour=FALSE)+guides(colour = guide_legend(override.aes = list(shape = unique(shape_vec1))))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
            }
            
            
        
            
            
            
            
            
            
            
          }else{
            #ylim(0,max_yval)
            options(repr.plot.width = 3, repr.plot.height = 2)
            
            plot_res<-ggplot(df.summary, aes(x = as.factor(x), y = Intensity,color = Class)) + geom_point(size = pca.cex.val,shape=shape_vec1) + geom_errorbar(aes(ymin = ymin, ymax = ymax),size=sizeval) + scale_color_manual(values=col_vec) + ylab(ylabel) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                                                                                                                                                                                                                                         axis.text= element_text(size=14*cex.val), axis.title=element_text(size=16*cex.val,face="bold"),
                                                                                                                                                                                                                                                                         strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                                                                                                                                                         strip.text = element_text(face="bold")) + scale_y_continuous(name=ylabel,breaks = scales::pretty_breaks(n = 10))
            
            plot_res<-plot_res+scale_x_discrete(name ="Factor",limits=unique(df_time))
            
          }
          #+ scale_colour_manual(values=col_vec)
        }
        
        plot_res<-plot_res + ggtitle(mzname)
        
        plot_res<-plot_res + theme(plot.title = element_text(size=8))
        # x axis title
        plot_res<-plot_res + theme(axis.title.x = element_text(size=14*cex.val))
        # y axis title
        plot_res<-plot_res + theme(axis.title.y = element_text(size=14*cex.val))+theme(axis.text.x = element_text(angle = 45, hjust = 1,size=11*cex.val))
        
        
        # print(plot_res)
        
        
        
        plot_res<-plot_res+ggtitle(mzname)
        
        #return(plot_res)
        #print(plot_res + ggtitle(mzname))
        return(list("df_write_temp"=df_write_temp,"plot_res"=plot_res))
      })
      
   #   save(plot_res,file="plot_res.Rda")
      
      
      
      if(multiple.figures.perpanel==TRUE){
        
          res<-lapply(seq(1,length(plot_res),4),function(i){
            p1=plot_res[[i]]$plot_res
            p2={}
            p3={}
            p4={}
            if((i+1)<length(plot_res)){
              p2=plot_res[[i+1]]$plot_res
            }
            if((i+2)<length(plot_res)){
              p3=plot_res[[i+2]]$plot_res
            }
            if((i+3)<length(plot_res)){
              p4=plot_res[[i+3]]$plot_res
            }
            if(length(p2)>0 && length(p3)>0 && length(p4)>0){
              
                figure<-ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2,heights=c(4,4),width=c(6,6),legend=FALSE,align = c("hv"))
            }else{
              if(length(p2)>0 && length(p3)>0){
                
                figure<-ggarrange(p1,p2,p3,ncol = 2, nrow = 2,heights=c(4,4),width=c(6,6),legend=FALSE,align = c("hv"))
              }else{
                
                if(length(p2)>0){
                  
                  figure<-ggarrange(p1,p2,ncol = 2, nrow = 1,heights=c(4,4),width=c(6,6),legend=FALSE,align = c("hv"))
                }else{
                  
                  
                    figure<-ggarrange(p1,ncol = 1, nrow = 1,heights=c(4,4),width=c(6,6),legend=FALSE,align = c("hv"))
                  
                }
              }
              
            }
            
            
                   return(figure)
          })
      
      for(i in 1:length(res)){
        
        print(res[[i]][[1]])
      }
      
      legend <- ggpubr::get_legend(plot_res[[i]]$plot_res)
      
      suppressMessages(library(grid))
      
      grid.newpage()
      grid.draw(legend)
      }else{
      
        res<-lapply(seq(1,length(plot_res),1),function(i){
          p1=plot_res[[i]]$plot_res
         
          figure<-ggarrange(p1,ncol = 1, nrow = 1,heights=plot.height,width=plot.width,legend=TRUE,align = c("hv"))
          return(figure)
        })
        
        for(i in 1:length(res)){
          
          print(res[[i]][[1]])
        }
      
        legend <- ggpubr::get_legend(plot_res[[i]]$plot_res)
        
        suppressMessages(library(grid))
        
        grid.newpage()
        grid.draw(legend)
    }
      
     
    }else{
      
      
      #one-way anova or one-way with repeated measures
      class_levels<-levels(as.factor(classgroup))
      t1<-table(classgroup)
      
      
      
      myData<-cbind(classgroup,classlabels_orig[,2],t(data_m))
      
      myData<-as.data.frame(myData)
      
      if(alphabetical.order==FALSE){
        myData[,2] <- factor(myData[,2], levels=unique(myData[,2]))
        
        
      }
      
      myData_sum <- do.call(data.frame,aggregate(list(myData[,-c(1:2)]),
                                                 by = list(myData[,2]),
                                                 FUN = function(x){
                                                   
                                                   x<-as.numeric(as.character(x))
                                                   c(mean = mean(x), sd = sd(x),n = length(x),se=sd(x)/sqrt(length(x)))
                                                   
                                                 }))
      
      
      label_inc_list<-seq(2,dim(myData_sum)[2],4)
      
      ###savelabel_inc_list,file="label_inc_list.Rda")
      ###savernames,file="rnames.Rda")
      
      #savemyData_sum,classgroup,classlabels_orig,label_inc_list,rnames,file="myDatasum.Rda")
      plot_res<-lapply(seq(2,dim(myData_sum)[2],4),function(pc)
        #for(pc in seq(2,dim(myData_sum)[2],4))
      {
        
        #print("here 0705")
        #  print(pc)
        df.summary<-myData_sum[,c(pc:(pc+3),1)] #cbind(xvec,classgroup,classlabels_orig[,2],classlabels_orig[,3])
        
        get_label_ind<-which(label_inc_list==pc)
        if(pairedanalysis==FALSE){
          
          
          if(error.bar==FALSE){
            mzname<-paste(rnames[get_label_ind]," distribution \nin each group using ",filename," feats vs factors",sep="")
          }else{
            mzname<-paste(rnames[get_label_ind]," distribution \nwith 95% confidence interval in each group using ",filename," feats vs factors",sep="")
          }
        }else{
          
          if(error.bar==FALSE){
            mzname<-paste(rnames[get_label_ind]," distribution \nin each group using ",filename," feats vs time",sep="")
          }else{
            
            mzname<-paste(rnames[get_label_ind]," distribution \n95% confidence interval in each group using ",filename," feats vs time",sep="")
          }
          
        }
        
        colnames(df.summary)<-c("Intensity","sd","number","se","Class")
        df.summary<-as.data.frame(df.summary)
        
        df.summary$Class<-as.numeric(as.factor(df.summary$Class))
        
        
        ymax = df.summary$Intensity + 1.96*df.summary$se
        
        ymin = df.summary$Intensity - 1.96*df.summary$se
        
        df_write_temp<-cbind(mzname,df.summary[,c(5,3,1,2,4)],ymin,ymax)
        
        colnames(df_write_temp)<-c("Name","Class","Number of subjects","mean","Std.deviation","Std.error","lower.limit.95%CI","upper.limit.95%CI")
        
        
        max_yval<-ceiling(max(df.summary$Intensity+(2.5* df.summary$se),na.rm=TRUE)) #round(max( df.summary$Intensity+(4* df.summary$se),na.rm=TRUE))
        
        
        min_yval<-floor(min(df.summary$Intensity-(2.5* df.summary$se),na.rm=TRUE))
        
        
        Class<-{}
        for(cnum in 1:length(class_levels)){
          
          Class<-c(Class,rep(class_levels[cnum],length(t1)))
          
          df.summary$Class[which(df.summary$Class==cnum)]<-class_levels[cnum]
          
        }
        
        
        Class<-unique(Class)
        
        t1<-table(df.summary$Class)
        
        
        #print(df.summary)
        #savedf.summary,file="df.summarydotp1.Rda")
        # ##save(df.summary,class_levels,col_vec,ymin,ymax,pairedanalysis,ylabel,file="lp1.Rda")
        
        
        ##savelist=ls(),file="debuglinep1.Rda")
        time.hour<- c(unique(as.character(classlabels_orig[,2])),unique(as.character(classlabels_orig[,2])))
        #Score<-df.summary$Score
        #df.summary$x<-time.hour
        
        #print(Class)
        # scale_color_hue(l=40)
        sizeval=1.5
        if(pairedanalysis==FALSE){
          
          if(timeseries.lineplots==FALSE){
            #ylim(0,max(ymax,na.rm=TRUE)) +
            options(repr.plot.width = 3, repr.plot.height = 2)
            plot_res<-ggplot(df.summary, aes(x = Class, y = Intensity,color = Class)) + geom_point(size=pca.cex.val)  + geom_errorbar(aes(ymin = ymin, ymax = ymax),size=0.3,width=0.1) + ylab(ylabel) + scale_x_continuous(breaks=seq(1,length(class_levels))) + scale_color_manual(values=col_vec) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                                                                                                                                                                                                                                                                                          axis.text= element_text(size=14*cex.val), axis.title=element_text(size=16*cex.val,face="bold"),
                                                                                                                                                                                                                                                                                                                          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                                                                                                                                                                                                          strip.text = element_text(face="bold")) + scale_y_continuous(name=ylabel,breaks = scales::pretty_breaks(n = 10))
            
            
            #plot_res<-plot_res+scale_x_discrete(name ="Factor",limits=unique(df.summary$time))
          }else{
            
            #ylim(0,max(ymax,na.rm=TRUE)) +
            options(repr.plot.width = 3, repr.plot.height = 2)
            plot_res<-ggplot(df.summary, aes(x = Class, y = Intensity,group=1)) + geom_point(size=pca.cex.val)  + geom_line(size=sizeval) + geom_errorbar(aes(ymin = ymin, ymax = ymax),size=0.3,width=0.1) + ylab(ylabel) +  scale_x_continuous(breaks=seq(1,length(class_levels))) + scale_color_manual(values=col_vec) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                                                                                                                                                                                                                                                                                                               axis.text= element_text(size=14*cex.val), axis.title=element_text(size=16*cex.val,face="bold"),
                                                                                                                                                                                                                                                                                                                                               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                                                                                                                                                                                                                               strip.text = element_text(face="bold")) + scale_y_continuous(name=ylabel,breaks = scales::pretty_breaks(n = 10))
            
          }
          
          
          #+ scale_colour_manual(values=col_vec)
        }else{
          
          
          #plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Class)) + geom_point(size = 2) + geom_line()  #+ geom_errorbar(aes(ymin = ymin, ymax = ymax))
          
          #ylim(0,max(ymax,na.rm=TRUE)) +
          options(repr.plot.width = 3, repr.plot.height = 2)
          plot_res<-ggplot(df.summary, aes(x = Class, y = Intensity,group=1)) +  geom_point(size = pca.cex.val) + geom_line(size=sizeval)  + geom_errorbar(aes(ymin = ymin, ymax = ymax),width=0.2,size=sizeval) + ylab(ylabel) +  scale_color_manual(values=col_vec) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),axis.text= element_text(size=14*cex.val), axis.title=element_text(size=16*cex.val,face="bold"),
                                                                                                                                                                                                                                                                                           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                                                                                                                                                                           strip.text = element_text(face="bold")) + scale_y_continuous(name=ylabel,breaks = scales::pretty_breaks(n = 10))
          
          #+ scale_colour_manual(values=col_vec)
          
          #plot_res<-ggplot(df.summary, aes(x = x, y = Intensity,color = Class)) +  geom_point(size = pca.cex.val,colour=col_vec) + geom_line(aes(group =Class),colour=col_vec)  + geom_errorbar(aes(ymin = ymin, ymax = ymax),colour=col_vec) + scale_color_hue(l=40)
          
          #plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Class)) +  geom_point() + geom_line(aes(group =Class)) # + geom_line()
          
        }
        
        
        plot_res<-plot_res + ggtitle(mzname)
        
        plot_res<-plot_res + theme(plot.title = element_text(size=8))
        # x axis title
        plot_res<-plot_res + theme(axis.title.x = element_text(size=14*cex.val))
        # y axis title
        plot_res<-plot_res + theme(axis.title.y = element_text(size=14*cex.val))+theme(axis.text.x = element_text(angle = 45, hjust = 1,size=11*cex.val))
        
        #  print(plot_res) #
        
        
      #  return(plot_res)
       
        return(list("plot_res"=plot_res,"df_write_temp"=df_write_temp)) 
        
        
      }) #return plot
      ##save(plot_res,file="plot_res.Rda")
      
      #pdf("t.pdf",onefile=TRUE)
      
      
      if(multiple.figures.perpanel==TRUE){
      
      res<-lapply(seq(1,length(plot_res),4),function(i){
        #for(i in seq(1,length(plot_res),4)){
        p1=plot_res[[i]]$plot_res
        p2=""
        p3=""
        p4=""
        
        if((i+1)<length(plot_res)){
          p2=plot_res[[i+1]]$plot_res
        }
        if((i+2)<length(plot_res)){
          p3=plot_res[[i+2]]$plot_res
        }
        if((i+3)<length(plot_res)){
          
          p4=plot_res[[i+3]]$plot_res
        }
        
        figure<-ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2,heights=c(4,4),width=c(6,6),
                          legend=FALSE,align = c("hv"))
        return(figure)
      })
      
      
      for(i in 1:length(res)){
        
        print(res[[i]][[1]])
      }
      
      legend <- ggpubr::get_legend(plot_res[[i]]$plot_res)
      
      suppressMessages(library(grid))
      
      grid.newpage()
      grid.draw(legend)
      
      }else{
        
        
        res<-lapply(seq(1,length(plot_res),1),function(i){
          #for(i in seq(1,length(plot_res),4)){
          p1=plot_res[[i]]$plot_res
          p2=""
          p3=""
          p4=""
          
        
        
          
          figure<-ggarrange(p1,ncol = 1, nrow = 1,heights=c(4,4),width=c(6,6),
                            legend=FALSE,align = c("hv"))
          return(figure)
        })
        
        
        for(i in 1:length(res)){
          
          print(res[[i]][[1]])
        }
        
        legend <- ggpubr::get_legend(plot_res[[i]]$plot_res)
        
        suppressMessages(library(grid))
        
        grid.newpage()
        grid.draw(legend)
        
      }
      
    }
  }else{
    
    
    
  }
  #df_fname<-paste("PC_score_distribution_matrix_",filename,"features.txt",sep="")
  #write.table(df_matrix,file=df_fname,sep="\t",row.names=TRUE)
  if(newdevice==TRUE){
    
    try(dev.off(),silent=TRUE)
  }
  try(dev.off(),silent=TRUE)
  return(plot_res)
  
}

get_data_summary <-
function(X=NA,Y=NA,feature_table_file=NA,parentoutput_dir=NA,
                               class_labels_file=NA,alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,
                               point.cex.val=4,legendlocation="topright",pca.ellipse=TRUE,
                               ellipse.conf.level=0.95,filename="all",newdevice=FALSE,
                               lineplot.col.opt=c("grey57"),ylabel="Intensity",error.bar=TRUE,
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
        classlabels_orig[,2] <- factor(classlabels_orig[,2], levels=unique(classlabels_orig[,2]))
        classlabels_orig[,3] <- factor(classlabels_orig[,3], levels=unique(classlabels_orig[,3]))
        #t1=table(classlabels_orig[,3])
      }
      
      
      classgroup<-paste(classlabels_orig[,2],":",classlabels_orig[,3],sep="") #classlabels_orig[,2]:classlabels_orig[,3]
      do_pca_anova=FALSE
    }else{
      
      if(alphabetical.order==FALSE){
        classlabels_orig[,2] <- factor(classlabels_orig[,2], levels=unique(classlabels_orig[,2]))
       
      }
      
      classgroup<-classlabels_orig[,2]
      do_pca_anova=TRUE
    }
    
    col_class_levels<-levels(as.factor(classlabels_orig[,2]))
  }else{
    
    if(alphabetical.order==FALSE){
      classlabels_orig[,2] <- factor(classlabels_orig[,2], levels=unique(classlabels_orig[,2]))
     
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
                      
                     # if(length(sample.col.opt)==1){
                      #  col_vec <-rep(sample.col.opt,length(col_class_levels))
                      #}else{
                        
                       # colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(col_class_levels))
                        
                        ##savecolfunc,file="colfunc.Rda")
                        
                      #}
                      
                      if(length(sample.col.opt)==1){
                        col_vec <-rep(sample.col.opt,length(class_labels_levels))
                      }else{
                        
                        if(length(sample.col.opt)>=length(class_labels_levels)){
                          
                          col_vec <-sample.col.opt
                          col_vec <- rep(col_vec,length(class_labels_levels))
                          
                          
                        }else{
                          colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
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
      
      
      myData<-cbind(as.data.frame(classgroup),as.data.frame(classlabels_orig[,2]),as.data.frame(classlabels_orig[,3]),t(data_m))
      
      myData<-as.data.frame(myData)
      
      #
      if(alphabetical.order==FALSE){
        myData[,2] <- factor(myData[,2], levels=unique(myData[,2]))
        myData[,3] <- factor(myData[,3], levels=unique(myData[,3]))
        t1=table(myData[,3])
      }
      
     
      
      
      
      
      myData_sum <- do.call(data.frame,aggregate(list(myData[,-c(1:3)]),
                                                 by = list(myData[,2],myData[,3]),
                                                 FUN = function(x){
                                                   
                                                   x<-as.numeric(as.character(x))
                                                   
                                                   c(mean = mean(x), sd = sd(x),n = length(x),se=sd(x)/sqrt(length(x)))
                                                   
                                                 }))
      
      
      
      label_inc_list<-seq(3,dim(myData_sum)[2],4)
     
      pch_vec<-c(19,17,23,22,13,0:12)
     
      if(output.format=="pdf"){
        
        par(mfrow=c(2,4),family="sans",cex=cex.val)
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
        
        return(df_write_temp)    
      })
    }
    else
    {
      
      
      #one-way anova or one-way with repeated measures
      class_levels<-levels(as.factor(classgroup))
      t1<-table(classgroup)
      
      
      myData<-cbind(as.data.frame(classgroup),as.data.frame(classlabels_orig[,2]),t(data_m))
      #myData<-cbind(as.factor(classgroup),as.factor(classlabels_orig[,2]),t(data_m))
      
      myData<-as.data.frame(myData)
      
      #save(myData,classgroup,classlabels_orig,data_m,file="myData.Rda")
      
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
   
        
        
        mzname<-rnames[get_label_ind]
        colnames(df.summary)<-c("Intensity","sd","number","se","Class")
        df.summary<-as.data.frame(df.summary)
        
        #df.summary$Class<-as.numeric(as.factor(df.summary$Class))
        
        
        ymax = df.summary$Intensity + 1.96*df.summary$se
        
        ymin = df.summary$Intensity - 1.96*df.summary$se
        
        df_write_temp<-cbind(mzname,df.summary[,c(5,3,1,2,4)],ymin,ymax)
        
        colnames(df_write_temp)<-c("Name","Class","Number of subjects","mean","Std.deviation","Std.error","lower.limit.95%CI","upper.limit.95%CI")
        return(df_write_temp)
      })
      
    }
    
  #  save(plot_res,file="plot_res.Rda")
    var_sum_mat<-{}
    #for(i in 1:length(plot_res))
  #  {
   #   var_sum_mat<-rbind(var_sum_mat,plot_res[[i]])
      
    #}
    var_sum_mat<-ldply(plot_res,rbind)
    #print(getwd())
    
    if(dir.exists("Tables")){
    write.table(var_sum_mat,file="Tables/data_summary.txt",sep="\t",row.names=FALSE)
    }else{
      
      write.table(var_sum_mat,file="data_summary.txt",sep="\t",row.names=FALSE)
    }
  }
  
}

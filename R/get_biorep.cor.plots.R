get_biorep.cor.plots <-
function(X=NA,Y=NA,feature_table_file=NA,parentoutput_dir=NA,
                               class_labels_file=NA,alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,
                               point.cex.val=4,legendlocation="topright",pca.ellipse=TRUE,
                               ellipse.conf.level=0.95,filename="QC_Correlation.Intragroup.Replicates",newdevice=TRUE,
                               lineplot.col.opt=c("grey57"),ylabel="Intensity",error.bar=TRUE,
                               cex.plots=0.8,
                               lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                               timeseries.lineplots=FALSE,name=NA,study.design="oneway",
                               alphabetical.order=TRUE,output.format="pdf",multiple.figures.perpanel=TRUE)
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
      
      classgroup<-paste(classlabels_orig[,2],":",classlabels_orig[,3],sep="") #classlabels_orig[,2]:classlabels_orig[,3]
      do_pca_anova=FALSE
    }else{
      classgroup<-classlabels_orig[,2]
      do_pca_anova=TRUE
    }
    
    col_class_levels<-levels(as.factor(classlabels_orig[,2]))
  }else{
    
    classgroup<-classlabels_orig[,2]
    col_class_levels<-levels(as.factor(classlabels_orig[,2]))
    
    do_pca_anova=TRUE
  }
  
  
  
  
  
  
  
  
  
  class_labels_levels<-levels(as.factor(classgroup))
  ordered_labels<-classgroup
  
  class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
  
  if(newdevice==TRUE){
    
    fname<-paste("",filename,".pdf",sep="")
    pdf(fname)
  }
  
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  
  
 # text(5,9,"Description:",font=2,col="blue")
  text(5, 8, "This PDF includes QC plots to evaluate intra-group heterogeneity\n based on the mean pairwise Pearson correlation between \nthe biological replicates in each group.\n\n\n Note: Wider error bars indicate higher intra-group heterogeneity.",cex=0.9,font=2)
 # text(5, 7, "The figures include: ")
  #text(5, 6, "a. pairwise PC score plots ")
  #text(5, 5, "b. scores for individual samples on each PC")
  #text(5, 4, "c. Lineplots using PC scores for data with repeated measurements")
  
  
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
                        
                    #  }
                      
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
      
      
      
      
    #  save(classlabels_orig,classgroup,data_m,file="myData.Rda")
      
      
      myData<-cbind(classgroup,as.data.frame(classlabels_orig[,2]),as.data.frame(classlabels_orig[,3]),t(data_m))
      
      myData<-as.data.frame(myData)
      
      #
      if(alphabetical.order==FALSE){
        myData[,2] <- factor(myData[,2], levels=unique(myData[,2]))
        myData[,3] <- factor(myData[,3], levels=unique(myData[,3]))
        
      }
      
      
      
      
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
      
      write.table(myData_sum_cor,file="myData_sum_cor.txt",sep="\t")
      myData_sum_cor$mean.Pearson.cor<-as.numeric(as.character(myData_sum_cor$mean.Pearson.cor))
      myData_sum_cor$lower.conf.limit<-as.numeric(as.character(myData_sum_cor$lower.conf.limit))
      myData_sum_cor$upper.conf.limit<-as.numeric(as.character(myData_sum_cor$upper.conf.limit))
      
     # save(myData_sum_cor,sizeval,cex.val,class_col_vec,shape_vec1,file="mysamp.Rda")
      group_num<-gsub(myData_sum_cor$group,pattern=":[a-z|A-Z|0-9|_|-]*",replacement="")
      myData_sum_cor<-cbind(myData_sum_cor,group_num)
      Factor2<-gsub(myData_sum_cor$group,pattern="[a-z|A-Z|0-9|-|_]*:",replacement="")
      myData_sum_cor<-cbind(myData_sum_cor,Factor2)
      
      plot_sample_cor<-ggplot(myData_sum_cor,aes(y=mean.Pearson.cor,x=Factor2,color=Factor2))+
        geom_point()+facet_wrap(~group_num, scale="free_x")+
        geom_errorbar(aes(ymin=lower.conf.limit,ymax=upper.conf.limit))+ 
        theme_bw() + theme(panel.border = element_blank(), 
                                                                                                                                       panel.grid.major = element_blank(),
                                                                                                                                       panel.grid.minor = element_blank(),
                                                                                                                                       axis.line = element_line(colour = "black",size=1),
                                                                                                                                       axis.text= element_text(size=14*cex.val), 
                                                                                                                                       axis.title=element_text(size=16*cex.val,face="bold"),
                                                                                                                                       strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                       strip.text = element_text(face="bold"))+scale_y_continuous(name="mean Pearson correlation\nbetween biological replicates in each group",breaks = scales::pretty_breaks(n = 10))
      #plot_sample_cor<-plot_sample_cor+scale_color_manual(values=unique(class_col_vec))
      #plot_sample_cor<-plot_sample_cor+scale_x_discrete(name ="Class")
      plot_sample_cor<-plot_sample_cor+theme(axis.text.x = element_text(angle = 45, hjust = 1))
      #guides(fill = guide_legend(title = FALSE),colour=FALSE)+
      #  guides(colour = guide_legend(override.aes = list(shape = unique(shape_vec1))))+
      
      print(plot_sample_cor)
      if(newdevice==TRUE){
        try(dev.off())
        }
      #
    }
  }
}

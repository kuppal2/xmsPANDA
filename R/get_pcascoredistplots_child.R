get_pcascoredistplots_child <-
function(X,Y,feature_table_file,parentoutput_dir,class_labels_file,
                                      sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300,
                                      alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,pca.cex.val=4,legendlocation="topright",pca.ellipse=TRUE,ellipse.conf.level=0.95,filename="all",paireddesign=NA,
                                      error.bar=TRUE,lineplot.col.opt="black",lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),newdevice=FALSE,timeseries.lineplots=FALSE,alphabetical.order=FALSE,pcascale=TRUE,pcacenter=TRUE,study.design="oneway",lme.modeltype="RI",cex.plots=0.8,ypos.adj.factor=0.5,...)
{
  
  analysistype=study.design
  
  if(is.na(X)==TRUE){
    data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
    
  }else{
    
    data_matrix<-X
  }
  
  
  
  if(typeof(Y)=="logical"){
    classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
    
  }else{
    
    classlabels<-Y
  }
  
  dir.create(parentoutput_dir,showWarnings = FALSE)
  setwd(parentoutput_dir)
  
  if(newdevice==TRUE){
    
    fname<-paste("pcaplots",filename,".pdf",sep="")
    pdf(fname)
    
    
  }
  
  
  
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
        
        rownames(data_matrix)<-Name
        
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
  
  data_m<-data_matrix[,-c(1:2)]
  
  data_m<-as.matrix(data_m)
  
  check_rnames<-rownames(data_matrix)
  
  
  mzvec<-data_matrix[,1]
  rtvec<-data_matrix[,2]
  
  if(length(check_rnames)<1){
    if(length(check_names)<1){
    rnames<-paste(mzvec,rtvec,sep="_")
    
    rownames(data_m)<-as.character(rnames)
    }
    
  }else{
  #if(length(check_names)<1){ 
   rownames(data_m)<-rownames(data_matrix)
    
  }
  #rnames<-paste(mzvec,rtvec,sep="_")
  
  #rownames(data_m)<-as.character(rnames)
  
  
  classlabels_orig<-classlabels
  
  #classlabelsorig<-classlabelsorig[match(rownames(X),classlabelsorig[,1]),]
  
 #save(data_m,classlabels,file="pcaclasslabels.Rda")
  
  if(analysistype=="twowayrepeat" | analysistype=="2wayrepeat" | analysistype=="onewayrepeat" | analysistype=="1wayrepeat"){
    
    pairedanalysis=TRUE
  }
  
  if(dim(classlabels)[2]>2){
    
    if(pairedanalysis==TRUE){
      
      if(is.na(paireddesign)==TRUE){
        paireddesign=classlabels_orig[,2]
      }
      classlabels_orig<-classlabels_orig[,-c(2)]
    }
    
   # print("here0")
    if(analysistype=="twoway" | analysistype=="2way" | analysistype=="twowayrepeat" | analysistype=="2wayrepeat"){
      if(dim(classlabels_orig)[2]>2){
        classgroup<-paste(classlabels_orig[,2],":",classlabels_orig[,3],sep="") #classlabels_orig[,2]:classlabels_orig[,3]
      }
      
     # print("here1")
     # print(head(classlabels_orig))
     # print(head(classgroup))
      
    }else{
      
      classgroup<-classlabels_orig[,2]
    }
    
    
    
    do_pca_anova=FALSE
  }else{
    
    if(analysistype=="regression"){
      
      classgroup<-rep("A",nrow(classlabels_orig))
      do_pca_anova=FALSE
      
    }else{
      classgroup<-classlabels_orig[,2]
      
      do_pca_anova=TRUE
    }
  }
  
  ##save(classgroup,file="classgroup.Rda")
  if(alphabetical.order==FALSE){
    classgroup <- factor(classgroup, levels=unique(classgroup))
  }
  
  # do_pca_anova=FALSE
  
  class_labels_levels<-levels(as.factor(classgroup))
  
  
  ordered_labels<-classgroup
  
  class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
  
  classlabels<-classlabels_orig
  
  if(is.na(col_vec)==TRUE)
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
        
        col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
      }else{
        if(sample.col.opt=="heat"){
          #col_vec<-heat.colors(256) #length(class_labels_levels))
          
          col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
        }else{
          if(sample.col.opt=="rainbow"){
            #col_vec<-heat.colors(256) #length(class_labels_levels))
            col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
            
            #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
          }else{
            
            if(sample.col.opt=="terrain"){
              #col_vec<-heat.colors(256) #length(class_labels_levels))
              #col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
              
              col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
            }else{
              
              if(sample.col.opt=="colorblind"){
                #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                
                if(length(class_labels_levels)<9){
                  
                  col_vec <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "#E64B35FF", "grey57")
                  
                }else{
                  
                  col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                             "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                  
                  #colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels)))
                }
                
                
              }else{
                
                check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                
                if(length(check_brewer)>0){
                  
                  sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                  
                  col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                  
                }else{
                  
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
                 #   if(length(sample.col.opt)==1){
                  #    col_vec <-rep(sample.col.opt,length(class_labels_levels))
                   # }else{
                      
                    #  colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                      
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
  }else{
    
    colfunc <-colorRampPalette(col_vec);col_vec<-colfunc(length(class_labels_levels))
  }
  
  class_col_vec=rep("red",nrow(classlabels))
  #    print(class_labels_levels)
  
  ordered_labels={}
  num_samps_group<-new("list")
  num_samps_group[[1]]<-0
  groupwiseindex<-new("list")
  groupwiseindex[[1]]<-0
  
  S<-new("list")
  
# save(class_labels_levels,col_vec,sample.col.opt,lineplot.col.opt,classlabels,file="pcaclass_labels_levels.Rda")
  for(c in 1:length(class_labels_levels))
  {
    
    classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
    #ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
    num_samps_group[[c]]<-length(classlabels_index)
    groupwiseindex[[c]]<-classlabels_index
    
    class_col_vec[which(classlabels[,2]==class_labels_levels[c])]<-col_vec[c]
    # S[[c]]<-cov(t(data_m[,c(classlabels_index)]))
    
  }
  
  
  
  sampleclass<-{}
  patientcolors<-{}
  
  classlabels<-as.data.frame(classlabels)
  
  
  for(c in 1:length(class_labels_levels)){
    
    num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
    
    
    sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
    
    patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group[[c]]))
  }
  
  if(length(mzvec)>4){
    max_per_row<-3
    
    
    par_rows<-ceiling(9/max_per_row)
    
  }else{
    max_per_row<-length(mzvec)
    par_rows<-1
  }
  
  file_ind<-0
  boxplots_fname<-paste("xyplots.pdf",sep="")
  
  suppressWarnings(dir.create("Tables",showWarnings = FALSE))
  t1<-table(classgroup)
  
  l1<-levels(as.factor(classgroup))
  
  
  patientcolors <- rep(col_vec[1:length(t1)], t1)
  
 # save(data_m,classgroup,legendlocation,filename,pcacenter,paireddesign,pcascale,col_vec,
  #     sample.col.opt,pca.cex.val,pca.ellipse,do_pca_anova,paireddesign,classlabels_orig,
   #   alphabetical.order,analysistype,file="pcad1A.Rda")
  
  #save(data_m,file="data_m.Rda")
  
  if(dim(data_m)[1]>2){
    
    
    res<-get_pca(X=data_m,samplelabels=classgroup,legendlocation=legendlocation,filename=filename,ncomp=10,pcacenter=pcacenter,pcascale=pcascale,legendcex=0.5,outloc=getwd(),
                 col_vec=col_vec,sample.col.opt=sample.col.opt,alphacol=0.3,class_levels=NA,pca.cex.val=pca.cex.val,pca.ellipse=pca.ellipse,do_pca_anova=do_pca_anova,
                 paireddesign=paireddesign,classlabelsorig=classlabels_orig,alphabetical.order=alphabetical.order,analysistype=analysistype,pairedanalysis=pairedanalysis,lme.modeltype=lme.modeltype)
    numcomp_plot<-res$numcomp_plot
    pcnum_limit=numcomp_plot
   
    pc_pval_vec<-res$pca_pval_vec
    res<-res$result
   
   
    
    fname<-paste("Tables/PCAloadings_",filename,"features.txt",sep="")
    
    loadings_res<-res$rotation
    scores_res<-res$x
    
    loadings_res<-round(loadings_res,2)
    scores_res<-round(scores_res,2)
    
    if(dim(loadings_res)[2]>10){
      loadings_res<-loadings_res[,c(1:10)]
      scores_res<-scores_res[,c(1:10)]
    }
    
    write.table(loadings_res,file=fname,sep="\t")
    
    fname<-paste("Tables/PCAscores_",filename,"features.txt",sep="")
    
    #save(classlabels_orig,scores_res,file="c2.Rda")
    
    classlabels_orig<-classlabels_orig[match(rownames(scores_res),classlabels_orig[,1]),]
    scores_res1<-cbind(classlabels_orig,scores_res)
    
    
    write.table(scores_res1,file=fname,sep="\t",row.names=FALSE)
    
    #pcnum_limit<-res$numcomp_plot #min(2,dim(scores_res)[2])
    
    get_pooled_sp<-function(n1,n2,S1,S2){a<-(n1-1)*S1;b<-(n2-1)*S2;c<-(n1+n2-2);return((a+b)/c)}
    
    
    if(analysistype=="regression"){
      
      
      df.sub1<-cbind(classlabels_orig[,2],scores_res[,1])
      
      
      colnames(df.sub1)<-c("Outcome","PCscore")
      df.sub1<-as.data.frame(df.sub1)
      
      s1<-summary(df.sub1$PCscore)
      sadj=(s1[5]-s1[3])*0.5
      
      plot_res<-suppressMessages(ggscatter(df.sub1,x="Outcome",y="PCscore",xlab="Outcome",ylab="PC1score",
                          title="Outcome vs PC1 scores scatter plot",col="darkblue",
                          palette="jco", shape = 20, size = 3, # Points color, shape and size
                          add = "reg.line",  # Add regressin line
                          add.params = list(color = "#0072B2", fill = "lightgray"), # Customize reg. line
                          conf.int = TRUE))+theme(plot.title = element_text(hjust = 0.5,size=10))+stat_cor(method = "spearman",
                                                                                                          label.y=max(df.sub1$PCscore+sadj))
      
      suppressMessages(print(plot_res))
      if(pcnum_limit>1){
      df.sub1<-cbind(classlabels_orig[,2],scores_res[,2])
      df.sub1<-as.data.frame(df.sub1)
      colnames(df.sub1)<-c("Outcome","PCscore")
      
      s1<-summary(df.sub1$PCscore)
      sadj=(s1[5]-s1[3])*0.5
      
      plot_res<-suppressMessages(ggscatter(df.sub1,x="Outcome",y="PCscore",xlab="Outcome",ylab="PC2score",
                          title="Outcome vs PC2 scores scatter plot",col="darkblue",
                          palette="jco", shape = 20, size = 3, # Points color, shape and size
                          add = "reg.line",  # Add regressin line
                          add.params = list(color = "#0072B2", fill = "lightgray"), # Customize reg. line
                          conf.int = TRUE))+theme(plot.title = element_text(hjust = 0.5,size=10))+stat_cor(method = "spearman",
                                                                                                          label.y=max(df.sub1$PCscore+sadj))
      
      
      suppressMessages(print(plot_res))
      }
      
      if(pcnum_limit>2){
      df.sub1<-cbind(classlabels_orig[,2],scores_res[,3])
      
      df.sub1<-as.data.frame(df.sub1)
      colnames(df.sub1)<-c("Outcome","PCscore")
      
      s1<-summary(df.sub1$PCscore)
      sadj=(s1[5]-s1[3])*0.5
      
      plot_res<-suppressMessages(ggscatter(df.sub1,x="Outcome",y="PCscore",xlab="Outcome",ylab="PC3score",
                          title="Outcome vs PC3 scores scatter plot",col="darkblue",
                          palette="jco", shape = 20, size = 3, # Points color, shape and size
                          add = "reg.line",  # Add regressin line
                          add.params = list(color = "#0072B2", fill = "lightgray"), # Customize reg. line
                          conf.int = TRUE))+theme(plot.title = element_text(hjust = 0.5,size=10))+stat_cor(method = "spearman",
                                                                                                          label.y=max(df.sub1$PCscore+sadj))
      
      suppressMessages(print(plot_res))
      }
      
      
    }else{
      
      class_levels<-levels(as.factor(classlabels_orig[,2]))
      
      
      if(alphabetical.order==FALSE){
        classlabels_orig[,2]<- factor(classlabels_orig[,2], levels=unique(classlabels_orig[,2]))
        
        if(dim(classlabels_orig)[2]>2){
          
          classlabels_orig[,3]<- factor(classlabels_orig[,3], levels=unique(classlabels_orig[,3]))
        }
      }
      
      
      
      #sample.col.opt=lineplot.col.opt
      
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
            
            col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
          }else{
            if(sample.col.opt=="heat"){
              #col_vec<-heat.colors(256) #length(class_labels_levels))
              
              col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
            }else{
              if(sample.col.opt=="rainbow"){
                #col_vec<-heat.colors(256) #length(class_labels_levels))
                col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
                
                #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
              }else{
                
                if(sample.col.opt=="terrain"){
                  #col_vec<-heat.colors(256) #length(class_labels_levels))
                  #col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
                  
                  col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
                }else{
                  
                  if(sample.col.opt=="colorblind"){
                    
                    if(length(class_labels_levels)<9){
                      
                      col_vec <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "#E64B35FF", "grey57")
                      
                    }else{
                      
                      col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                 "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                      
                      
                    }
                    
                    
                  }else{
                    
                    check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                    
                    if(length(check_brewer)>0){
                      
                      sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                      
                      col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                      
                    }else{
                      
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
                          
                          ###save(sample.col.opt,file="sample.col.opt.Rda")
                          ###save(class_labels_levels,file="class_labels_levels.Rda")
                       #   col_vec <-rep(sample.col.opt,length(class_labels_levels))
                      #  }else{
                          #fixthis2
                          #colfunc <-colorRampPalette(sample.col.opt)
                          
                       #   col_vec<-sample.col.opt #colfunc(length(class_labels_levels))
                          
                        #}
                        
                        if(length(sample.col.opt)==1){
                          col_vec <-rep(sample.col.opt,length(class_labels_levels))
                        }else{
                          
                          #if(length(sample.col.opt)>=length(class_labels_levels)){
                            
                            col_vec <-sample.col.opt
                            #col_vec <- rep(col_vec,length(class_labels_levels))
                            
                            
                          #}else{
                          #  colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
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
      df_matrix<-{}
      
      
      #lineplot.lty.option={}
      
      
      
      class_levels<-levels(as.factor(classlabels_orig[,2]))
      
      
      
      if(length(class_levels)>6){
        
        lineplot.lty.option=c(lineplot.lty.option,rep("solid",(length(class_levels)-6)))
      }
      #if(filename=="all")
      {
        
        
        fname<-paste("PCAloadings_",filename,"features.txt",sep="")
        
        loadings_res<-res$rotation
        scores_res<-res$x
        
        loadings_res<-round(loadings_res,2)
        scores_res<-round(scores_res,2)
        
        if(dim(loadings_res)[2]>10){
          loadings_res<-loadings_res[,c(1:10)]
          scores_res<-scores_res[,c(1:10)]
        }
        #   ##save(res,file="res.Rda")
        # write.table(loadings_res,file=fname,sep="\t")
        
        fname<-paste("PCAscores_",filename,"features.txt",sep="")
        
        scores_res1<-cbind(classlabels_orig,scores_res)
        
        
        # write.table(scores_res1,file=fname,sep="\t",row.names=FALSE)
        
        
      }
      
      #pcnum_limit<-res$numcomp_plot #min(5,dim(scores_res)[2])
      
      get_pooled_sp<-function(n1,n2,S1,S2){a<-(n1-1)*S1;b<-(n2-1)*S2;c<-(n1+n2-2);return((a+b)/c)}
      #S2<-cov(scores_res[16:46,1:2])
      #S1<-cov(scores_res[1:15,1:2])
      
      if(alphabetical.order==FALSE){
        
        
      }
      class_levels<-levels(as.factor(classlabels_orig[,2]))
      
      #if(length(class_labels_levels)==2)
      if(do_pca_anova==TRUE)
      {
        if(FALSE){
          h1<-hotelling.test(scores_res[1:num_samps_group[[1]],c(1)],scores_res[(1+num_samps_group[[1]]):(num_samps_group[[1]]+num_samps_group[[2]]),c(1)],shrinkage=FALSE)
          print("Hotelling test using PC1: ")
          print(h1)
          
          h2<-hotelling.test(scores_res[1:num_samps_group[[1]],c(2)],scores_res[(1+num_samps_group[[1]]):(num_samps_group[[1]]+num_samps_group[[2]]),c(2)],shrinkage=FALSE)
          print("Hotelling test using PC2: ")
          print(h2)
          
          h3<-hotelling.test(scores_res[1:num_samps_group[[1]],c(3)],scores_res[(1+num_samps_group[[1]]):(num_samps_group[[1]]+num_samps_group[[2]]),c(3)],shrinkage=FALSE)
          print("Hotelling test using PC3: ")
          print(h3)
        }
        
        pc_pval_vec<-{}
        
        
        
        for(pcnum in 1:pcnum_limit){
          
          pc1_pval<-anova(lm(cbind(scores_res[,pcnum])~classlabels_orig[,2]))
          pc1_pval<-round(pc1_pval[[5]][1],3)
          pc_pval_vec<-c(pc_pval_vec,pc1_pval)
        }
        
        
        
        
      }
      # #savenum_samps_group,scores_res,file="HotellingTestInput.Rda")
      #h1<-hotelling.test(scores_res[1:15,c(4:5)],scores_res[16:46,c(4:5)],shrinkage=FALSE)
      #h1<-hotelling.test(scores_res[1:15,c(2:3)],scores_res[16:46,c(2:3)],shrinkage=FALSE)
      
      df_matrix<-{}
      
      #keeps the class order same as in the input file; avoids arrangement by alphabetical order
      if(alphabetical.order==FALSE){
        classlabels_orig[,2] <- factor(classlabels_orig[,2], levels=unique(classlabels_orig[,2]))
      }
      
      class_col_vec=col_vec
      
     # save(res,classlabels_orig,file="df1.Rda")
      #save(pc_pval_vec,file="pc_pval_vec.Rda")
      
      #if(pairedanalysis==TRUE)
      {
        if(dim(classlabels_orig)[2]>2){
          class_levels<-levels(as.factor(classlabels_orig[,2]))
          t1<-table(classlabels_orig[,3])
          
          for(c in 1:length(class_levels))
          {
            
            
            #  class_col_vec[which(classlabels[,2]==class_levels[c])]<-col_vec[c]
            
            
          }
          #print("Doing two factors")
          
          # for(pc in 1:pcnum_limit)
          lapply(1:pcnum_limit,function(pc)
          {
            
            ylabel_text=paste("PC",pc,"score",sep="")
            
            classlabels_orig<-as.data.frame(classlabels_orig)
            classgroup<-as.data.frame(classgroup)
            res$x<-as.data.frame(res$x)
            if(pairedanalysis==FALSE){
              df<-cbind(res$x[,pc],classgroup,classlabels_orig[,2],classlabels_orig[,3])
              
              mzname<-paste("PC",pc," scores distribution with 95% confidence interval \n in each group using ",filename," feats vs factors",sep="")
              
            }else{
              
              
              df<-cbind(res$x[,pc],classgroup,classlabels_orig[,2],classlabels_orig[,3])
              
              mzname<-paste("PC",pc," scores distribution with 95% confidence interval \n in each group using ",filename," feats vs time",sep="")
              
            }
            
            fname<-paste("pc_",pc,"_scoreplot",".tiff",sep="")
            par_rows=3
            
            colnames(df)<-c("y","x","Class","time")
            df=as.data.frame(df)
            df$y<-as.numeric(as.character(df$y))
            
            
            
           #save(df,file="df.Rda")
            
            df_fname<-paste("dftable_pc",pc,".txt",sep="")
            
            if(pc>1){
              if(pc==2){
                PC2Score=df[,1]
                df_matrix<-cbind(df_matrix,PC2Score)
              }else{
                
                
                PC3Score=df[,1]
                df_matrix<-cbind(df_matrix,PC3Score)
                
                
              }
            }else{
              PC1Score<-df[,1]
              df_matrix<-cbind(df[,-c(1)],PC1Score)
              
            }
            
            if(alphabetical.order==FALSE){
              df$Class <- factor(df$Class, levels=unique(df$Class))
              df$time <- factor(df$time, levels=unique(df$time))
              
            }
            
            
            df.summary <- aggregate(df$y,
                                    by = list(df$Class,df$time),
                                    FUN = function(y) c(ysd = sd(y),
                                                        yse =sd(y)/sqrt(length(y)),
                                                        Score = mean(y),number=length(y)))
           # save(df.summary,file="df.summary.Rda")
            df.summary<-do.call(data.frame, df.summary)
            
            df.summary<-cbind(df.summary[,c(1,3,4,5,6,1,2)])
            
            colnames(df.summary)<-c("x","sd","se","PCscore","number","Class","time")
            df.summary<-as.data.frame(df.summary)
            
            ymax = df.summary$PCscore + 1.96* df.summary$se
            
            ymin =  df.summary$PCscore - 1.96* df.summary$se
            
            max_yval<-ceiling(max((df.summary$PCscore + (6* df.summary$se)),na.rm=TRUE)) #round(max( df.summary$Intensity+(4* df.summary$se),na.rm=TRUE))
            
            
            min_yval<-floor(min((df.summary$PCscore - (6*df.summary$se)),na.rm=TRUE))
            
            class_levels_time<-levels(as.factor(classlabels_orig[,3]))
            
            df.summary$Class<-as.numeric(as.factor(df.summary$Class))
            
            Class<-{}
            for(cnum in 1:length(class_levels)){
              
              Class<-c(Class,rep(class_levels[cnum],length(t1)))
              
              df.summary$Class[which(df.summary$Class==cnum)]<-class_levels[cnum]
              
            }
            
            df_time<-df.summary$time
            df.summary$time<-as.numeric(as.factor(df.summary$time))
            
            
            Class<-unique(Class)
            if(alphabetical.order==FALSE){
              df.summary$Class <- factor(df.summary$Class, levels=unique(df.summary$Class))
            }
            ###save(df.summary,file="df.summary.Rda")
            df.summary$x<- df.summary$time
            
            #savelist=ls(),file="pcadebugfactor.Rda")
            
            shape_vec=c(0:2,5:25)
            t1shape=levels(factor(df.summary$Class))
            
            shape_vec1<-df.summary$Class
            
            for(i in 1:length(t1shape)){
              shape_vec1[which(shape_vec1==t1shape[i])]<-shape_vec[i]
            }
            
            
            pca_plot_fname<-paste("pca_plots_",pc,".Rda",sep="")
            
            unique_class_col_vec=unique(class_col_vec)
            
            #save(df.summary,shape_vec1,df_time,col_vec,classlabels,
            #  classgroup,classlabels_orig,unique_class_col_vec,class_col_vec,ymin,ymax,pca.cex.val,df_time,file=pca_plot_fname)
            
            
            if(pairedanalysis==TRUE){
              
              
              sizeval=1.5
              #plot.margin=unit(c(10,8,8,8),"mm"),
              
              #+ labs(color = "Class")
              shape_vec1<-as.numeric(as.character(shape_vec1))
              
              
              
              
              plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(x), y = PCscore,
                                                                color = as.factor(Class),linetype=Class)) + 
                                           labs(color = "Class",linetype="Class") + 
                                           geom_point(size = pca.cex.val,shape=shape_vec1) + 
                                           geom_line(aes(group =as.factor(Class)),size=sizeval) + 
                                           geom_errorbar(aes(ymin = ymin, ymax = ymax),size=1,width=0.1) + 
                                           xlab("TimePoint") + theme_bw() + 
                                           theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(),
                                                 axis.line = element_line(colour = "black",size=sizeval),
                                                 axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=16*cex.plots,face="bold"),
                                                 
                                                 strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                 strip.text = element_text(face="bold")) + 
                                           scale_y_continuous(breaks = scales::pretty_breaks(n = 10)))
              
              
              
              if(length(unique_class_col_vec)==1){
                plot_res<-plot_res+scale_color_manual(values=rep(unique(class_col_vec),length(class_labels_levels)))
                
              }else{
                plot_res<-plot_res+scale_color_manual(values=unique(class_col_vec))
              }
              
              plot_res<-plot_res+scale_x_discrete(name ="Time",labels=as.factor(unique(df_time)))
              
              plot_res<-plot_res+theme(legend.text = element_text(size = 11*cex.plots))+guides(fill = guide_legend(title = NULL))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
              
              #if(FALSE)
              {
                plot_res<-plot_res+guides(fill = guide_legend(title = FALSE),colour=FALSE)+guides(colour = guide_legend(override.aes = list(shape = unique(shape_vec1))))
              }
              
              #  plot_res<-plot_res+guides(fill = guide_legend(title = FALSE),colour=FALSE)+guides(colour = guide_legend(override.aes = list(shape = unique(shape_vec1))))
              
              # plot_res+guides(color=FALSE,linetype=FALSE)
              
              
            }else{
              
              
              
              if(timeseries.lineplots==TRUE){
                
                
                shape_vec1<-as.numeric(as.character(shape_vec1))
                sizeval=1.5
                
                
                if(is.na(lineplot.lty.option)==TRUE){
                  plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(x), y = PCscore,colour = as.factor(Class),linetype=Class)) + geom_point(size = pca.cex.val,shape=shape_vec1) + geom_line(aes(group =Class),size=sizeval)  + labs(colour = "Class") + geom_errorbar(aes(ymin = ymin, ymax = ymax),size=1,width=0.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                                                                                                                                                                                                                                                                                                                           axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=16*cex.plots,face="bold"),
                                                                                                                                                                                                                                                                                                                                                           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                                                                                                                                                                                                                                           strip.text = element_text(face="bold")) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)))
                }else{
                  plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(x), y = PCscore,colour = as.factor(Class))) + geom_point(size = pca.cex.val,shape=shape_vec1) + geom_line(aes(group =Class),size=sizeval)  + labs(colour = "Class") + geom_errorbar(aes(ymin = ymin, ymax = ymax),size=1,width=0.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                                                                                                                                                                                                                                                                                                            axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=16*cex.plots,face="bold"),
                                                                                                                                                                                                                                                                                                                                            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                                                                                                                                                                                                                            strip.text = element_text(face="bold")) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)))
                }
                
                
                if(length(unique_class_col_vec)==1){
                  plot_res<-plot_res+scale_color_manual(values=rep(unique(class_col_vec),length(class_labels_levels)))
                  
                }else{
                  plot_res<-plot_res+scale_color_manual(values=unique(class_col_vec))
                }
                plot_res<-plot_res+scale_x_discrete(name ="Time",labels=unique(df_time))
                plot_res<-plot_res+theme(legend.text = element_text(size = 11*cex.plots))+guides(fill = guide_legend(title = NULL))
                
                #legend.labs=levels(factor(df.summary$Class)), color = "strata")
                #plot.margin = unit(c(25,25,5.5,28), "pt"),
                
                uniq_class<-unique(df.summary$Class)
                plot_res<-plot_res+guides(fill = guide_legend(title = FALSE),colour=FALSE)+guides(colour = guide_legend(override.aes = list(shape = unique(shape_vec1))))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                
                
                #    plot_res<-plot_res+guides(linetype=FALSE)
                
              }else{
                sizeval=1.5
                
                
                
                plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(x), y = PCscore,color = as.factor(Class))) + geom_point(size = pca.cex.val,shape=shape_vec1) + labs(color = "Class") +  geom_errorbar(aes(ymin = ymin, ymax = ymax),size=1,width=0.1) + xlab("") + scale_color_manual(values=unique(class_col_vec)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                                                                                                                                                                                                                                                                                                                         axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=16*cex.plots,face="bold"),
                                                                                                                                                                                                                                                                                                                                                         strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                                                                                                                                                                                                                                         strip.text = element_text(face="bold")) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_fill_manual(values=unique(class_col_vec)))
                
                
                plot_res<-plot_res+theme(legend.text = element_text(size = 11*cex.plots))+guides(fill = guide_legend(title = NULL))
                
                
                plot_res<-plot_res+scale_x_discrete(name ="Factor2",labels=as.factor(unique(df_time)))
                
                plot_res<-plot_res+guides(fill = guide_legend(title = FALSE),colour=FALSE)+guides(colour = guide_legend(override.aes = list(shape = unique(shape_vec1))))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                
              }
              
              
              
            }
            
            plot_res<-plot_res + ggtitle(mzname)
            
            plot_res<-plot_res + theme(plot.title = element_text(size=8))
            # x axis title
            plot_res<-plot_res + theme(axis.title.x = element_text(size=14*cex.plots))
            # y axis title
            plot_res<-plot_res + theme(axis.title.y = element_text(size=14*cex.plots))+theme(axis.text.x = element_text(angle = 45, hjust = 1,size=11**cex.plots))
            
            # plot_res<-plot_res + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
            
            
            #plot_res
            print(plot_res)
            
          })
          
          ##############
        }else{
          
          if(alphabetical.order==FALSE){
            classgroup <- factor(classgroup, levels=unique(classgroup))
          }
          class_levels<-levels(as.factor(classgroup))
          t1<-table(classgroup)
          
          
          for(c in 1:length(class_levels))
          {
            
            
            class_col_vec[which(as.factor(classgroup)==class_levels[c])]<-col_vec[c]
            
            
          }
          
          # for(pc in 1:pcnum_limit)
          
          
          #save(classgroup,class_levels,pcnum_limit,classlabels_orig,res,pairedanalysis,col_vec,pc_pval_vec,file="debuglpA.Rda")
          lapply(1:pcnum_limit,function(pc)
          {
            ylabel_text<-paste("PC",pc,"score",sep="")
            df<-cbind(res$x[,pc],as.character(classlabels_orig[,2]),as.character(classlabels_orig[,2]))
            
            df=as.data.frame(df)
            
            if(pairedanalysis==FALSE){
              #  mzname<-paste("PC",pc," scores distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs factors (p=",pc_pval_vec[pc],")",sep="")
              
              mzname<-paste("PC",pc," scores distribution with 95% confidence interval in each group using ",filename," feats vs factors (p=",pc_pval_vec[pc],")",sep="")
            }else{
              
              #mzname<-paste("PC",pc," scores distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs time",sep="")
              mzname<-paste("PC",pc," scores distribution with 95% confidence interval in each group using ",filename," feats vs time (p=",pc_pval_vec[pc],")",sep="")
              
            }
            
            fname<-paste("pc_",pc,"_scoreplot",".tiff",sep="")
            par_rows=3
            
            colnames(df)<-c("y","x","Class")
            df=as.data.frame(df)
            df$y<-as.numeric(as.character(df$y))
            
            
            
            #savedf,file="pcadf.Rda")
            
            df_fname<-paste("dftable_pc",pc,".txt",sep="")
            # write.table(df,file=df_fname,sep="\t",row.names=FALSE)
            
            if(pc>1){
              if(pc==2){
                PC2Score=df[,1]
                df_matrix<-cbind(df_matrix,PC2Score)
              }else{
                
                
                PC3Score=df[,1]
                df_matrix<-cbind(df_matrix,PC3Score)
                
                
              }
            }else{
              PC1Score<-df[,1]
              df_matrix<-cbind(df[,-c(1)],PC1Score)
              
            }
            
            
            if(FALSE){
              df.summary <- aggregate(df$y,
                                      by = list(df$Class),
                                      FUN = function(y) c(ymin = quantile(y,0.25),
                                                          ymax = quantile(y,0.75),
                                                          Score = median(y),number=length(y)))
            }
            
            if(alphabetical.order==FALSE){
              df$Class <- factor(df$Class, levels=unique(df$Class))
              
              
            }
            df.summary <- aggregate(df$y,
                                    by = list(df$Class),
                                    FUN = function(y) c(ysd = sd(y),
                                                        yse =sd(y)/sqrt(length(y)),
                                                        Score = mean(y),number=length(y)))
            
            df.summary<-do.call(data.frame, df.summary)
            
            df.summary<-cbind(df.summary[,c(1,2:5,1)])
            
            colnames(df.summary)<-c("x","sd","se","PCscore","number","Class")
            df.summary<-as.data.frame(df.summary)
            ymax = df.summary$PCscore + 1.96* df.summary$se
            
            df.summary<-as.data.frame(df.summary)
            
            ymin =  df.summary$PCscore - 1.96* df.summary$se
            
            max_yval<-ceiling(max((df.summary$PCscore + (6* df.summary$se)),na.rm=TRUE)) #round(max( df.summary$Intensity+(4* df.summary$se),na.rm=TRUE))
            
            
            min_yval<-floor(min((df.summary$PCscore - (6*df.summary$se)),na.rm=TRUE))
            
            # ##save(df.summary,ymin,ymax,classlabels_orig,class_levels,file="df.summary.Rda")
            
            time.hour<- c(unique(as.character(classlabels_orig[,2])),
                          unique(as.character(classlabels_orig[,2])))
            
            
            ##save(df.summary,class_col_vec,pca.cex.val,file="pca_plots.Rda")
            
            
            sizeval=1.5
            
            #save(df.summary,df,pca.cex.val,class_levels,class_labels_levels,pairedanalysis,timeseries.lineplots,ymin,ymax,
            #class_col_vec,sizeval,classlabels_orig,file="pcad1plotB.Rda")
            
            
            
            {
              if(pairedanalysis==FALSE){
                
                if(timeseries.lineplots==TRUE){
                  sizeval=1.5
                  plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(Class), y = PCscore,group=1)) 
                                             + geom_point(size = cex.plots) + geom_line(size=sizeval)
                                             + geom_errorbar(aes(ymin = ymin, ymax = ymax),size=1,width=0.1) + xlab("") 
                                             + scale_color_manual(values=unique(class_col_vec)) + scale_size_manual(values=c(sizeval)) 
                                             + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                                  axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=16*cex.plots,face="bold"),
                                                                  strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                  strip.text = element_text(face="bold")) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)))
                  
                  plot_res<-plot_res+scale_x_discrete(name ="Time",labels=unique(df.summary$Class))
                  plot_res<-plot_res+theme(legend.text = element_text(size = 12))+guides(fill = guide_legend(title = NULL))
                  
                  
                }else{
                  
                  sizeval=1.5
                  # #save(df.summary,pca.cex.val,ymin,ymax,sizeval,class_col_vec,cex.plots,file="debug1.Rda")
                  if(length(unique(class_col_vec))==1){  
                    if(length(unique(class_col_vec))==1 && analysistype=="regression"){
                      
                      df.sub1<-df.summary[,c("x","PCscore")]
                      df.sub1<-apply(df.sub1,2,as.numeric)
                      df.sub1<-as.data.frame(df.sub1)
                      
                      s1=summary(df.sub1$PCscore)
                      
                      sadj=(s1[5]-s1[3])*ypos.adj.factor
                      
                      plot_res<-suppressMessages(ggscatter(df.sub1,x="x",y="PCscore",xlab="Outcome",col=class_col_vec,
                                          palette="jco", shape = 20, size = 3, # Points color, shape and size
                                          add = "reg.line",  # Add regressin line
                                          add.params = list(color = "#0072B2", fill = "lightgray"), # Customize reg. line
                                          conf.int = TRUE))+stat_cor(method = "spearman",
                                 #aes(label = paste(..r.label.., ..p.label..),
                                 label.x = 3,label.y=max(df.sub1$PCscore+sadj)
                      )
                      
                      
                      
                      
                    }else{
                      plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(Class), y = PCscore)
                                                        + geom_point(size = pca.cex.val) +
                                                          geom_errorbar(aes(ymin = ymin, ymax = ymax),size=1,width=0.1) + xlab("") +
                                                          theme_bw() +
                                                          scale_size_manual(values=c(sizeval)) +
                                                          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                                axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=16*cex.plots,face="bold"),
                                                                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                strip.text = element_text(face="bold")) +
                                                          scale_y_continuous(breaks = scales::pretty_breaks(n = 10))))
                      
                    }
                    
                    #   if(length(unique_class_col_vec)==1){
                    plot_res<-plot_res+scale_color_manual(values=class_col_vec) #rep(unique(class_col_vec),length(class_labels_levels)))
                    
                    plot_res<-plot_res+theme(legend.title = element_blank())
                    
                    # plot_res<-plot_res+theme(legend=FALSE)
                  }else{
                    #   print("here1")
                    
                    plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(Class), y = PCscore,color = as.factor(Class))) 
                                               + geom_point(size = pca.cex.val) + labs(color = "Class") + 
                                                 geom_errorbar(aes(ymin = ymin, ymax = ymax),size=1,width=0.1) + xlab("") +
                                                 scale_color_manual(values=unique(class_col_vec)) + theme_bw() +
                                                 scale_size_manual(values=c(sizeval)) + 
                                                 theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                       axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=16*cex.plots,face="bold"),
                                                       strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                       strip.text = element_text(face="bold")) + 
                                                 scale_y_continuous(breaks = scales::pretty_breaks(n = 10)))
                    plot_res<-plot_res+scale_x_discrete(name ="",labels=as.factor(class_levels))
                    plot_res<-plot_res+theme(legend.text = element_text(size = 11*cex.plots))+guides(fill = guide_legend(title = NULL))
                    
                  }
                  
                  
                  #,legend.position="none",
                }
              }else{
                #print("here2")
                
                 #save(df.summary,pca.cex.val,ymin,ymax,class_col_vec,class_levels,classgroup,file="t1.Rda")
                sizeval=1.5
                plot_res<-suppressWarnings(ggplot(df.summary, aes(x = as.factor(x), y = PCscore,group=1)) +  geom_point(size = pca.cex.val) + geom_line(size=sizeval) + geom_errorbar(aes(ymin = ymin, ymax = ymax),size=1,width=0.1) + scale_color_manual(values=unique(class_col_vec)) + scale_size_manual(values=c(sizeval)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval),
                                                                                                                                                                                                                                                                                                                                                     axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=16*cex.plots,face="bold"),
                                                                                                                                                                                                                                                                                                                                                     strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                                                                                                                                                                                                                                     strip.text = element_text(face="bold")) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)))
                
                #plot_res<-plot_res+suppressWarnings(scale_x_discrete(name ="Time",labels=unique(Class)))
                plot_res<-plot_res+scale_x_discrete(name ="",labels=as.factor(class_levels))
                plot_res<-plot_res+theme(legend.text = element_text(size = 11*cex.plots))+guides(fill = guide_legend(title = NULL))
                
                
                
              }
            }
            
            
            plot_res<-plot_res + ggtitle(mzname)
            
            plot_res<-plot_res + theme(plot.title = element_text(size=8))
            # x axis title
            plot_res<-plot_res + theme(axis.title.x = element_text(size=14*cex.plots))
            # y axis title
            plot_res<-plot_res + theme(axis.title.y = element_text(size=14*cex.plots))+theme(axis.text.x = element_text(angle = 45, hjust = 1,size=11*cex.plots))
            
            #    plot_res<-plot_res + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
            
            suppressMessages(print(plot_res))
            
          })
        }
        
      }
    }
  }
  
  df_fname<-paste("PC_score_distribution_matrix_",filename,"features.txt",sep="")
  
  if(newdevice==TRUE){
    
    dev.off()
  }
  
  return(res)
}

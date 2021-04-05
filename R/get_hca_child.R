get_hca_child <-
function(feature_table_file,parentoutput_dir,class_labels_file,X=NA,Y=NA,heatmap.col.opt="RdBu",cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
                        sample.col.opt="rainbow",plots.width=8,plots.height=8,plots.res=600, plots.type="cairo", alphacol=0.3, 
                        hca_type,newdevice=FALSE,input.type="intensity",mainlab="",cexRow=0.5, cexCol=0.5,
                        plot.bycluster=FALSE,color.rows=FALSE,similarity.matrix="correlation",deepsplit=2,minclustsize=10,mergeCutHeight=0.1,
                        num_nodes=2,alphabetical.order=FALSE,pairedanalysis=FALSE,cutree.method="halfheight",study.design="multiclass",labRow.value = FALSE,
                        labCol.value = FALSE,power_val=6,row.col.opt="journal",show.silhouette=FALSE,cexLegend=0.7,ylab_text="",xlab_text="")
{
  
  
  suppressMessages(library(WGCNA))
  suppressMessages(library(flashClust))
  suppressMessages(library(gplots))
  suppressMessages(library(cluster))
  suppressMessages(library(mclust))
  
  suppressMessages(library(RColorBrewer))
  
  analysistype=study.design
  
  if(typeof(X)=="logical"){
    data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
 
    
   }else{
    data_matrix<-X
    rm(X)
    
  }
  
  
  mycl_metabs={}
  mycl_samples={}    
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
        
        rownames(data_matrix)<-as.character(Name)
      }else{
        
        if(length(check_names1)>0 & length(check_names2)>0){
          
          check_ind<-gregexpr(cnames,pattern="^name$")
          check_ind<-which(check_ind>0)
          Name<-as.character(data_matrix[,check_ind])
          data_matrix<-data_matrix[,-check_ind]
          names_with_mz_time=cbind(Name,data_matrix$mz,data_matrix$time)
          colnames(names_with_mz_time)<-c("Name","mz","time")
          names_with_mz_time<-as.data.frame(names_with_mz_time)
          data_matrix<-as.data.frame
          rownames(data_matrix)<-as.character(Name)
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
  col_metabs=color.rows
  
  data_m<-data_matrix[,-c(1:2)]
  
  data_m<-as.matrix(data_m)
  
  rownames(data_m)<-rownames(data_matrix)
  col_samples<-TRUE
  
  suppressWarnings(dir.create("Tables",showWarnings = FALSE))
  
  
  if(labRow.value==TRUE){
    
    labRow.value=rownames(data_m)
  }else{
    
    labRow.value=NA
    
  }
  
  if(labCol.value==TRUE){
    
    labCol.value=colnames(data_m)
  }else{
    
    labCol.value=NA
  }
  
  if(is.na(cexLegend)==TRUE){
    
    hca.show.legend=FALSE
  }else{
    hca.show.legend=TRUE
    
  }
  
  if(typeof(Y)=="logical"){
    if(is.na(class_labels_file)==FALSE){
      classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
    }else{
      
      classlabels<-rep("classA",dim(data_m)[2])
      classlabels<-cbind(classlabels,classlabels)
      
      col_samples<-FALSE
      
    }
  }else{
    classlabels<-Y
    
  }
  
  
  classlabels_orig<-classlabels
  
  if(analysistype=="twowayrepeat" | analysistype=="2wayrepeat" | analysistype=="onewayrepeat" | analysistype=="1wayrepeat"){
    
    pairedanalysis=TRUE
  }
  
  
  if(dim(classlabels)[2]>2){
    
    if(pairedanalysis==TRUE){
      
      paireddesign=classlabels_orig[,2]
      
      classlabels_orig<-classlabels_orig[,-c(2)]
      
      classlabels<-classlabels[,-c(2)]
    }
    
    if(dim(classlabels_orig)[2]>2){
      
      if(analysistype=="twowayrepeat" | analysistype=="2wayrepeat" | analysistype=="twoway" | analysistype=="2way" | analysistype=="twowayanova" | analysistype=="twowayanovarepeat"){
        
        classgroup<-paste(classlabels_orig[,2],":",classlabels_orig[,3],sep="")
        
      }else{
        
        classgroup<-classlabels_orig[,2]
      }
      
    }else{
      
      classgroup<-classlabels_orig[,2]
    }
    
  }else{
    
    classgroup<-classlabels_orig[,2]
    
    
  }
  classlabels[,2]<-classgroup

  
  names(classgroup)<-classlabels[,1]
  
  #patientcolors<-rep("green",dim(data_m)[2])
  
  #keeps the class order same as in the input file; avoids arrangement by alphabetical order
  if(alphabetical.order==FALSE){
    classlabels[,2] <- factor(classlabels[,2], levels=unique(classlabels[,2]))
  }
  
  
  class_labels_levels<-levels(as.factor(classlabels[,2]))
  ordered_labels<-classlabels[,2]
  
  #class_label_alphabets<-c("A","B","C","D","E","F","G","H","I","J","K","L","M")
  class_label_alphabets<-class_labels_levels #paste("C",1:length(class_labels_levels),sep="")
  
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
                
                # col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
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
                  
                 # if(length(sample.col.opt)==1){
                  #  col_vec <-rep(sample.col.opt,length(class_labels_levels))
                  #}else{
                    
                   # colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                    
                   # col_vec <-sample.col.opt
                    #col_vec <- rep(col_vec,length(class_labels_levels))
                    
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
  col_samples=FALSE
  
  
  if(analysismode=="classification")
  {
    #save(classlabels,class_labels_levels,class_label_alphabets,ordered_labels,file="hca2.Rda")
    
    sampleclass<-{}
    patientcolors<-rep("green",nrow(classlabels))
    #print(classlabels)
    classlabels<-as.data.frame(classlabels)
    f<-factor(classlabels[,1])
    
    col_samples=TRUE
    for(c in 1:length(class_labels_levels)){
      
      num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
      
      #classlabels<-c(classlabels,rep(paste("Class",class_label_alphabets,sep=""),num_samps_group_cur))
      #,rep("ClassB",num_samps_group[[2]]),rep("ClassC",num_samps_group[[3]]))
      sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
      
      #patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group_cur))
      patientcolors[which(ordered_labels==class_labels_levels[c])]<-col_vec[c]
      
    }
    
    
    
  }
  
  if(heatmap.col.opt=="RdBu"){
    
    heatmap.col.opt="redblue"
  }
  
  heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
  heatmap_cols<-rev(heatmap_cols)
  
  if(heatmap.col.opt=="topo"){
    heatmap_cols<-topo.colors(256)
    heatmap_cols<-rev(heatmap_cols)
  }else{
    if(heatmap.col.opt=="heat"){
      heatmap_cols<-heat.colors(256)
      heatmap_cols<-rev(heatmap_cols)
    }else{
      
      if(heatmap.col.opt=="yellowblue"){
        
        heatmap_cols<-colorRampPalette(c("yellow","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
        #heatmap_cols<-blue2yellow(256) #colorRampPalette(c("yellow","blue"))(256)
        heatmap_cols<-rev(heatmap_cols)
      }else{
        
        if(heatmap.col.opt=="redblue"){
          
          heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
          heatmap_cols<-rev(heatmap_cols)
        }else{
          
          #my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
          if(heatmap.col.opt=="redyellowgreen"){
            
            heatmap_cols <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
            heatmap_cols<-rev(heatmap_cols)
          }else{
            if(heatmap.col.opt=="yellowwhiteblue"){
              
              heatmap_cols<-colorRampPalette(c("yellow2","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
              heatmap_cols<-rev(heatmap_cols)
            }else{
              
              if(heatmap.col.opt=="redwhiteblue"){
                
                heatmap_cols<-colorRampPalette(c("red","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                heatmap_cols<-rev(heatmap_cols)
              }else{
                
                
                if(length(grep(heatmap.col.opt,pattern = "brewer."))>0){
                 
                  heatmap.col.opt<-gsub(heatmap.col.opt,pattern="brewer.",replacement="")
                
                heatmap_cols <- colorRampPalette(brewer.pal(10, heatmap.col.opt))(256)
                heatmap_cols<-rev(heatmap_cols)
                }else{
                  
                  if(heatmap.col.opt=="bluered"){
                    
                    heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                   # heatmap_cols<-rev(heatmap_cols)
                  }else{
                    
                    if(heatmap.col.opt=="blueorange"){
                      
                      heatmap_cols <- colorRampPalette(c("orange","blue"))(256)
                      # heatmap_cols<-rev(heatmap_cols)
                    }else{
                      
                      if(heatmap.col.opt=="orangeblue"){
                        
                        heatmap_cols <- colorRampPalette(c("orange","blue"))(256)
                        heatmap_cols<-rev(heatmap_cols)
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
  
  # print(classlabels)
  #print(patientcolors)
  col_vec2 <- topo.colors(255, alpha=alphacol)
  
  mergeCutheight=mergeCutHeight 
  
  
#  save(data_m,heatmap_cols,patientcolors,file="hcaDfalse.Rda")
  
  if(input.type=="intensity"){
    
    
    if(similarity.matrix=="TOM"){
      ##save(data_m,cor.method,num_nodes,file="debugcorselfeats.Rda")
      simmat<-WGCNA::cor(t(data_m),nThreads=num_nodes,method=cor.method,use = 'p')
      write.table(simmat,file="Tables/pairwisecorrelation_selectedfeatures.txt",sep="\t",row.names=TRUE)
      
      if(is.na(power_val)==TRUE){
              #data_m<-log2(data_m+1)
              powers = c(c(1:10), seq(from = 12, to=20, by=2))
              
              
              sft = try(pickSoftThreshold.fromSimilarity(similarity=simmat, powerVector = powers, verbose = 0),silent=TRUE)
              #power_val=sft$powerEstimate
              
              if(is(sft,"try-error")){
                power_val=6
              }else{
                power_val=sft$powerEstimate
              }
              if(is.na(power_val)==TRUE){
                power_val=6 
              }

      }else{
        if(is.na(power_val)==TRUE){
          power_val=6 
        }
        
      }
      
      
      ADJdataOne<-adjacency.fromSimilarity(similarity=simmat,power=power_val)
      
      TOM = TOMsimilarity(ADJdataOne, TOMType="signed") # specify network type
      dissTOMCormat = 1-TOM
      #dissTOMCormat=TOMdist(ADJdataOne) #(1-global_cor)
      hr = fastcluster::hclust(as.dist(dissTOMCormat),method="complete");
      
      
      
      pdf("metabplot.pdf")
      plot(hr,labels=F,main="Dendrogram")
      dev.off()
      
      mycl_metabs <-cutreeDynamic(hr,distM= dissTOMCormat,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = FALSE, pamStage=TRUE,verbose=0)
      
      m2=try(mergeCloseModules(t(data_m),colors=mycl_metabs,cutHeight=mergeCutheight),silent=TRUE)
      
     # save(mycl_metabs,m2,data_m,mergeCutheight,file="d1.Rda")
      
      if(is(m2,"try-error")){
        
        mod_list<-mycl_metabs
      }else{
        
        mod_list<-as.numeric(m2$colors)
      }
      
      mycl_metabs<-mod_list
      
      s2.metab<- silhouette(mycl_metabs,dmatrix=dissTOMCormat)
      
      
      ###samples
      simmat<-WGCNA::cor((data_m),nThreads=num_nodes,method=cor.method,use = 'p')
      
      write.table(simmat,file="Tables/pairwisecorrelation_samples.txt",sep="\t",row.names=TRUE)
      
      
      if(is.na(power_val)==TRUE){
        #data_m<-log2(data_m+1)
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        
        
        sft = try(pickSoftThreshold.fromSimilarity(similarity=simmat, powerVector = powers, verbose = 0),silent=TRUE)
        #power_val=sft$powerEstimate
        
        if(is(sft,"try-error")){
          power_val=6
        }else{
          power_val=sft$powerEstimate
        }
        if(is.na(power_val)==TRUE){
          power_val=6 
        }
        
      }else{
        if(is.na(power_val)==TRUE){
          power_val=6 
        }
        
      }
      
      
      ADJdataOne<-adjacency.fromSimilarity(similarity=simmat,power=power_val)
      
      TOM = TOMsimilarity(ADJdataOne, TOMType="signed") # specify network type
      dissTOMCormat = 1-TOM
     # dissTOMCormat=TOMdist(ADJdataOne) #(1-global_cor)
      hc = fastcluster::hclust(as.dist(dissTOMCormat),method="complete");
      
      pdf("sampleplot.pdf")
      plot(hc,labels=F,main="Dendrogram")
      dev.off()
      
      
      mycl_samples <-cutreeDynamic(hc,distM= dissTOMCormat,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = FALSE, pamStage=FALSE,verbose=0)
      
      m3=try(mergeCloseModules((data_m),colors=mycl_samples,cutHeight=mergeCutheight),silent=TRUE)
      
      if(is(m3,"try-error")){
        
        mod_list2<-mycl_samples
      }else{
        
        mod_list2<-as.numeric(m3$colors)
      }
      
      mycl_samples<-mod_list2
      
      s1.samp<- silhouette(mycl_samples,dmatrix=dissTOMCormat)
      
      
      
      
    }else{
      
      
      simmat=WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs")
      distc_m<-1-simmat
      
      write.table(simmat,file="Tables/pairwisecorrelation_samples.txt",sep="\t",row.names=TRUE)
      
      #save(data_m,file="data_m.Rda")
      simmat=WGCNA::cor(t(data_m),method=cor.method,use="pairwise.complete.obs")
      distr_m<-1-simmat
      write.table(simmat,file="Tables/pairwisecorrelation_selectedfeatures.txt",sep="\t",row.names=TRUE)
      
      
      
      distc<-as.dist(distc_m)
      distr<-as.dist(distr_m)
      
      
      
      if(cutree.method=="dynamic"){
        
       # save(distc,distc_m,deepsplit,minclustsize,mergeCutHeight,data_m,file="debugdynamic.Rda")
        hc = fastcluster::hclust(distc,method="complete")
        
        mycl_samples <-cutreeDynamic(hc,distM= distc_m,deepSplit=deepsplit, method="hybrid",minClusterSize=minclustsize, pamRespectsDendro = FALSE, pamStage=FALSE,verbose=0)
        
        m3=try(mergeCloseModules((data_m),colors=mycl_samples,cutHeight=mergeCutHeight,verbose = 0),silent=TRUE)
        
        if(is(m3,"try-error")){
          
          mod_list2<-mycl_samples
        }else{
          
          mod_list2<-as.numeric(m3$colors)
        }
        
        mycl_samples<-mod_list2
        names(mycl_samples)<-hc$labels
        hr = fastcluster::hclust(distr,method="complete")
        
        mycl_metabs <-cutreeDynamic(hr,distM=distr_m,deepSplit=deepsplit, method="hybrid",minClusterSize=minclustsize, pamRespectsDendro = FALSE, pamStage=FALSE,verbose=0)
        
        m3=try(mergeCloseModules(t(data_m),colors=mycl_metabs,cutHeight=mergeCutHeight,verbose = 0),silent=TRUE)
        
        if(is(m3,"try-error")){
          
          mod_list2<-mycl_metabs
        }else{
          
          mod_list2<-as.numeric(m3$colors)
        }
        
        mycl_metabs<-mod_list2
        
        names(mycl_metabs)<-hr$labels
      }else{
        
        hc <- try(hclust(distc),silent=TRUE) #samples
        hr <- try(hclust(distr),silent=TRUE) #metabolites
        
        mycl_samples <- cutree(hc, h=max(hc$height)/2)
        mycl_metabs <- cutree(hr, h=max(hr$height)/2)
        
      }
      
   # save(hr,hc,deepsplit,minclustsize,mergeCutheight,mycl_metabs,distr_m,data_m,distr,patientcolors,file="hr_hcC.Rda")
      
      s1.samp<- silhouette(mycl_samples,distc)
      s2.metab<- silhouette(mycl_metabs,distr)
      
    }
    
    
  }else{
    
    
  }
  

  if(FALSE)
  {
    ##save(hc,file="hc.Rda")
    ##save(hr,file="hr.Rda")
    ##save(distc,file="distc.Rda")
    ##save(distr,file="distr.Rda")
    ##save(s1.samp,file="s1samp.Rda")
    ##save(s2.metab,file="s2metab.Rda")
  }
  
  if(is.na(s1.samp)==FALSE){
    s1.samp=round(mean(s1.samp[,3],na.rm=TRUE),2)
  }
  if(is.na(s2.metab)==FALSE){
    s2.metab=try(round(mean(s2.metab[,3],na.rm=TRUE),2),silent=TRUE)
  }
  
  classgroup_levels<-levels(as.factor(classgroup))
  
  classgroup<-as.data.frame(classgroup)
  colnames(classgroup)<-"factor_levels"
  
  class_df <- data.frame(factor_levels = classgroup_levels,
                         Class = seq(1,length(classgroup_levels)),
                         stringsAsFactors = FALSE)
  
  
  dt3 <- merge(classgroup, class_df, by = "factor_levels", all.x = TRUE)
  
  #save(classgroup,mycl_samples,file="hcad.Rda")
  ari_val<-try(round(adjustedRandIndex(x=classgroup[,1], y=mycl_samples),2),silent=TRUE)
  
  
  if(show.silhouette==FALSE){
  
  mainlab1<-"" #paste("HCA using ",mainlab," selected features",sep="")
  }else{
  mainlab1<-"" #paste("Average  width samples: ",s1.samp,"\n ","Average Silhouette width features: ", s2.metab," \n ", "Adjusted Rand index (comparison with true class labels): ",ari_val,sep="")
  }
  
  #rownames(data_m)<-NULL
  #colnames(data_m)<-NULL
  
  
  sample.col.opt=row.col.opt
  
  if(is(hr,"try-error") || is(hc,"try-error")){
    
    print("Hierarchical clustering can not be performed. ")
  }else{
    heatmap_file<-paste("heatmap.jpeg",sep="")
    
    if(newdevice==TRUE){
      png(heatmap_file,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      #png("test_png.png",res=600,width=8,height=8,units="in",type="cairo")
    }
    
    if(hca_type=="two-way"){
      
      
      
      
      if(col_metabs==TRUE){
        
        col_vec2 <- rainbow(length(unique(mycl_metabs)), alpha=alphacol)
        
        
        if(sample.col.opt=="topo"){
          
          col_vec2 <- topo.colors(length(unique(mycl_metabs)), alpha=alphacol)
        }else{
          if(sample.col.opt=="heat"){
            
            col_vec2 <- heat.colors(length(unique(mycl_metabs)), alpha=alphacol)
          }else{
            if(sample.col.opt=="rainbow"){
              col_vec2<-rainbow(length(unique(mycl_metabs)), start = 0, end = alphacol)
              
              
            }else{
              
              if(sample.col.opt=="terrain"){
                
                col_vec2 <- cm.colors(length(unique(mycl_metabs)), alpha=alphacol)
              }else{
                
                if(sample.col.opt=="colorblind"){
                  #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                  # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                  
                  if(length(unique(mycl_metabs))<9){
                    
                    col_vec2 <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "#E64B35B2", "grey57")
                    
                  }else{
                    
                    #col_vec2<-colorRampPalette(brewer.pal(10, "RdBu"))(length(unique(mycl_metabs)))
                    col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                    
                  }
                  
                  
                }else{
                  
                  check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                  
                  if(length(check_brewer)>0){
                    
                    sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                    col_vec2 <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(unique(mycl_metabs)))
                    
                  }else{
                    
                    if(sample.col.opt=="journal"){
                      
                      col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                 "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                 "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                 
                                 "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                 "#E64B3519","#4DBBD519","#631879E5","grey75")
                      if(length(class_labels_levels)<8){
                        col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75")
                        
                        #col_vec2<-brewer.pal(n = 8, name = "Dark2")
                        
                      }else{
                        if(length(class_labels_levels)<=28){
                          # col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "grey75","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666","#1B9E77", "#7570B3", "#E7298A", "#A6761D", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
                          
                          col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                     "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                     "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                     
                                     "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF", "#8BD76BFF",
                                     "#E64B3519","#9DBBD0FF","#631879E5","#666666","grey75")
                          
                        }else{
                          
                          
                          
                          
                          colfunc <-colorRampPalette(c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75"));col_vec<-colfunc(length(class_labels_levels))
                          
                          col_vec2<-col_vec[sample(col_vec)]
                          
                          
                        }
                      }
                    }else{
                      col_vec2 <-sample.col.opt
                    }
                    
                  }
                  
                }
              }
              
              
            }
            
          }
          
        }
        
        #col_vec2 <- topo.colors(length(unique(mycl_metabs)), alpha=alphacol)
        
        col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                    "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                    "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                    
                    "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                    "#E64B3519","#4DBBD519","#631879E5","grey75")
        
       
        col_vec2<-brewer.pal(10,"Set2")
        colfunc <-colorRampPalette(c(col_vec2))
        
        col_vec2<-colfunc(length(unique(mycl_metabs)))
        if(min(as.numeric(mycl_metabs),na.rm=TRUE)==0){
          
          rowcolors=col_vec2[as.numeric(mycl_metabs)+1] #+1]
        }else{
          rowcolors=col_vec2[as.numeric(mycl_metabs)] #+1]
        }
        
      }else{
        rowcolors=NA #rep("",length(mycl_metabs))
      }
      
      #save(mycl_metabs,rowcolors,file="d2.Rda")
      if(plot.bycluster==TRUE){
        
        #col_vec2 <- rainbow(length(unique(mycl_samples)), alpha=alphacol)
        if(sample.col.opt=="topo"){
          
          col_vec2 <- topo.colors(length(unique(mycl_samples)), alpha=alphacol)
        }else{
          if(sample.col.opt=="heat"){
            
            col_vec2 <- heat.colors(length(unique(mycl_samples)), alpha=alphacol)
          }else{
            if(sample.col.opt=="rainbow"){
              col_vec2<-rainbow(length(unique(mycl_samples)), start = 0, end = alphacol)
              
              
            }else{
              
              if(sample.col.opt=="terrain"){
                
                col_vec2 <- cm.colors(length(unique(mycl_samples)), alpha=alphacol)
              }else{
                
                if(sample.col.opt=="colorblind"){
                  #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                  # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                  
                  if(length(unique(mycl_samples))<9){
                    
                    col_vec2 <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "#E64B35B2", "grey57")
                    
                  }else{
                    
                    # col_vec2<-colorRampPalette(brewer.pal(10, "RdBu"))(length(unique(mycl_samples)))
                    col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                    
                  }
                  
                  
                }else{
                  
                  check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                  
                  if(length(check_brewer)>0){
                    
                    sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                    
                    col_vec2 <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(unique(mycl_samples)))
                    
                  }else{
                    
                    if(sample.col.opt=="journal"){
                      
                      col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                  "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                  "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                  
                                  "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                  "#E64B3519","#4DBBD519","#631879E5","grey75")
                      if(length(class_labels_levels)<8){
                        col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75")
                        
                        #col_vec2<-brewer.pal(n = 8, name = "Dark2")
                        
                      }else{
                        if(length(class_labels_levels)<=28){
                          # col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "grey75","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666","#1B9E77", "#7570B3", "#E7298A", "#A6761D", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
                          
                          col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                      "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                      "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                      
                                      "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF", "#8BD76BFF",
                                      "#E64B3519","#9DBBD0FF","#631879E5","#666666","grey75")
                          
                        }else{
                          
                          
                          
                          
                          colfunc <-colorRampPalette(c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75"));col_vec2<-colfunc(length(class_labels_levels))
                          
                          col_vec2<-col_vec2[sample(col_vec2)]
                          
                          
                        }
                      }
                      
                    }else{
                      col_vec2 <-sample.col.opt
                    }
                    
                  }
                  
                }
              }
              
              
            }
            
          }
          
        }
        
        patientcolors<-col_vec2[as.numeric(mycl_samples)]
        
      }
      
 # save(list=c("data_m","hr","labRow.value","labCol.value","cexLegend","hc","heatmap_cols","mainlab1","rowcolors","patientcolors","cexRow","cexCol","col_vec","class_labels_levels","labRow.value","labCol.value"),file="hcadebug.Rda")

  
      if(col_samples==FALSE){
        if(is.data.znorm==FALSE){
          
          if(is.na(rowcolors)==FALSE){
            
            w <- 0.1
            par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
            
            
            
            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, 
                           symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab=xlab_text,ylab=ylab_text,
                           main=mainlab1,RowSideColors=rowcolors,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
            
            
            if(hca.show.legend==TRUE){
            le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend))
            }
            
            
          }else{
            w <- 0.1
            par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
            
            
            
            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",
                           key=TRUE, symkey=FALSE, density.info="none", trace="none", 
                           cexRow=cexRow, cexCol=cexCol,xlab=xlab_text,ylab=ylab_text, main=mainlab1,
                           labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
            
            if(hca.show.legend==TRUE){
              le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend))
            }
            
          }
          
        }else{
          
          if(is.na(rowcolors)==FALSE){
            w <- 0.1
            par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
            
            
            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols,
                           scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none",
                           cexRow=cexRow, cexCol=cexCol,xlab=xlab_text,ylab=ylab_text, main=mainlab1,
                           RowSideColors=rowcolors,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
            
           
            if(hca.show.legend==TRUE){
             le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend))
            }
            
          }else{
            w <- 0.1
            par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
            
            
            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",
                           key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                           xlab=xlab_text,ylab=ylab_text, main=mainlab1,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
            
            if(hca.show.legend==TRUE){
              le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend))
            }
            
          }
        }
        
      }else{
        
        
        if(is.data.znorm==FALSE){
          
          if(is.na(rowcolors)==FALSE){
            w <- 0.1
            par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
            #print(rowcolors)
            
            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, 
                           symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol, xlab=xlab_text,ylab=ylab_text,
                           main=mainlab1, ColSideColors=patientcolors,RowSideColors=rowcolors,labRow = labRow.value, labCol = labCol.value,
                           cex.main=0.8)
            
            if(hca.show.legend==TRUE){
            (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
            }
            
          }else{
            
            w <- 0.1
            par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
            
           #print("DOING THIS HCA")
            
          #  save(data_m,hr,hc,heatmap_cols,cexRow,cexCol,mainlab1,patientcolors,labRow.value,labCol.value,file="thishca.Rda")
            
            
            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, 
                           density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol, xlab=xlab_text,ylab=ylab_text,
                           main=mainlab1, ColSideColors=patientcolors,labRow = labRow.value, labCol = labCol.value)
            
            # h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol, xlab="Samples",ylab="", main=mainlab1, ColSideColors=patientcolors) #,labRow = FALSE, labCol = FALSE)
            
            if(hca.show.legend==TRUE){
            (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
            }
            
          }
          
        }else{
          
          if(is.na(rowcolors)==FALSE){
            
            w <- 0.1
            par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
            
            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",
                           key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                           xlab=xlab_text,ylab=ylab_text, main=mainlab1, ColSideColors=patientcolors,
                           RowSideColors=rowcolors,labRow = labRow.value, labCol = labCol.value)
            
            if(hca.show.legend==TRUE){
            (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
            }
            
          }else{
            
            w <- 0.1
            par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
            
            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",
                           key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                           xlab=xlab_text,ylab=ylab_text, main=mainlab1, ColSideColors=patientcolors,
                           labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
            
            if(hca.show.legend==TRUE){
            (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
            }
            
          }
        }
      }
      
      
      if(newdevice==TRUE){
        try(dev.off(),silent=TRUE)
      }
      
      
      
      ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
    }
    else{
      if(hca_type=="one-way" | hca_type=="oneway_features" | hca_type=="features"){
        hc<-seq(1,dim(data_m)[2])
        
        
        # mycl_metabs <- cutree(hr, h=max(hr$height)/2)
        
        if(col_metabs==TRUE){
          
          #col_vec2 <- rainbow(length(unique(mycl_metabs)), alpha=alphacol)
          
          if(sample.col.opt=="topo"){
            
            col_vec2 <- topo.colors(length(unique(mycl_metabs)), alpha=alphacol)
          }else{
            if(sample.col.opt=="heat"){
              
              col_vec2 <- heat.colors(length(unique(mycl_metabs)), alpha=alphacol)
            }else{
              if(sample.col.opt=="rainbow"){
                col_vec2<-rainbow(length(unique(mycl_metabs)), start = 0, end = alphacol)
                
                
              }else{
                
                if(sample.col.opt=="terrain"){
                  
                  col_vec2 <- cm.colors(length(unique(mycl_metabs)), alpha=alphacol)
                }else{
                  
                  if(sample.col.opt=="colorblind"){
                    #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                    # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                    
                    if(length(unique(mycl_metabs))<9){
                      
                      col_vec2 <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "#E64B35B2", "grey57")
                      
                    }else{
                      
                      #col_vec2<-colorRampPalette(brewer.pal(10, "RdBu"))(length(unique(mycl_metabs)))
                      col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                  "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                      
                    }
                    
                    
                  }else{
                    
                    check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                    
                    if(length(check_brewer)>0){
                      
                      sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                      col_vec2 <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(unique(mycl_metabs)))
                      
                    }else{
                      
                      if(sample.col.opt=="journal"){
                        
                        col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                    "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                    "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                    
                                    "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                    "#E64B3519","#4DBBD519","#631879E5","grey75")
                        if(length(class_labels_levels)<8){
                          col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75")
                          
                          #col_vec2<-brewer.pal(n = 8, name = "Dark2")
                          
                        }else{
                          if(length(class_labels_levels)<=28){
                            # col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7", "grey75","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666","#1B9E77", "#7570B3", "#E7298A", "#A6761D", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
                            
                            col_vec2<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                        "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                        "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                        
                                        "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF", "#8BD76BFF",
                                        "#E64B3519","#9DBBD0FF","#631879E5","#666666","grey75")
                            
                          }else{
                            
                            
                            
                            
                            colfunc <-colorRampPalette(c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","grey75"));col_vec2<-colfunc(length(class_labels_levels))
                            
                            col_vec2<-col_vec2[sample(col_vec2)]
                            
                            
                          }
                        }
                        
                      }else{
                        col_vec2 <-sample.col.opt
                      }
                      
                    }
                    
                  }
                }
                
                
              }
              
            }
            
          }
          
          
          rowcolors=col_vec2[as.numeric(mycl_metabs)+1]
        }else{
          rowcolors=NA
        }
        
        
        
        if(col_samples==FALSE){
          if(is.data.znorm==FALSE){
            
            if(is.na(rowcolors)==FALSE){
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, 
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,RowSideColors=rowcolors,
                             labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }else{
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE,
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }
          }else{
            
            if(is.na(rowcolors)==FALSE){
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",
                             key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, 
                             cexCol=cexCol,xlab=xlab_text,ylab=ylab_text, main=mainlab1,RowSideColors=rowcolors,
                             labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }else{
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE,
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }
          }
          
        }else{
          if(is.data.znorm==FALSE){
            
            if(is.na(rowcolors)==FALSE){
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE,
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,dendrogram = c("row"),
                             ColSideColors=patientcolors,RowSideColors=rowcolors,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }else{
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, 
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,dendrogram = c("row"),
                             ColSideColors=patientcolors,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }
          }else{
            
            if(is.na(rowcolors)==FALSE){
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, 
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,dendrogram = c("row"),
                             ColSideColors=patientcolors,RowSideColors=rowcolors,labRow = labRow.value, labCol = labCol.value,
                             cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }else{
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE,
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,dendrogram = c("row"),
                             ColSideColors=patientcolors,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }
          }
        }
        
        if(newdevice==TRUE){
          try(dev.off(),silent=TRUE)
        }
        
        mycl_samples<-seq(1,dim(data_m)[2])
        #mycl_samples <- cutree(hc, h=max(hc$height)/2)
        #mycl_metabs <- cutree(hr, h=max(hr$height)/2)
        
        ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
        
      }else{
        if(hca_type=="samples" | hca_type=="oneway_samples"){
          hr<-seq(1,dim(data_m)[1])
          
          
          
          
          if(col_samples==FALSE){
            if(is.data.znorm==FALSE){
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",
                             key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, 
                             cexCol=cexCol,xlab=xlab_text,ylab=ylab_text, main=mainlab1,
                             labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
              
            }else{
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE,
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol, 
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }
            
          }else{
            if(is.data.znorm==FALSE){
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              
              h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE,
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab=xlab_text,ylab=ylab_text,
                             main=mainlab1,dendrogram = c("col"),ColSideColors=patientcolors,
                             labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
            }else{
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              
              h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE,
                             symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,dendrogram = c("col"),
                             ColSideColors=patientcolors,labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }
          }
          
          
          
          if(newdevice==TRUE){
            try(dev.off(),silent=TRUE)
          }
          
          #mycl_samples<-seq(1,dim(data_m)[2])
          # mycl_samples <- cutree(hc, h=max(hc$height)/2)
          mycl_metabs <- seq(1,dim(data_m)[1]) #cutree(hr, h=max(hr$height)/2)
          
          ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
          
        }else{
          
          
          #only plot heatmap; no clustering
          hr<-seq(1,dim(data_m)[1])
          hc<-seq(1,dim(data_m)[2])
          
          
          if(col_samples==FALSE){
            if(is.data.znorm==FALSE){
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, 
                             density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab=xlab_text,ylab=ylab_text,
                             main=mainlab1,dendrogram = c("none"),labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }else{
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE,
                             density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,
                             xlab=xlab_text,ylab=ylab_text, main=mainlab1,dendrogram = c("none"),
                             labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }
            
          }else{
            if(is.data.znorm==FALSE){
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, 
                             density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab=xlab_text,ylab=ylab_text,
                             main=mainlab1,ColSideColors=patientcolors,dendrogram = c("none"),
                             labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
              
            }else{
              
              w <- 0.1
              par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
              
              h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE,
                             density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab=xlab_text,ylab=ylab_text,main=mainlab1,ColSideColors=patientcolors,dendrogram = c("none"),
                             labRow = labRow.value, labCol = labCol.value,cex.main=0.8)
              if(hca.show.legend==TRUE){
              (le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.8, title = "Class",cex=cexLegend)))
              }
            }
          }
          
          #  par(xpd=TRUE)
          #legend("bottomleft",legend=levels(ordered_labels),text.col=unique(patientcolors),pch=13,cex=0.4)
          #par(xpd=FALSE)
          
          if(newdevice==TRUE){
            try(dev.off(),silent=TRUE)
          }
          
          mycl_samples<-seq(1,dim(data_m)[2])
          #mycl_samples <- cutree(hc, h=max(hc$height)/2)
          mycl_metabs <- seq(1,dim(data_m)[1]) #cutree(hr, h=max(hr$height)/2)
          
          ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
          
          
        }
      }
    }
    
    cnames1<-colnames(ord_data)
    cnames1[1]<-"mz_cluster_label"
    colnames(ord_data)<-cnames1
    fname1<-paste("Tables/Clustering_based_sorted_intensity_using_",mainlab,"features.txt",sep="")
    write.table(ord_data,file=fname1,sep="\t",row.names=FALSE)
    
    fname2<-paste("Tables/Sample_clusterlabels_using_",mainlab,"features.txt",sep="")
    
    sample_clust_num<-mycl_samples[h73$colInd]
    temp1<-classlabels[h73$colInd,1]
    temp2<-classlabels
    
    #print(head(temp2))
    #temp1<-as.data.frame(temp1)
    
    #print(dim(temp1))
    match_ind<-match(temp1,temp2[,1])
    
    temp3<-temp2[match_ind,]
    
    #print(head(temp3))
    temp4<-cbind(temp1,temp3,sample_clust_num)
    
    #write.table(temp4,file="s1.txt",sep="\t",row.names=FALSE)
    #   print(head(temp1))
    
    rnames1<-rownames(temp4)
    #temp4<-cbind(rnames1,temp4)
    temp4<-as.data.frame(temp4)
    temp4<-temp4[,-c(1)]
    # print(temp4[,1:4])
    
    
    
    
    if(analysismode=="regression"){
      
      
      
      temp3<-temp4 #[,-c(1)]
      temp3<-as.data.frame(temp3)
      temp3<-apply(temp3,2,as.numeric)
      
      
      temp_vec<-as.vector(temp3[,2])
      
      
      
      names(temp_vec)<-as.character(temp4[,1])
      
      
      # if(output.device.type!="pdf"){
      
      if(newdevice==TRUE){
        
        temp_filename_1<-"Figures/Barplot_dependent_variable_ordered_by_HCA.png"
        
        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      }
      
      
      #print(temp_vec)
      #tiff("Barplot_sample_cluster_ymat.tiff", width=plots.width,height=plots.height,res=plots.res, compression="lzw")
      barplot(temp_vec,col="brown",ylab="Y",cex.axis=0.5,cex.names=0.5,main="Dependent variable levels in samples; \n ordered based on hierarchical clustering")
      #dev.off()
      
      
      
      if(newdevice==TRUE){
        try(dev.off(),silent=TRUE)
      }
      
    }
    
    write.table(temp4,file=fname2,sep="\t",row.names=FALSE)
    
    
    
    fname3<-paste("Tables/Metabolite_clusterlabels_for_",mainlab,"features.txt",sep="")
    
    mycl_metabs_ord<-mycl_metabs[rev(h73$rowInd)]
    write.table(mycl_metabs_ord,file=fname3,sep="\t",row.names=TRUE)
    
  }
  return(list("h73"=h73,"classgroup"=classgroup,"Silhouette.sample"=s1.samp,"Silhouette.feature"=s2.metab,"adj.rand.index"=ari_val,"metab.clusters"=mycl_metabs,"metab.samples"=mycl_samples))
}

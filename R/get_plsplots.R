get_plsplots <-
function(X,plsres,plsvar,samplelabels,filename=NA,ncomp=5,center=TRUE,scale=TRUE,legendcex=0.5,
                       outloc=getwd(),col_vec=NA,sample.col.opt="default",alphacol=0.3,legendlocation="topright",
                       class_levels=NA,pca.cex.val=0.8,pls.ellipse=TRUE,ellipse.conf.level=0.95,main_text="PLS-DA score plots",alphabetical.order=FALSE){
  
  result<-plsres
  r1<-plsvar
  
  pch.val=19
  legendlocation="bottomleft"
  
  samplelabels<-as.data.frame(samplelabels)
  samplelabels<-as.factor(samplelabels[,1])
  if(alphabetical.order==FALSE){
    
    samplelabels=factor(samplelabels,levels=unique(samplelabels))
  }
  l2<-levels(as.factor(samplelabels))
  col_all=topo.colors(256)
  
  t1<-table(samplelabels)
  t1=t1[which(t1>0)]
  
  if(is.na(class_levels)==TRUE){
    
    l1<-levels(as.factor(samplelabels))
  }else{
    l1<-class_levels
    
    
  }
  
  class_labels_levels<-l1
  
  ncomp=min(dim(X)[1],dim(X)[2])
  
  if(is.na(col_vec)==TRUE){
    col_vec<-c("mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
               "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
               "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
               "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
    
  }
  
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
                
                #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                col_vec<-c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                           "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                
              }
              
              
            }else{
              
              check_brewer<-grep(pattern="brewer",x=sample.col.opt)
              
              if(length(check_brewer)>0){
                
                sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                
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
                  
                 # if(length(sample.col.opt)==1){
                  #  col_vec <-rep(sample.col.opt,length(class_labels_levels))
                  #}else{
                    
                   # colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                    
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
  #col_vec<-col_vec[sample(1:length(col_vec),length(col_vec))]
  
  l1<-gsub(x=l1,pattern="Class",replacement="")
  
  dir.create(outloc,showWarnings=FALSE)
  setwd(outloc)
  
  ## 1) raw data
  #tiff(fname,res=300, width=2000,height=2000)
  col <- rep(col_vec[1:length(t1)], t1)
  #col<-rep(col_all[1:length(l1)],t1)
  ## Choose different size of points
  cex <- rep(pca.cex.val, dim(X)[1])
  ## Choose the form of the points (square, circle, triangle and diamond-shaped
  
  
  
  #  pch_vec=c(3,5,7,9,12,13,2,17,21,23) 
  #}else{
   # set.seed(999)
  pch_vec<-rep(seq(0,25),2) #sample(seq(1,50),50) #c(3,5,7,9,12,13,2,17,21) # #
  #}
  #pch_vec <- rep(pch.val,dim(X)[1])
  #pch=pch_vec
  
  pch <- rep(pch_vec[1:length(t1)], t1)
  
  cex <- rep(pca.cex.val, dim(X)[1])
  col_per_group<-{}
  pch_per_group<-{}
  for(p1 in 1:length(l2)){
    
    pch[which(samplelabels==l2[p1])]=pch_vec[p1]
    col[which(samplelabels==l2[p1])]=col_vec[p1]
    col_per_group<-c(col_per_group,col_vec[p1])
    pch_per_group<-c(pch_per_group,pch_vec[p1])
  }
  
  
  #pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
  #pch_vec <- rep(pch.val,dim(X)[1])
  cex <- rep(pca.cex.val, dim(X)[1])
  
  
  ####savelist=ls(),file="plsdebug.Rda")
  
  #print(plotIndiv(result, comp = c(1,2),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=FALSE, ellipse.level = 0.95, X.label=paste("PLS1 (",r1[1],"% variation)",sep=""),Y.label=paste("PLS2 (",r1[2],"% variation)",sep=""),add.legend=TRUE))
  
  #save(result,samplelabels,col_per_group,col_vec,main_text,col,cex,pch,pls.ellipse,r1,file="plsdebug.Rda")
  if(result$ncomp>1){
    
   # w <- 0.1
    #par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
    
   # print(plotIndiv(result, comp = c(1,2),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=pls.ellipse, style="lattice",col.per.group=col_per_group,
    #                ellipse.level = 0.95, X.label=paste("PLS1 (",r1[1],"% variation)",sep=""),Y.label=paste("PLS2 (",r1[2],"% variation)",sep=""),legend=TRUE,title=main_text))
    plotIndiv(result,comp=c(1,2),group=as.factor(samplelabels),legend = TRUE,
              ellipse=pls.ellipse,style="lattice",
              title = main_text,col.per.group=col_per_group,
              X.label=paste("PLS1 (",r1[1],"% variation)",sep=""),Y.label=paste("PLS2 (",r1[2],"% variation)",sep=""))
    
   
  }
  if(result$ncomp>2){
    
   # w <- 0.1
    #par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
   # print(plotIndiv(result, comp = c(1,3),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=pls.ellipse, style="lattice",col.per.group=col_per_group,
    #                ellipse.level = 0.95, X.label=paste("PLS1 (",r1[1],"% variation)",sep=""),Y.label=paste("PLS3 (",r1[3],"% variation)",sep=""),legend=TRUE,title=main_text))
    plotIndiv(result,comp=c(1,3),group=as.factor(samplelabels),legend = TRUE,
              ellipse=pls.ellipse,style="lattice",
              title = main_text,col.per.group=col_per_group,
              X.label=paste("PLS1 (",r1[1],"% variation)",sep=""),Y.label=paste("PLS3 (",r1[3],"% variation)",sep=""))
    
   
   # w <- 0.1
    #par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
   
   # print(plotIndiv(result, comp = c(2,3),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=pls.ellipse, style="lattice",col.per.group=col_per_group,
    #                ellipse.level = 0.95, X.label=paste("PLS2 (",r1[2],"% variation)",sep=""),Y.label=paste("PLS3 (",r1[3],"% variation)",sep=""),legend=TRUE,title=main_text))
    
    plotIndiv(result,comp=c(2,3),group=as.factor(samplelabels),legend = TRUE,
              ellipse=pls.ellipse,style="lattice",
              title = main_text,col.per.group=col_per_group,
              X.label=paste("PLS2 (",r1[2],"% variation)",sep=""),Y.label=paste("PLS3 (",r1[3],"% variation)",sep=""))
    
  }
  
  
  
}

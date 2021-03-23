get_individualsampleplots <-
function(feature_table_file,class_labels_file,X=NA,Y=NA,parentoutput_dir,newdevice=FALSE,
                                    ylabel="Intensity",bar.colors=NA,cex.plots=0.75,sample.col.opt=c("grey57"),alphabetical.order=FALSE,name=NA){
  
  cex.val=cex.plots
  
  plottype="barplot"
  
  if(typeof(X)=="logical"){
    data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  }else{
    data_matrix<-X
    rm(X)
    
  }
  
  par(mfrow=c(1,1),family="sans",cex=cex.val)
  alphacol=0.3
  
  dir.create(parentoutput_dir,showWarnings = FALSE)
  setwd(parentoutput_dir)
  
  mzlabels<-data_matrix[,1]
  
  timelabels<-data_matrix[,2]
  
  
  
  data_m<-data_matrix[,-c(1:2)]
  
  data_m<-as.matrix(data_m)
  
  col_samples<-TRUE
  
  
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
  
  
  
  classlabels<-as.data.frame(classlabels)
  
  
  
  
  if(dim(classlabels)[2]>2){
    
    Class<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
  }else{
    
    Class<-classlabels[,2]
  }
  
  #keeps the class order same as in the input file; avoids arrangement by alphabetical order
  if(alphabetical.order==FALSE){
    Class <- factor(Class, levels=unique(Class))
  }
  
  class_labels_levels<-levels(as.factor(Class))
  
  
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
                  #col_vec <-rep(sample.col.opt,length(class_labels_levels))
                #  if(length(sample.col.opt)==1){
                 #   col_vec <-rep(sample.col.opt,length(class_labels_levels))
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
  
  mtcars<-t(data_m)
  Class<-as.character(Class)
  
  barplot.col.opt=col_vec
  
  mtcars<-cbind(Class,mtcars)
  mtcars<-as.data.frame(mtcars)
  
  mtcars<-mtcars[order(mtcars$Class),]
  
  if(newdevice==TRUE){
    
    pdf("individual.sample.plots.pdf")
    
  }
  
  col_per_group<-{}
  pch_per_group<-{}
  
  l2<-levels(as.factor(Class))
  col1<-rep("red",length(Class))
  
  
  
  par_rows=2
  max_per_row=2
  
  # par(mfrow=c(par_rows,max_per_row))
  
  myData_list<-new("list")
  
  mtcars <- do.call(data.frame, mtcars)
  
  ####savemtcars,file="mtcars.Rda")
  
  
  
  #for(i in 2:dim(mtcars)[2]){
  lapply(2:dim(mtcars)[2],function(i){
    x<-as.numeric(as.character(mtcars[,i]))
    
    myData <- x
    
    myData<-cbind(mtcars$Class,myData)
    colnames(myData)<-c("Class","Intensity")
    #myData<-myData[c(2:3,1,4:5),]
    
    myData <- as.data.frame(myData)
    myData<-myData[order(myData$Class),]
    
    
    
    round_mzval<-mzlabels[(i-1)] #sprintf("%.4f",mzlabels[(i-1)])
    round_timeval<-timelabels[(i-1)] #sprintf("%.1f",timelabels[(i-1)])
    
    #mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
    
    if(is.na(name[1])==TRUE){
      
      mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
    }else{
      
      mzlabel_cur<-name[i-1]
    }
    
    t1<-table(myData$Class)
    if(length(barplot.col.opt)<2){
      
      barplot.col.opt1=rep(barplot.col.opt,length(myData$Class))
    }else{
      
      
      if(length(barplot.col.opt)==length(levels(factor(myData$Class)))){
        
        #barplot.col.opt1=rep(barplot.col.opt,t1)
        
        barplot.col.opt1=rep(barplot.col.opt[1:length(t1)],t1)
        
        
      }else{
        
        # print("Number of classes is greater than the length of the color vector. Using default colors.")
        #col_clust<-topo.colors(length(t1))
        #barplot.col.opt1=col1
        #t1<-table(mtcars$Class)
        #barplot.col.opt1=rep(barplot.col.opt,t1)
        barplot.col.opt1=rep(barplot.col.opt[1:length(t1)],t1)
      }
      
      
    }
    
    #print(barplot.col.opt1)
    #barplot.col.opt1=barplot.col.opt
    
    # print(barplot.col.opt1)
    
    if(alphabetical.order==FALSE){
      myData$Class<-factor(myData$Class,levels=unique(myData$Class))
      
    }else{
      myData$Class<-factor(myData$Class)
    }
    Samples<-seq(1,nrow(myData))
    Var2<-rep(1,nrow(myData))
    myData<-cbind(myData,Samples,Var2)
    
    
  #  save(myData,file="myData.Rda")
    # ###savelist=ls(),file="barplotdebug.Rda")
    
    if(plottype=="barplot"){
      #barplot(as.vector(myData[,2]),col=c(barplot.col.opt1),main=mzlabel_cur, ylab="Intensity",xlab="Sample") #,ylim=c(min(myData[,2])-1,max(myData[,2])+1),xpd=FALSE)
      
      
      w <- 0.1
      par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
      
      min_yval=min(myData[,2],na.rm=TRUE)
      ymax=max(myData[,2],na.rm=TRUE)
      ylim=range(pretty(c(min_yval,ymax)))
      
      plot(as.vector(myData[,2]),col=c(barplot.col.opt1),main=mzlabel_cur, ylab=ylabel,xlab="Sample",type="h",lwd=2,ylim=ylim)
      
      
      
      legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,levels(as.factor(mtcars$Class)), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])), pt.cex = 0.6, title = "Class",cex=0.8)
      
      
      if(FALSE){
        barplot1 <- ggplot(data=myData) + geom_bar(aes(x=Samples, y=Intensity, fill=Class), stat="identity") + facet_wrap(~Var2) + scale_fill_manual(values=unique(barplot.col.opt1)) + ggtitle(mzlabel_cur) + ylab(ylabel) + xlab("Samples") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=sizeval), plot.margin=unit(c(10,5,5,5),"mm"),
                                                                                                                                                                                                                                                                   strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                                                                                                                                                                                   strip.text = element_text(face="bold"))
      }
      
      #print(barplot1)
      
    }else{
      if(plottype=="point"){
        plot(as.vector(myData[,2]),col=c(barplot.col.opt1),main=mzlabel_cur, ylab=ylab,xlab="Samples",type="p",ylim=range(pretty(c(0,max(myData[,2],na.rm=TRUE)))))
      }else{
        
        if(plottype=="ggplot"){
          p <- ggplot(data = myData, aes(x = Class, y = Intensity,fill=Class)) + scale_fill_manual(values = c(barplot.col.opt1)) + ylab(ylabel)
          p=p+geom_point(size=cex.plots)+ggtitle(mzlabel_cur) #,aes(colour=barplot.col.opt))
          print(p)
        }
      }
    }
    
    
  })
  if(newdevice==TRUE){
    try(dev.off(),silent=TRUE)
  }
  
  
  
  
}

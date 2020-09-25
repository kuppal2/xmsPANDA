get_barplots <-
function(feature_table_file,class_labels_file,X=NA,Y=NA,parentoutput_dir,newdevice=FALSE,ylabel="Intensity",bar.colors=NA,cex.plots=0.6,barplot.col.opt="journal",error.bar=TRUE,barplot.xaxis="Factor2",alphabetical.order=FALSE,name=NA,study.design="oneway"){
    
    analysistype=study.design
    cex.val=cex.plots
    xaxis=barplot.xaxis
    if(is.na(X[1])==TRUE){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
    }else{
        data_matrix<-X
        rm(X)
        
    }
    
    dir.create(parentoutput_dir)
    setwd(parentoutput_dir)
    
    
    mzlabels<-data_matrix[,1]
    
    timelabels<-data_matrix[,2]

    data_m<-data_matrix[,-c(1:2)]
    
    data_m<-as.matrix(data_m)
    
    col_samples<-TRUE
    
   
    
    if(is.na(Y)==TRUE){
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
    
    pairedanalysis=FALSE
    if(analysistype=="twowayrepeat" | analysistype=="2wayrepeat" | analysistype=="onewayrepeat" | analysistype=="1wayrepeat"){
                    
                    pairedanalysis=TRUE
                }
           
    
    if(dim(classlabels)[2]>2){
        if(pairedanalysis==TRUE){
            
            #    paireddesign=classlabels[,2]
            
            #classlabels_orig<-classlabels_orig[,-c(2)]
        }
    
        if(analysistype=="twowayrepeat" | analysistype=="2wayrepeat" | analysistype=="twoway" | analysistype=="2way"){
            Class<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
            
            tabbed_groups<-TRUE
        }else{
            Class<-classlabels[,2]
            
            tabbed_groups=FALSE
        }
    
    }else{
        
        Class<-classlabels[,2]
        
        tabbed_groups=FALSE
    }
    
    if(alphabetical.order==FALSE){
        Class <- factor(Class, levels=unique(Class))
    }

    
    class_labels_levels<-levels(as.factor(Class))
    
    
    if(is.na(barplot.col.opt)==FALSE)
    {
        if(barplot.col.opt=="default"){
            
            col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
            "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
            "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
            "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
            
        }else{
            if(barplot.col.opt=="topo"){
                #col_vec<-topo.colors(256) #length(class_labels_levels))
                
                #col_vec<-col_vec[seq(1,length(col_vec),)]
                
                col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
            }else{
                if(barplot.col.opt=="heat"){
                    #col_vec<-heat.colors(256) #length(class_labels_levels))
                    
                    col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
                }else{
                    if(barplot.col.opt=="rainbow"){
                        #col_vec<-heat.colors(256) #length(class_labels_levels))
                        col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
                        
                        #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
                    }else{
                        
                        if(barplot.col.opt=="terrain"){
                            
                            col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
                        }else{
                            if(is.na(barplot.col.opt)==TRUE){
                                col_vec<-c("black")
                            }else{
                                
                                if(barplot.col.opt=="colorblind"){
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
                                    
                                    check_brewer<-grep(pattern="brewer",x=barplot.col.opt)
                                    
                                    if(length(check_brewer)>0){
                                        barplot.col.opt=gsub(x=barplot.col.opt,pattern="brewer.",replacement="")
                                        col_vec <- colorRampPalette(brewer.pal(10, barplot.col.opt))(length(class_labels_levels))
                                        
                                    }else{
                                        
                                        if(barplot.col.opt=="journal"){
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
                                            col_vec <-barplot.col.opt
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
        
        col_vec<-c("grey57")
        
    }

    
    mtcars<-t(data_m)
    Class<-as.character(Class)
    
  
    
    
    
    if(newdevice==TRUE){
    
    pdf("barplots.pdf")
    
    plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")


                           text(5, 8, "Barplots using the\n normalized intensities/adundance levels")
                 
    
    }
    
    #par(mfrow=c(2,2),family="sans",cex=cex.plots)
    #par(mfrow=c(2,2),family="sans",cex=cex.plots,pty="s")
    # par(mfrow=c(par_rows,max_per_row))
    
    myData_list<-new("list")
    
   
    
   time_start<-Sys.time()
   ##savelist=ls(),file="debugbarplot.Rda")
   if(tabbed_groups==FALSE){
       
       mtcars<-cbind(Class,mtcars)
       mtcars<-as.data.frame(mtcars)
       mtcars <- do.call(data.frame, mtcars)
       
       if(alphabetical.order==FALSE){
           mtcars$Class<- factor(mtcars$Class, levels=unique(mtcars$Class))
           
           
       }
        
       mtcars_sum<-do.call(data.frame,aggregate(list(mtcars[-c(1)]),by=list(mtcars$Class),FUN = function(x) {x<-as.numeric(as.character(x));c(mean = mean(x), sd = sd(x),n = length(x),se=sd(x)/sqrt(length(x)))}))
       
       label_inc_list<-seq(2,dim(mtcars_sum)[2],4) #seq(1,dim(mtcars_sum)[2],3)
       # #save(list=ls(),file="debugbarplot.Rda")
       #  for(i in 2:dim(mtcars)[2]){
       myData_list<-lapply(seq(2,dim(mtcars_sum)[2],4),function(i){
           
           myData<-mtcars_sum[,c(1,i:(i+3))]
           colnames(myData) <- c("Class", "Intensity", "sd", "n", "se")
           
           get_label_ind<-which(label_inc_list==i)
           round_mzval<-sprintf("%.4f",mzlabels[get_label_ind])
           round_timeval<-sprintf("%.1f",timelabels[get_label_ind])
           
           if(is.na(name)==TRUE){
       
           mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
           }else{
               
               mzlabel_cur<-name[get_label_ind]
               
           }
           
        
           ymax = myData$Intensity + 1.96*myData$se
           
           ymin = myData$Intensity - 1.96*myData$se
           
           max_yval<-ceiling(max((myData$Intensity + (2.5*myData$se)),na.rm=TRUE)) #round(max(myData$Intensity+(4*myData$se),na.rm=TRUE))
           
           
           min_yval<-max(0,floor(min((myData$Intensity - (2.5*myData$se)),na.rm=TRUE)))
           
           below_zero_check<-which(ymin<0)
           if(length(below_zero_check)>0){
               
               ymin[below_zero_check]<-0
           }
           
           
           
           myData$Class<-as.factor(myData$Class)
           
           colnames(myData) <- c("Class", "Intensity", "sd", "n", "se")
           
           t1<-table(myData$Class)
           
           barplot.col.opt1=rep(col_vec[1:length(t1)],t1)
           
           
           
           w <- 0.1
           par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
           ylim=range(pretty(c(min_yval,ymax)))
           rnum=length(unique(myData$Class))
           if(rnum<7)
           {
               
               space_vec=c(2,2,1.5,1,0.5,0.25)
               name_vec=unique(myData$Class)
              # barCenters=barplot(myData$Intensity,ylab=ylabel,xlab="",main=mzlabel_cur,col=barplot.col.opt1,width=0.125,xlim=c(0,1),names.arg=name_vec[1:rnum],xpd=FALSE,ylim=ylim,space=space_vec[rnum])
              barCenters=barplot(myData$Intensity,ylab=ylabel,xlab="",main=mzlabel_cur,col=barplot.col.opt1,width=0.125,xlim=c(0,1),names.arg=NULL,xpd=FALSE,space=space_vec[rnum], ylim=ylim)
            
           }else{
            #   barCenters=barplot(myData$Intensity,ylab=ylabel,xlab="",main=mzlabel_cur,col=barplot.col.opt1,space=0.1,names.arg=unique(myData$Class),xpd=FALSE,ylim=ylim)
            barCenters=barplot(myData$Intensity,ylab=ylabel,xlab="",main=mzlabel_cur,col=barplot.col.opt1,space=0.1,names.arg=NULL,xpd=FALSE,ylim=ylim)
           }
           
           if(error.bar==TRUE){
               arrows(barCenters, ymin,barCenters, ymax,angle=90,code=3,lty=1,lwd=1.25,length=0.05)
           }
           
           
           
          legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(myData$Class), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])), pt.cex = 0.6, title = "Class",cex=0.8)
           
           return(myData)
           
       })
       
   }else{
       
       mtcars<-cbind(classlabels[,2:3],mtcars)
       mtcars<-as.data.frame(mtcars)
       cnames_mtcars<-colnames(mtcars)
       cnames_mtcars[1:2]<-c("Factor1","Factor2")
       colnames(mtcars)<-cnames_mtcars
       mtcars <- do.call(data.frame, mtcars)
       
       
       if(alphabetical.order==FALSE){
          mtcars$Factor1<- factor(mtcars$Factor1, levels=unique(mtcars$Factor1))
          mtcars$Factor2<- factor(mtcars$Factor2, levels=unique(mtcars$Factor2))
           
       }
       
       mtcars_sum<-do.call(data.frame,aggregate(list(mtcars[-c(1:2)]),by=list(Factor1=mtcars$Factor1,Factor2=mtcars$Factor2),FUN = function(x) {x<-as.numeric(as.character(x));c(mean = mean(x), sd = sd(x),n = length(x),se=sd(x)/sqrt(length(x)))}))
        
      
       
       
       label_inc_list<-seq(3,dim(mtcars_sum)[2],4) #seq(1,dim(mtcars_sum)[2],3)
       ####savelist=ls(),file="debugbarplot.Rda")
       #  for(i in 2:dim(mtcars)[2]){
       myData_list<-lapply(seq(3,dim(mtcars_sum)[2],4),function(i){
           
           myData<-mtcars_sum[,c(1:2,i:(i+3))]
           colnames(myData) <- c("Factor1", "Factor2", "Intensity", "sd", "n", "se")
           
           if(xaxis=="Factor1"){
               tabbedMeans <- tapply(myData$Intensity, list(myData$Factor2,myData$Factor1),
               function(x) c(x = x))
               
               tabbedSE <- tapply(myData$se, list(myData$Factor2,myData$Factor1),
               function(x) c(x = x))
               
           }else{
               
               if(xaxis=="Factor2"){
                   tabbedMeans <- tapply(myData$Intensity, list(myData$Factor1,myData$Factor2),
                   function(x) c(x = x))
                   
                   tabbedSE <- tapply(myData$se, list(myData$Factor1,myData$Factor2),
                   function(x) c(x = x))
               }
               
               
           }
           get_label_ind<-which(label_inc_list==i)
           round_mzval<-sprintf("%.4f",mzlabels[get_label_ind])
           round_timeval<-sprintf("%.1f",timelabels[get_label_ind])
           
           # mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
           
           if(is.na(name)==TRUE){
                
                    mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
                    }else{
                        
                        mzlabel_cur<-name[get_label_ind]
                        
                    }
                    
           
           ymax = tabbedMeans + 1.96*tabbedSE
           
           ymin = tabbedMeans - 1.96*tabbedSE
           
           max_yval<-ceiling(max((tabbedMeans + (2.5*tabbedSE)),na.rm=TRUE)) #round(max(myData$Intensity+(4*myData$se),na.rm=TRUE))
           
           min_yval<-max(0,floor(min((tabbedMeans - (2.5*tabbedSE)),na.rm=TRUE)))
           
           below_zero_check<-which(ymin<0)
           
           if(length(below_zero_check)>0){
               
               ymin[below_zero_check]<-0
           }
           
           if(xaxis=="Factor2"){
               
               t1<-table(factor(myData$Factor1))
               t2<-table(factor(myData$Factor2))
           }else{
               t1<-table(factor(myData$Factor2))
               t2<-table(factor(myData$Factor1))
           }
           
           barplot.col.opt1=rep(col_vec[1:length(t1)],length(t2))
           
           w <- 0.1
           par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
           ylim=range(pretty(c(min_yval,ymax)))
           
          
          barCenters=barplot(tabbedMeans,beside=TRUE, ylab=ylabel,xlab="",main=mzlabel_cur,col=barplot.col.opt1,las=1,names.arg=unique(myData$Class),xpd=FALSE,border = "black",ylim=ylim)
          
          segments(barCenters, ymin, barCenters,ymax, lwd = 1.25)
          
           if(error.bar==TRUE){
               # arrows(barCenters, ymin,barCenters, ymax,angle=90,code=3,lty=1,lwd=1.25,length=0.05)
               arrows(barCenters, ymin, barCenters, ymax, lwd = 1.25, angle = 90, code = 3, length = 0.05)
           }
           
           if(xaxis=="Factor1"){
               legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(myData$Factor2), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])), pt.cex = 0.6, title = "",cex=0.8)
           }else{
           legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(myData$Factor1), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])), pt.cex = 0.6, title = "",cex=0.8)
               
           }
           return(myData)
           
       })
       
       
       
       
       
       
   }
   
    
    
    if(newdevice==TRUE){
        try(dev.off(),silent=TRUE)
    }
    
    # suppressWarnings(dir.create("Tables"))
   
    ###savemyData_list,file="Tables/barplots_data.Rda")

    
}

get_pca <-
function(X,samplelabels,legendlocation="topright",filename=NA,ncomp=5,pcacenter=TRUE,pcascale=TRUE,legendcex=0.5,
                  outloc=getwd(),col_vec=NA,sample.col.opt="default",alphacol=0.3,class_levels=NA,pca.cex.val=3,pca.ellipse=TRUE,
                  ellipse.conf.level=0.95,samplenames=FALSE,do_pca_anova=FALSE,paireddesign=NA,pairedanalysis=FALSE,classlabelsorig=NA,alphabetical.order=FALSE,analysistype="oneway",lme.modeltype="RI"){
    
    suppressMessages(library(mixOmics))
    suppressMessages(library(car))
    
    X<-as.matrix(t(X))
    
    par(mfrow=c(1,1),family="sans",cex=0.9)
    pch.val<-19
 
    samplelabels<-as.data.frame(samplelabels)
    
    samplelabels<-paste("",as.factor(samplelabels[,1]),sep="")
    
    if(alphabetical.order==FALSE){
    samplelabels <- factor(samplelabels, levels=unique(samplelabels))
    }
    
    if(analysistype=="twowayrepeat" | analysistype=="2wayrepeat" | analysistype=="onewayrepeat" | analysistype=="1wayrepeat"){
        
        pairedanalysis=TRUE
        
        
    }
    
    
 
    
    if(dim(classlabelsorig)[2]==2){
        X=X[order(classlabelsorig[,2]),]
        samplelabels<-samplelabels[order(classlabelsorig[,2])]
        classlabelsorig<-classlabelsorig[order(classlabelsorig[,2]),]
    
    }else{
        
        
        if(dim(classlabelsorig)[2]==3){
            
            
            if(analysistype=="twowayrepeat" | analysistype=="2wayrepeat" | analysistype=="2way" | analysistype=="twoway"){
            
                X=X[order(classlabelsorig[,2],classlabelsorig[,3]),]
               samplelabels<-samplelabels[order(classlabelsorig[,2],classlabelsorig[,3])]
               classlabelsorig<-classlabelsorig[order(classlabelsorig[,2],classlabelsorig[,3]),]
               
            }else{
                X=X[order(classlabelsorig[,2]),]
                samplelabels<-samplelabels[order(classlabelsorig[,2])]
                classlabelsorig<-classlabelsorig[order(classlabelsorig[,2]),]
                
            }
            
        }else{
                if(analysistype=="twowayrepeat" | analysistype=="2wayrepeat" | analysistype=="2way" | analysistype=="twoway"){
                
                    X=X[order(classlabelsorig[,2],classlabelsorig[,3]),]
                   samplelabels<-samplelabels[order(classlabelsorig[,2],classlabelsorig[,3])]
                   classlabelsorig<-classlabelsorig[order(classlabelsorig[,2],classlabelsorig[,3]),]
                   
                }else{
                    X=X[order(classlabelsorig[,2]),]
                    samplelabels<-samplelabels[order(classlabelsorig[,2])]
                    classlabelsorig<-classlabelsorig[order(classlabelsorig[,2]),]
                    
                }
        }
    }
    
    if(alphabetical.order==FALSE){
     samplelabels <- factor(samplelabels, levels=unique(samplelabels))
    }
    l2<-levels(as.factor(samplelabels))
    col_all=topo.colors(256)
    
    t1<-table(samplelabels)
    if(is.na(class_levels)==TRUE){
        
        l1<-levels(as.factor(samplelabels))
    }else{
        l1<-class_levels
        
        
    }
    l1<-levels(as.factor(samplelabels))
    ##save(X,file="pcaX.Rda")
    #print(dim(X))
    if(pcascale=="pareto"){
        
        X<-apply(X,2,function(x){y<-(x-mean(x,na.rm=TRUE))/sqrt(sd(x,na.rm=TRUE));return(y)})
       
        pcascale=FALSE
    }else{
        
        if(pcascale=="uv" || pcascale=="autoscale"){
            
            #X<-apply(X,2,function(x){y<-(x-mean(x,na.rm=TRUE))/(sd(x,na.rm=TRUE));return(y)})
            
            pcascale=TRUE
        }
        
    }
    
    ##save(X,file="pcaX1.Rda")
    class_labels_levels<-l1
    
    ncomp=min(c(10,dim(X)[1],dim(X)[2]))
    
    #   p1<-pcaMethods::pca(t(X),method="svd",center=TRUE,scale="uv",cv="q2",nPcs=10)
    if(is.na(paireddesign)==TRUE){
    metabpcaresultlog2allmetabs5pcs<-mixOmics::pca(X,ncomp=ncomp,center=pcacenter,scale=pcascale)
    }else{
        metabpcaresultlog2allmetabs5pcs<-mixOmics::pca(X,ncomp=ncomp,center=pcacenter,scale=pcascale) #,multilevel=paireddesign)
          
    }
    result<-metabpcaresultlog2allmetabs5pcs
    

  fname1<-paste("pcares",filename,".Rda",sep="")
  #save(result,file=fname1)
  if(analysistype=="regression"){
    
  }
  #save(X,result,samplelabels,ellipse.conf.level,classlabelsorig,analysistype,file="debug3.Rda")
 
    
    s1<-summary(result)
    r1<-s1$importance[2,]
    r1<-round(r1,2)*100
    
      barplot(r1,beside=TRUE,main="% variation explained by each component",ylab="% variation",col=c("#0072B2"),cex.main=0.7,ylim=range(pretty(c(0,max(r1)))))
    
    
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
                                    if(length(sample.col.opt)==1){
                                        col_vec <-rep(sample.col.opt,length(class_labels_levels))
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
    #col_vec<-col_vec[sample(1:length(col_vec),length(col_vec))]
    
    l1<-gsub(x=l1,pattern="Class",replacement="")
    
    dir.create(outloc,showWarnings=FALSE)
    setwd(outloc)
    # print(paste("Generating PCA plots",sep=""))
    
    fname<-paste("PCA_eval",filename,".tiff",sep="")
    ## 1) raw data
    #tiff(fname,res=300, width=2000,height=2000)
    col <- rep(col_vec[1:length(t1)], t1)
    
    #col<-rep(col_all[1:length(l1)],t1)
    ## Choose different size of points
    cex <- rep(pca.cex.val, dim(X)[1])
    ## Choose the form of the points (square, circle, triangle and diamond-shaped
    
    pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
    pch <- rep(pch.val,dim(X)[1])
    #pch_vec <- rep(pch.val,dim(X)[1])
    
    
    cex <- pca.cex.val #rep(pca.cex.val, dim(X)[1])
    col_per_group<-{}
    pch_per_group<-{}
    for(p1 in 1:length(l1)){
        col[which(samplelabels==l1[p1])]=col_vec[p1]
        pch[which(samplelabels==l1[p1])]=pch_vec[p1]
        
        col_per_group<-c(col_per_group,col_vec[p1])
        pch_per_group<-c(pch_per_group,pch_vec[p1])
    }
   pch=pch_per_group
   # print(table(pch))
   
   main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing",sep="")
   
       legendcex<-0.7 #0.5*pca.cex.val

       #if(pca.ellipse=="car")
   {
      
      
      pc_pval_single<-{}
     
     #    do_pca_anova=FALSE
      
     
      if(do_pca_anova==TRUE)
      {
          scores_res<-result$x
           
          
          if(dim(classlabelsorig)[2]==2){
              dtemp<-cbind(classlabelsorig,scores_res)
              dtemp<-as.data.frame(dtemp)
              
              pc1_pval<-anova(lm(cbind(scores_res[,1],scores_res[,2])~samplelabels))
              
              pc1_pval<-pc1_pval[[6]][2]
              pc2_pval<-anova(lm(cbind(scores_res[,1],scores_res[,3])~samplelabels))
              
              pc2_pval<-pc2_pval[[6]][2]
              
              pc3_pval<-anova(lm(cbind(scores_res[,2],scores_res[,3])~samplelabels))
              pc3_pval<-pc3_pval[[6]][2]
              ##save(dtemp,samplelabels,classlabelsorig,scores_res,paireddesign,file="pcdtemp.Rda")
              if(is.na(paireddesign)==TRUE){
                  
                  
                 testname="one-way ANOVA"
                  pc_pval_single<-lapply(1:min(5,ncol(scores_res)),function(pcn1){
                      pc1_only<-anova(lm(scores_res[,pcn1]~samplelabels))
                      pc1_only<-pc1_only[[5]][1]
                      return(pc1_only)
                      
                  })
                  
                  pc_pval_single<-unlist(pc_pval_single)
              }else{
                  
                   testname="one-way ANOVA with repeated measures"
                   
                     #one-way ANOVA repeat
                      pc_pval_single<-lapply(3:ncol(dtemp),function(pcn1){
                          dataA<-cbind(dtemp[,pcn1],dtemp[,c(2)])
                          
                          colnames(dataA)<-c("Response","Factor1")
                          
                          #pcalmonewayanova
                          pc1_res<-diffexplmonewayanovarepeat(dataA=dataA,subject_inf=paireddesign,modeltype=lme.modeltype)
                          
                          return(pc1_res$mainpvalues)
                          
                      })
                      
                      pc_pval_single<-unlist(pc_pval_single)
                   
              }
          
          }else{
              
              dtemp<-cbind(classlabelsorig,scores_res)
              dtemp<-as.data.frame(dtemp)
              
              #save(dtemp,classlabelsorig,paireddesign,file="pcdtemp2.Rda")
              
              if(dim(classlabelsorig)[2]>=2){
                  
                  if(is.na(paireddesign)==FALSE){
                      
                      testname="two-way ANOVA with repeated measures in one factor"
                          #two-way ANOVA repeat
                          pc_pval_single<-lapply(4:ncol(dtemp),function(pcn1){
                              dataA<-cbind(dtemp[,pcn1],dtemp[,c(2:3)])
                              
                              colnames(dataA)<-c("Response","Factor1","Factor2")
                              
                              dataA$Response<-dtemp[,pcn1]
                              
                              #save(dataA,paireddesign,file="pcdataA.rda")
                              pc1_res<-diffexplmtwowayanovarepeat(dataA=dataA,subject_inf=paireddesign,modeltype=lme.modeltype)

                               return(pc1_res$mainpvalues)
                              
                          })
                          pc_pval_single<-unlist(pc_pval_single)
                  }else{
                      
                      testname="two-way ANOVA"
                      #two-way ANOVA
                      pc_pval_single<-lapply(4:ncol(dtemp),function(pcn1){
                          dataA<-cbind(dtemp[,pcn1],dtemp[,c(2:3)])
                          
                          colnames(dataA)<-c("Response","Factor1","Factor2")
                          
                          dataA$Response<-dtemp[,pcn1]
                          
                          ###savedataA,file="pcdataA.rda")
                          pc1_res<-diffexplmtwowayanova(dataA=dataA)
                          
                          return(pc1_res$mainpvalues)
                          
                      })
                      pc_pval_single<-unlist(pc_pval_single)
                      
                  }
                  
              }
              
          }
          
          if(dim(classlabelsorig)[2]==2){
          pc_pval_vec<-c(pc1_pval,pc2_pval,pc3_pval)
          
          }
      
        pc_pval_single<-round(pc_pval_single,3)
      }
      
      if(pca.ellipse==TRUE){
      
            suppressMessages(library(car))
      }
      
      
     
     if(dim(classlabelsorig)[2]==2){
         if(do_pca_anova==TRUE){
             main_text=paste("Pairwise PC score plots using ",filename," features\np-value for differences between groups using PC1 and PC2 in a multivariate\n ANOVA model=",round(pc_pval_vec[1],3),sep="")
         }
         
     }
    
    
      
      # plot(c(1,1),plot=FALSE)
      #l <- legend(0, 0, bty='n', l1,plot=FALSE, pch = pch_per_group, pt.cex = 0.6)
      # calculate right margin width in ndc
      w <- 0.1 #grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
      par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
     
     iqr_xlim=2*sd(result$variates$X[,1],na.rm=TRUE) #4*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
     iqr_ylim=2*sd(result$variates$X[,2],na.rm=TRUE) #4*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
     
    # save(X,result,samplelabels,ellipse.conf.level,iqr_ylim,iqr_xlim,r1,classlabelsorig,col_per_group,main_text,analysistype,l1,pch,file="debug3.Rda")
         
      if(pca.ellipse==TRUE){
          
          suppressMessages(library(car))
                      
                 if(dim(classlabelsorig)[2]==2){
                     de1<-(dataEllipse(x=result$x[,1], y=result$x[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,center.pch=FALSE,ylim=range(pretty(c(floor(min(result$variates$X[,2])-iqr_ylim),ceiling(max(result$variates$X[,2])+iqr_ylim)))),
                     xlim=range(pretty(c(floor(min(result$variates$X[,1])-iqr_xlim),ceiling(max(result$variates$X[,1])+iqr_xlim)))),verbose=FALSE))
                 
                 }
                 else
                 {
                           
                  de1<-(dataEllipse(x=result$x[,1], y=result$x[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,
                                    levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),
                                    ylab=paste("PC2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.6,cex.axis=0.8,
                                    cex.lab=0.6,center.pch=FALSE,ylim=range(pretty(c(floor(min(result$variates$X[,2])-iqr_ylim),
                                                                                     ceiling(max(result$variates$X[,2])+iqr_ylim)))),
                                    xlim=range(pretty(c(floor(min(result$variates$X[,1])-iqr_xlim),ceiling(max(result$variates$X[,1])+iqr_xlim))))))
                     
                 }
                   le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
                   
      }
      else{
       par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
       
       if(dim(classlabelsorig)[2]==2){
     de1<-(dataEllipse(x=result$x[,1], y=result$x[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,ylim=range(pretty(c(floor(min(result$variates$X[,2])-iqr_ylim),ceiling(max(result$variates$X[,2])+iqr_ylim)))),xlim=range(pretty(c(floor(min(result$variates$X[,1])-iqr_xlim),ceiling(max(result$variates$X[,1])+iqr_xlim))))))
       }else{
           
           #de1<-(dataEllipse(x=result$x[,1], y=result$x[,2],groups=as.factor(classlabelsorig[,2]),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.7,ylim=range(pretty(c(floor(min(result$variates$X[,2])-iqr_ylim),ceiling(max(result$variates$X[,2])+iqr_ylim)))),xlim=range(pretty(c(floor(min(result$variates$X[,1])-iqr_xlim),ceiling(max(result$variates$X[,1])+iqr_xlim))))))
           
          de1<-(dataEllipse(x=result$x[,1], y=result$x[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,ylim=range(pretty(c(floor(min(result$variates$X[,2])-iqr_ylim),ceiling(max(result$variates$X[,2])+iqr_ylim)))),xlim=range(pretty(c(floor(min(result$variates$X[,1])-iqr_xlim),ceiling(max(result$variates$X[,1])+iqr_xlim))))))
           
           
       }
       
       #(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
       
      le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
       
   }
     
       if(ncomp>2){
           
           if(dim(classlabelsorig)[2]==2){
           if(do_pca_anova==TRUE){
               main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing\np-value for overall differences between groups using PC1 and PC3 in a multivariate\n ANOVA model=",round(pc2_pval,3),sep="")
           }
           }
           
          
          
            w <- 0.1
           par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
           
           iqr_xlim=2*sd(result$variates$X[,1],na.rm=TRUE) #(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
           iqr_ylim=2*sd(result$variates$X[,3],na.rm=TRUE) #(summary(result$variates$X[,3])[5]-summary(result$variates$X[,3])[3])
        
         if(pca.ellipse==TRUE){
           if(dim(classlabelsorig)[2]>=2){
      de1<-(dataEllipse(x=result$x[,1], y=result$x[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,center.pch=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
           }else{
               
               de1<-(dataEllipse(x=result$x[,1], y=result$x[,3],groups=as.factor(classlabelsorig[,2]),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,center.pch=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
               
               
               
           }
           
      #(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
      le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
      
      
         }else{
      
       w <- 0.1
      par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
      
      if(dim(classlabelsorig)[2]>=2){
      de1<-(dataEllipse(x=result$x[,1], y=result$x[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
      }else{
          
          de1<-(dataEllipse(x=result$x[,1], y=result$x[,3],groups=as.factor(classlabelsorig[,2]),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
          
          
      }
      # (legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
      le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
         }
      
      w <- 0.1
      par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
      
      iqr_xlim=2*sd(result$variates$X[,2],na.rm=TRUE) #2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
      iqr_ylim=2*sd(result$variates$X[,3],na.rm=TRUE) #2.5*(summary(result$variates$X[,3])[5]-summary(result$variates$X[,3])[3])
      
      if(dim(classlabelsorig)[2]==2){
      if(do_pca_anova==TRUE){
          
          main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing\np-value for overall differences between groups using PC2 and PC3 in a multivariate\n ANOVA model=",round(pc3_pval,3),sep="")
      }
      }
      
      if(dim(classlabelsorig)[2]==2){
          if(do_pca_anova==TRUE){
              
              main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing\np-value for overall differences between groups using PC2 and PC3 in a multivariate\n ANOVA model=",round(pc3_pval,3),sep="")
          }
          }
      
       if(pca.ellipse==TRUE){
      if(dim(classlabelsorig)[2]>=2){
      de1<-(dataEllipse(x=result$x[,2], y=result$x[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC2 (",r1[2],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,center.pch=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,2])-iqr_xlim,max(result$variates$X[,2])+iqr_xlim)))
      }else{
          
            de1<-(dataEllipse(x=result$x[,2], y=result$x[,3],groups=as.factor(classlabelsorig[,2]),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC2 (",r1[2],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,center.pch=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,2])-iqr_xlim,max(result$variates$X[,2])+iqr_xlim)))
          
      }
      
       #(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
       le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
       
       }else{
       
    
      
       
        w <- 0.1
       par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
       
        if(dim(classlabelsorig)[2]>=2){
       de1<-(dataEllipse(x=result$x[,2], y=result$x[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC2 (",r1[2],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,2])-iqr_xlim,max(result$variates$X[,2])+iqr_xlim)))
        }else{
            
            de1<-(dataEllipse(x=result$x[,2], y=result$x[,3],groups=as.factor(classlabelsorig[,2]),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC2 (",r1[2],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.6,cex.axis=0.8,cex.lab=0.6,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,2])-iqr_xlim,max(result$variates$X[,2])+iqr_xlim)))
            
            
        }
       #(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
     le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
       
       
       }
       
       }
       
       
      

}
   
   suppressMessages(library(reshape2))


   Class<-samplelabels
   
   if(ncol(result$x)>5){
       
       melted <- cbind(Class, melt(result$x[,1:5]))
       
   }else{
       
       melted <- cbind(Class, melt(result$x))
   }
   
   colnames(melted)<-c("Class","Samples","Var2","PCscore")
   melted$Samples<-seq(1,nrow(result$x))
   
   ##save(pc_pval_single,file="pc_pval_single.Rda")
   # for(i in 1:5)
   ##save(melted,result,samplelabels,Class,col,col_per_group,col_vec,file="pc_score_plots.Rda")
   testname=""
   lapply(1:min(5,ncol(result$x)),function(i){
       
       pcname<-paste("PC",i,sep="")
       melted_pc1<-melted[which(melted$Var2==pcname),]
       
       melted_pc1
       
       if(do_pca_anova==TRUE){
           
          
       main_text=paste("PC score plots using ",filename," features after preprocessing\np-value for differences between groups ",testname," using ",pcname,"\n",testname," p=",round(pc_pval_single[i],3),sep="")
       }else{
           
            main_text=paste("PC score plots using ",filename," features after preprocessing",sep="")
       }
       
   
       w <- 0.1
       par(omd=c(0, 1-w, 0, 1),cex.main=0.7)
       
       # barCenters=barplot(myData$Intensity,ylab="Intensity",xlab="Class",main=mzlabel_cur,col=barplot.col.opt1,ylim=c(0,max_yval))
       #arrows(barCenters, ymin,barCenters, ymax,angle=90,code=3,lty=1,lwd=1.25,length=0.05)
       ylab1<-paste(pcname,"score",sep="")
       (plot(as.vector(melted_pc1[,4]),col=c(col),main=main_text, ylab=ylab1,xlab="Sample",type="h",lwd=2,cex.main=0.7))
       
       le1<-(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(melted_pc1$Class), col = col_per_group,pch = rep(19,length(col_per_group)), pt.cex = 0.6, title = "Class",cex=0.8))
       
       
       #print(barplot1)
       
       
       
   })
   
    ##savelist=ls(),file="getpcadebug.Rda")

 
    return(list(result=result,pca_pval_vec=pc_pval_single))
}

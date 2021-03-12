get_survivalanalysis <-
function(Xmat,Ymat,GroupBy.variable="sex",cox.glmnet=FALSE,newdevice=TRUE,plot.width=8,plot.height=8,KMplot=TRUE)
{
  
  
  Xmat1<-cbind(colnames(Xmat[,-c(1)]),t(Xmat[,-c(1)]))
  Xmat1<-as.data.frame(Xmat1)
  cnames<-colnames(Xmat1)
  cnames[1]<-"SID"
  colnames(Xmat1)<-cnames
  
  cnames<-colnames(Ymat)
  cnames[1]<-"SID"
  cnames[2]<-"time"
  cnames[3]<-"status"
  colnames(Ymat)<-cnames
  
  
  sdata1<-merge(Ymat,Xmat1,by="SID")
  
  
  if(cox.glmnet==TRUE){
    cox.cv.glmnet<-cv.glmnet(x=data.matrix(sdata1[,-c(1:3)]),y=Surv(sdata1$time, sdata1$status),family="cox",alpha=1)
    
    cox.glmnet.res<-glmnet(x=data.matrix(sdata1[,-c(1:3)]),y=Surv(sdata1$time, sdata1$status),family="cox",lambda=cox.cv.glmnet$lambda.min)
    
    select_feats<-which(abs(cox.glmnet.res$beta)>0)
    
    
    if(length(select_feats)<0){
      stop("No features selected.")
    }
    #sdata<-sdata[,c(1,2,(select_feats+2))]
    #fit<-cox.glmnet.res	
    Xmat2<-Xmat1[,c(1,select_feats+1)]
    sdata1<-merge(Ymat,Xmat2,by="SID")
  }
  
  
  
  
 
  plist = list()
  
  if(newdevice==TRUE){
    
    pdf("CPH_results.pdf",width=plot.width,height=plot.height)
  }
  if(KMplot==TRUE){
    
    
    if(length(grep(colnames(sdata1),pattern=GroupBy.variable))>0){
      
      fit <- coxph(Surv(time, status) ~ ., data=sdata1[,-c(1)])
      p<- get_survplot_adjustedcurve(fit,data=sdata1,variable=GroupBy.variable,L_ratio_test=FALSE,waldtest=FALSE,
                                   logrank=TRUE,concordance = FALSE,add.risk.table=TRUE, add.cumulative.events=FALSE)
    
    }else{
      
      fit <- coxph(Surv(time, status) ~ ., data=sdata1[,-c(1)])
      p<- get_survplot_adjustedcurve(fit,data=sdata1,method = "single",L_ratio_test=FALSE,waldtest=FALSE,
                                     logrank=FALSE,concordance = FALSE,add.risk.table=FALSE, add.cumulative.events=FALSE)
      
    }
    
   # plist[[length(plist)+1]] <- cowplot::ggdraw(p)
    
    text_val <-"Kaplan-Meier survival curve"
    
    # Create a text grob
    tgrob <- text_grob(text_val,size = 16)
    # Draw the text
    plot_0 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,1,0,0, "cm"))
    
    print(ggarrange(plot_0,p,ncol = 1,nrow = 2,heights = c(0.75,5)))
    
  }
  
  
  
  #run univariate analysis
  res<-lapply((ncol(Ymat)+1):ncol(sdata1),function(j){
                
                intensity<-as.numeric(as.character(sdata1[,c(j)]))
                
                if(length(which(is.na(intensity)==FALSE))<1){
                  intensity<-as.factor(sdata1[,c(j)])
                  
                }
                sdata_temp<-cbind(Ymat[,-c(1)],intensity)
                
                sdata_temp<-as.data.frame(sdata_temp)
                
              
                
                cnames_all<-colnames(sdata_temp)
                
                
                cnames_other<-cnames_all[-which(cnames_all%in%c("time","status","intensity"))]
                
                
                if(length(cnames_other)>0){
                sdata_temp<-sdata_temp[,c("time","status","intensity",cnames_other)]
                }else{
                  sdata_temp<-sdata_temp[,c("time","status","intensity")]
                  
                }
                
                cnames_1<-colnames(sdata_temp)
                
                cnames_1[which(cnames_1=="intensity")]<-colnames(sdata1)[j]
                colnames(sdata_temp)<-cnames_1
                
                if(length(grep(colnames(Ymat),pattern=GroupBy.variable))>0){
                  
                  KM.method="conditional"
                 
                  
                  res<-get_survivalanalysis_univariate(sdata=sdata_temp,KMplot=TRUE,forestplot=TRUE,diagnosis=TRUE, maintitle="",
                                                       KM.variable = GroupBy.variable,
                                                       KM.reference = NULL,
                                                       KM.method = KM.method,
                                                       KM.plottitle = "",
                                                       KM.xlab = "Time", KM.ylab = "Survival probability",
                                                       KM.size = 1, KM.xincrement="auto", KM.yincrement=0.25,
                                                       KM.concordance=TRUE, KM.rsquare=TRUE, KM.L_ratio_test=FALSE, KM.waldtest=FALSE, KM.logrank=TRUE,
                                                       KM.linecolor="auto", KM.add.risk.table=TRUE, KM.add.cumulative.events=TRUE,
                                                       forest.xlab="",forest.ylab="Hazard Ratio (95% CI)",forest.sigcolor="red",
                                                       diagnosis.covariate=NULL,diagnosis.type="pha",cox.glmnet=FALSE,cv.glm.alpha=1)
                  
                  
                  text_val <- paste("Forest plot for variable: ",colnames(sdata1)[j],sep="")
                  
                  # Create a text grob
                  tgrob <- text_grob(text_val,size = 16)
                  # Draw the text
                  plot_0 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))
                  
                  print(ggarrange(plot_0,res$plots[[2]],
                            ncol = 1,nrow = 2,heights = c(1,5)))
                  
                  
                }else{
                  
                  KM.method="single"
                  
                 # print(head(sdata_temp))
                  
                  res<-get_survivalanalysis_univariate(sdata=sdata_temp,KMplot=TRUE,forestplot=TRUE,diagnosis=TRUE, maintitle="",
                                                       KM.variable = GroupBy.variable,
                                                       KM.reference = NULL,
                                                       KM.method = KM.method,
                                                       KM.plottitle = "",
                                                       KM.xlab = "Time", KM.ylab = "Survival probability",
                                                       KM.size = 1, KM.xincrement="auto", KM.yincrement=0.25,
                                                       KM.concordance=TRUE, KM.rsquare=TRUE, KM.L_ratio_test=FALSE, KM.waldtest=FALSE, KM.logrank=TRUE,
                                                       KM.linecolor="auto", KM.add.risk.table=TRUE, KM.add.cumulative.events=TRUE,
                                                       forest.xlab="",forest.ylab="Hazard Ratio (95% CI)",forest.sigcolor="red",
                                                       diagnosis.covariate=NULL,diagnosis.type="pha",cox.glmnet=FALSE,cv.glm.alpha=1)
                  
                  
                  text_val <- paste("Forest plot for CPH model for variable: ",colnames(sdata1)[j],sep="")
                  
                  # Create a text grob
                  tgrob <- text_grob(text_val,size = 16)
                  # Draw the text
                  plot_0 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))
                  
                  print(ggarrange(plot_0,res$plots[[2]],
                            ncol = 1,nrow = 2,heights = c(1,5)))
                }
               
                s1=summary(res$stats[[1]])
                
                resvec<-c(s1[[7]],s1[[7]][2]-1.96*s1[[7]][3],s1[[7]][2]+1.96*s1[[7]][3])
                
                return(resvec)
    })
  #save(res,file="res1.Rda")
  res<-ldply(res,rbind)
  rownames(res)<-colnames(sdata1)[c((ncol(Ymat)+1):ncol(sdata1))]
  try(dev.off(),silent=TRUE)
  colnames(res)<-c("coef","exp(coef)","se(coef)","z","pvalue","lower.confint.95pct","upper.confint.95pct")
  return(res)
}

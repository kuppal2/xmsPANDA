get_survivalanalysis_univariate <-
function(sdata,KMplot=TRUE,forestplot=TRUE,diagnosis=TRUE, maintitle="",
                                 KM.variable = NULL,
                                 KM.reference = NULL,
                                 KM.method = "conditional",
                                 KM.plottitle = "",
                                 KM.xlab = "Time", KM.ylab = "Survival probability",
                                 KM.size = 1, KM.xincrement="auto", KM.yincrement=0.25,
                                 KM.concordance=TRUE, KM.rsquare=TRUE, KM.L_ratio_test=TRUE, KM.waldtest=TRUE, KM.logrank=TRUE,
                                 KM.linecolor="auto", KM.add.risk.table=TRUE, KM.add.cumulative.events=TRUE,
                                 forest.xlab="",forest.ylab="Hazard Ratio (95% CI)",forest.sigcolor="red",
                                 diagnosis.covariate=NULL,diagnosis.type="pha",cox.glmnet=FALSE,cv.glm.alpha=1){
  
  suppressMessages(library(glmnet))
  suppressMessages(library(cowplot))
  
  
#  print("Survival analysis is running.")
  #print("Using the first three columns as follow-up time, status, intensity.")
  outlist = list()
  
  #subject_id=as.vector(data[,-1])

  
  if(class(sdata[,1])=="numeric" & class(sdata[,2])=="numeric"){ # & class(sdata[,3])=="numeric"){
    
    #colnames(sdata)[1:3]=c("time","status","intensity")
    
    #fit <- coxph(Surv(time, status) ~ ., data=sdata)
    
   
    
    #colnames(sdata)[1:3]=c("time","status","intensity")
    fit <- coxph(Surv(time, status) ~ ., data=sdata)
    
    outlist[[length(outlist)+1]] <- fit
    names(outlist)[length(outlist)] <- "fitted_model"
  
      
    
    #print(summary(fit))
    
    plist = list()
    
    if(KMplot==TRUE){
      
      p<- get_survplot_adjustedcurve(fit,data=sdata,variable=KM.variable,reference=KM.reference,method=KM.method,
                                     plottitle=KM.plottitle,xlab=KM.xlab,ylab=KM.ylab,size=KM.size,
                                     xincrement=KM.xincrement,yincrement=KM.yincrement,concordance=KM.concordance,
                                     rsquare=KM.rsquare,L_ratio_test=KM.L_ratio_test,waldtest=KM.waldtest,
                                     logrank=KM.logrank,linecolor=KM.linecolor,add.risk.table=KM.add.risk.table,add.cumulative.events=KM.add.cumulative.events)
      
      plist[[length(plist)+1]] <- cowplot::ggdraw(p)
    }
    
    
    if(forestplot==TRUE){
      
      p <- get_forestplot(fit,data=sdata,xlab=forest.xlab,ylab=forest.ylab,sigcolor=forest.sigcolor)
      
      plist[[length(plist)+1]] <- cowplot::ggdraw(p)
      
    }
    
    
    if(length(plist)!=0){
      
      pcombine <- cowplot::plot_grid(plotlist = plist ,ncol =2)
      title <- cowplot::ggdraw() + cowplot::draw_label(maintitle, fontface='bold', size=16)
      pcombine <- cowplot::plot_grid(title,pcombine,nrow =2,rel_heights=c(0.1, 1))
      
      outlist[[length(outlist)+1]] <- pcombine
      names(outlist)[length(outlist)] <- "analysisplots"
    }
    
    if(diagnosis==TRUE){
      
      pdiagnosis <- get_diagnosis(fit,data=sdata,covariate=diagnosis.covariate,type=diagnosis.type,maintitle=maintitle)
      
      outlist[[length(outlist)+1]] <- pdiagnosis
      names(outlist)[length(outlist)] <- "diagnosis"
      
    }
    
    names(outlist)
   # print(plist)
    return(list("stats"=outlist,"plots"=plist))
    
  }else{
    
    print("Error occured: There are something wrong in the first three columns. Please double check your input data and make sure there is no missing value and all numeric values.")
    
  }
  
  
}

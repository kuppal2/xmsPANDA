get_forestplot <-
function(fit,data,xlab="",ylab="Hazard Ratio (95% CI)",sigcolor="red"){
  
  indata <- .get_data(fit, data)
  #save(fit,data,file="fit.Rda")
  fit_summary <- summary(fit)
  coef_table <- fit_summary$coefficients
  plotdata <- data.frame(predictor=as.character(rownames(coef_table)),
                         hazard_ratio=as.numeric(round(coef_table[,2],3)),
                         upper=as.numeric(coef_table[,2])+1.96*as.numeric(coef_table[,3]),
                         lower=as.numeric(coef_table[,2])-1.96*as.numeric(coef_table[,3]),
                         pvalue=as.numeric(round(coef_table[,5],3)),
                         stringsAsFactors=FALSE
  )
  plotdata$predictor <- factor(plotdata$predictor,levels=plotdata$predictor)
  plotdata[,"sig"]<-ifelse(plotdata$pvalue<0.05,"Yes","No")
  plotdata[plotdata$pvalue<0.05,"annotation"] <- paste(round(plotdata[plotdata$pvalue<0.05,"hazard_ratio"],3)," (p value ",round(plotdata[plotdata$pvalue<0.05,"pvalue"],3),"*)",sep="")
  plotdata[plotdata$pvalue<0.01,"annotation"] <- paste(round(plotdata[plotdata$pvalue<0.01,"hazard_ratio"],3)," (p value ",round(plotdata[plotdata$pvalue<0.01,"pvalue"],3),"**)",sep="")
  plotdata[plotdata$pvalue<0.001,"annotation"] <- paste(round(plotdata[plotdata$pvalue<0.001,"hazard_ratio"],3)," (p value ",round(plotdata[plotdata$pvalue<0.001,"pvalue"],3),"***)",sep="")
  plotdata[is.na(plotdata$annotation),"annotation"] <- paste(round(plotdata[is.na(plotdata$annotation),"hazard_ratio"],3)," (p value ",round(plotdata[is.na(plotdata$annotation),"pvalue"],3),")",sep="")
  
  p <- ggplot(data=plotdata, aes(x=predictor, y=hazard_ratio, ymin=lower, ymax=upper, color=sig)) +
    geom_errorbar(width=0.2) +
    geom_point() +
    geom_text(aes(label=plotdata$annotation),vjust=-1, hjust=0.5) +
    geom_hline(yintercept=1, lty=2) +
    coord_flip() +
    xlab(xlab) + ylab(ylab) +
    scale_x_discrete(limits = rev(plotdata$predictor)) +
    scale_colour_manual(name ='Significant', values =c('No'='black','Yes'=sigcolor), labels = c('No','Yes')) +
    theme(axis.title=element_text(size=10),
          axis.text =element_text(size=15),
          axis.line = element_line(size=1, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.key=element_blank(),
          legend.title=element_text(size=13),
          legend.position="none",
          plot.title = element_text(size=18, hjust = 0),
          axis.text.x=element_text(colour="black", size = 12,angle=0),
          axis.text.y=element_text(colour="black", size = 12))
  return(p)
}

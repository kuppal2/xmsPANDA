get_survplot_adjustedcurve <-
function(fit,
                                       variable = NULL,
                                       data = NULL,
                                       reference = NULL,
                                       method = "conditional",
                                       plottitle = "",
                                       xlab = "Time", ylab = "Survival probability",
                                       size = 1, xincrement="auto", yincrement=0.25,
                                       concordance=TRUE, rsquare=FALSE, L_ratio_test=FALSE, waldtest=TRUE, logrank=FALSE,
                                       linecolor="auto", add.risk.table=TRUE, add.cumulative.events=FALSE
) {
  stopifnot(method %in% c("marginal", "average", "conditional", "single"))
  indata <- .get_data(fit, data)
  
  # deal with default arguments
  # reference = NULL
  if (is.null(reference))
    reference <- indata
  
  # variable = NULL
  if (is.null(variable)) {
    # is there a 'strata' component?
    term.labels <- attr(terms(fit$formula), "term.labels")
    strata.term.labels <- grep(term.labels, pattern = "strata(", fixed = TRUE, value = TRUE)
    if (length(strata.term.labels) > 0) {
      variable <- gsub(
        gsub(
          strata.term.labels,
          pattern = "strata(", replacement = "", fixed = TRUE)[1],
        pattern = "[\\) ]", replacement = "")
      cat("The variable argument is missing. Using", variable, "as extracted from strata\n")
    } else {
      # if not then leave variable = NULL
      method = "single"
    }
  }
  
  if(method=="single"){
    
    add.risk.table=FALSE; add.cumulative.events=FALSE
  }
  
  
  curve <- switch(method,
                  single = ggadjustedcurves.single(indata, fit),
                  average =  ggadjustedcurves.average(indata, fit, variable),
                  conditional = suppressWarnings(ggadjustedcurves.conditional(indata, fit, variable, reference)),
                  marginal = ggadjustedcurves.marginal(indata, fit, variable))
  
  timepoints <- unique(curve$time)
  
  mintime <- 0
  maxtime <- ceiling(max(timepoints)/100)*100
  if(xincrement=="auto"){
    xincrement <- maxtime/5
  }else{
    xincrement <- xincrement
  }
  
  sfit <- summary(fit)
  ydecrease <- 0
  if(concordance==TRUE){
    label_concordance <- str_pad(paste("Concordance = ",format(round(sfit$concordance[1],3),nsmall=3), sep=""),46,side="right",pad=" ")
    annotation_concordance <- annotate("text", label = label_concordance, x = (maxtime - 120), y = (1-ydecrease))
    ydecrease <- ydecrease +0.04
  }else{
    annotation_concordance = NULL
  }
  if(rsquare==TRUE){
    label_rsquare <- str_pad(paste("R-square = ",format(round(sfit$rsq[1],3),nsmall=3), sep=""),49,side="right",pad=" ")
    annotation_rsquare <- annotate("text", label = label_rsquare, x = (maxtime - 120), y = (1-ydecrease))
    ydecrease <- ydecrease +0.04
  }else{
    annotation_rsquare = NULL
  }
  if(L_ratio_test==TRUE){
    label_LRT <- str_pad(paste("Likelihood ratio test p-value = ",format(round(sfit$logtest[3],3),nsmall=3),sep=""),39,side="right",pad=" ")
    annotation_LRT <- annotate("text", label = label_LRT, x = (maxtime - 120), y = (1-ydecrease))
    ydecrease <- ydecrease +0.04
  }else{
    annotation_LRT = NULL
  }
  if(waldtest==TRUE){
    label_waldtest <- str_pad(paste("Wald test p-value = ",format(round(sfit$waldtest[3],3),nsmall=3),sep=""),43,side="right",pad=" ")
    annotation_waldtest <- annotate("text", label = label_waldtest, x = (maxtime - 120), y = (1-ydecrease))
    ydecrease <- ydecrease +0.04
  }else{
    annotation_waldtest = NULL
  }
  if(logrank==TRUE){
    label_logrank <- str_pad(paste("Log-rank test p-value = ",format(round(sfit$sctest[3],3),nsmall=3),sep=""),41,side="right",pad=" ")
    annotation_logrank <- annotate("text", label = label_logrank, x = (maxtime - 120), y = (1-ydecrease))
    ydecrease <- ydecrease +0.04
  }else{
    annotation_logrank = NULL
  }
  
  if( length(linecolor)==1 & linecolor[1]=="auto"){
    scale_color_value <- NULL
    table_group_color <- NULL
  }else {
    group <- levels(factor(curve$variable))
    names(linecolor) <- group
    scale_color_value <- scale_colour_manual(values =linecolor)
    table_group_color <- linecolor
  }
  
  plot <- list(); i<-0
  
  p1 <- ggplot(curve, aes(x = time, y = surv)) +
    geom_step(aes(color=factor(variable)), size=size) +
    labs(color=variable) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(plottitle) +
    scale_color_value +
    annotation_concordance + annotation_rsquare + annotation_LRT + annotation_waldtest + annotation_logrank +
    scale_x_continuous(breaks = seq(mintime, maxtime, by = xincrement), lim=c(0,maxtime)) +
    scale_y_continuous(breaks = seq(0, 1, by = yincrement)) +
    theme(axis.title=element_text(size=12),
          axis.text =element_text(size=12),
          axis.line = element_line(size=1, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.key=element_blank(),
          legend.title=element_text(size=13),
          legend.position="bottom",
          plot.title = element_text(size=18, hjust = 0),
          axis.text.x=element_text(colour="black", size = 11,angle=0),
          axis.text.y=element_text(colour="black", size = 11))
  
  if(method=="single"){
    
    p1 <- p1 + theme(legend.position = "none")
  }
  i <- i+1; plot[[i]]  <- p1
  
  if(add.risk.table==TRUE){
    
    p2 <- ggsurvtable(survfit(as.formula(paste(as.character(fit$call$formula)[2],"~",variable,sep=" ")),data=indata),
                      data=indata, survtable=c("risk.table"),
                      break.time.by = xincrement, xlim=c(0,maxtime),ylab="",xlsab="",
                      legend.labs=levels(factor(indata$variable)), color = "strata", palette = table_group_color)
    p2 <- p2 + theme(legend.position="none",
                     plot.margin = unit(c(5.5,5.5,5.5,28), "pt"),
                     axis.line = element_blank())
    i <- i+1; plot[[i]]  <- p2
    
  }
  
  if(add.cumulative.events==TRUE){
    
    p3 <- ggsurvtable(survfit(as.formula(paste(as.character(fit$call$formula)[2],"~",variable,sep=" ")),data=indata),
                      data=indata, survtable=c("cumevents"),
                      break.time.by = xincrement, xlim=c(0,maxtime),ylab="",xlsab="",
                      legend.labs=levels(factor(indata$variable)), color = "strata", palette = table_group_color)
    p3 <- p3 + theme(legend.position="none",
                     plot.margin = unit(c(5.5,5.5,5.5,28), "pt"),
                     axis.line = element_blank())
    i <- i+1; plot[[i]]  <- p3
    
  }
  
  if(length(plot)==1){
    
    #p <- grid.arrange(plot[[1]])
    p <- cowplot::plot_grid(plotlist = plot, nrow=1)
  } else if (length(plot)==2) {
    
    height1 <- 0.85 - (length(levels(curve$variable))-2)*0.3 ; height2 <- 0.15 + (length(levels(curve$variable))-2)*0.3
    p <- cowplot::plot_grid(plotlist = plot, nrow=2, rel_heights=c(height1, height2))
    #p <- grid.arrange(plot[[1]], plot[[2]], nrow = 2, heights = c(height1, height2))
  } else if (length(plot)==3) {
    
    height1 <- 0.7 - (length(levels(curve$variable))-2)*0.6
    height2 <- 0.15 + (length(levels(curve$variable))-2)*0.3
    height3 <- 0.15 + (length(levels(curve$variable))-2)*0.3
    p <- cowplot::plot_grid(plotlist = plot, nrow=3, rel_heights=c(height1, height2, height3))
    #p <- grid.arrange(plot[[1]], plot[[2]], plot[[3]], nrow = 3, heights = c(height1, height2, height3))
    
  
  }
  #print(plot)
  return(p)
}

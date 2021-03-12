get_diagnosis <-
function(fit,data,covariate=NULL,type="pha",maintitle=""){
  
  
  
  indata <- .get_data(fit, data)
  custom_theme <- theme(axis.title=element_text(size=15),
                        axis.text =element_text(size=12),
                        axis.line = element_line(size=1, colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        legend.key=element_blank(),
                        legend.title=element_text(size=13),
                        legend.position="none",
                        plot.title = element_text(size=18, hjust = 0),
                        strip.background =element_blank(),
                        strip.text = element_text(size=12),
                        axis.text.x=element_text(colour="black", size = 12,angle=0),
                        axis.text.y=element_text(colour="black", size = 12))
  
  if(type=="pha"){
    ## check the proportional hazards assumption for the Cox model
    ptmp <- ggcoxzph(cox.zph(fit), ggtheme = custom_theme)
    p <- cowplot::plot_grid(plotlist = ptmp, ncol=2)
    p <- cowplot::add_sub(p,"Notice: The proportional hazard assumption is supported by the non-significant p-value", y = 0.5, vjust = 0, colour = "red")
    p <- cowplot::ggdraw(p)
    title <- cowplot::ggdraw() + cowplot::draw_label(maintitle, fontface='bold', size=16)
    p <- cowplot::plot_grid(title,p,nrow =2,rel_heights=c(0.1, 1))
    
  } else if (type=="io"){
    ## check influential observations
    p <- ggcoxdiagnosis(fit, type = "dfbeta", linear.predictions = FALSE, ggtheme = custom_theme, ylab="Residuals",title = maintitle)
  } else if (type=="nl" & !is.null(covariate)){
    ## check non linearity in relationship between the log hazard and the covariates
    newfunction <- coxph(as.formula(paste(as.character(fit$call$formula)[2],"~",covariate,"+log(",covariate,")+sqrt(",covariate,")",sep="")),data=indata)
    p <- ggcoxfunctional(newfunction,data=indata)
    p <- cowplot::plot_grid(plotlist = p, ncol=2)
    title <- cowplot::ggdraw() + cowplot::draw_label(maintitle, fontface='bold', size=16)
    p <- cowplot::plot_grid(title,p,nrow =2,rel_heights=c(0.1, 1))
  } else {
    stop("This function only supports three types: 'pha', 'io' and 'nl'. If chossing 'nl', please make sure the covariate is continuous.")
  }
  
}

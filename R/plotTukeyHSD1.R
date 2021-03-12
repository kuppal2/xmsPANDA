plotTukeyHSD1 <-
function(tukey.res,
                                          x.axis.label = "Comparison",
                                          y.axis.label = "Effect Size",
                                          plot.cex.axis=0.7,var.name=NA){
  

  res1<-lapply(1:length(tukey.res),function(j)
  {
  tukey.out <- as.data.frame(tukey.res[[j]])
  means <- tukey.out$diff
  categories <- row.names(tukey.out)
 # categories<-factor(categories,levels=unique(categories))
  groups <- length(categories)
  ci.low <- tukey.out$lwr
  ci.up  <- tukey.out$upr  
  padj<-tukey.out$`p adj`
  padjbool<-rep(0,length(padj))
  padjbool<-replace(padjbool,which(padj<0.05),1)
  
  tmp1<-cbind(categories,means,ci.low,ci.up,padj,padjbool)
  tmp1<-as.data.frame(tmp1)
  tmp1[,-c(1)]<-apply(tmp1[,-c(1)],2,function(x){as.numeric(as.character(x))})
  cex.plots=plot.cex.axis
  tmp1$categories<-factor(tmp1$categories,levels=unique(tmp1$categories))
 
  p1=ggplot(tmp1,aes(y=means,x=categories,colour=factor(padjbool),show.legend = FALSE)) + geom_point() + 
    geom_errorbar(aes(ymin=ci.low,ymax=ci.up),width=0.1) + coord_flip()+ #scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    geom_hline(yintercept = 0,lty=2,col="blue") + theme_bw()+ ylab(as.character(y.axis.label[j])) + xlab("Comparison")
 
  if(is.na(var.name)==FALSE){
   
      title.sub<-paste(" (",var.name,")",sep="") 
  }else{
    
      title.sub=""
  }
       if(j==1){
        
        p1=p1+ggtitle(paste("TukeyHSD: Factor 1",title.sub,sep=""))
      }else{
        if(j==2){
          
          p1=p1+ggtitle(paste("TukeyHSD: Factor 2",title.sub,sep=""))
        }else{
          
          if(j==3){
            
            p1=p1+ggtitle(paste("TukeyHSD: Factor 1x2",title.sub,sep=""))
          }
        }
        
      }
  
  p1=p1 + labs(fill="p<0.05") + font("legend.text", size = 10*cex.plots, color = "black") + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                        panel.grid.minor = element_blank(),
                                                                                                        panel.spacing=unit(1,"lines"),
                                                                                                        axis.line = element_line(colour = "black",size=1),
                                                                                                        axis.text= element_text(size=9*cex.plots), axis.title=element_text(size=11*cex.plots,face="bold"),
                                                                                                        plot.title = element_text(hjust = 0.5,size=12*cex.plots),
                                                                                                        axis.ticks.length = unit(-0.05, "in"),
                                                                                                        axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
                                                                                                        axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
                                                                                                        axis.ticks.x = element_blank(),
                                                                                                        aspect.ratio = 1,
                                                                                                        legend.background = element_rect(color = "black", fill = "white"),
                                                                                                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                        strip.text = element_text(face="bold")) + scale_color_manual(values=c("blue","red"),labels=c("padj>=0.05","padj<0.05")) 

 #p1=p1+aes(show.lege)
   p1=p1+theme(legend.position="none")
  #}   
  
  return(p1) # + xlab(x.axis.label)
  })
  
  res2<-lapply(1:1,function(j)
  {
    tukey.out <- as.data.frame(tukey.res[[j]])
    means <- tukey.out$diff
    categories <- row.names(tukey.out)
    # categories<-factor(categories,levels=unique(categories))
    groups <- length(categories)
    ci.low <- tukey.out$lwr
    ci.up  <- tukey.out$upr  
    padj<-tukey.out$`p adj`
    padjbool<-rep(0,length(padj))
    padjbool<-replace(padjbool,which(padj<0.05),1)
    
    tmp1<-cbind(categories,means,ci.low,ci.up,padj,padjbool)
    tmp1<-as.data.frame(tmp1)
    tmp1[,-c(1)]<-apply(tmp1[,-c(1)],2,function(x){as.numeric(as.character(x))})
    cex.plots=plot.cex.axis
    tmp1$categories<-factor(tmp1$categories,levels=unique(tmp1$categories))
    
    p1=ggplot(tmp1,aes(y=means,x=categories,colour=factor(padjbool),show.legend = FALSE)) + geom_point() + 
      geom_errorbar(aes(ymin=ci.low,ymax=ci.up),width=0.1) + coord_flip()+ #scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
      geom_hline(yintercept = 0,lty=2,col="blue") + theme_bw()+ ylab(as.character(y.axis.label[j])) 
    
    
    
    p1=p1 + labs(fill="p<0.05") + font("legend.text", size = 8*cex.plots, color = "black") + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                   panel.grid.minor = element_blank(),
                                                                                                   panel.spacing=unit(1,"lines"),
                                                                                                   axis.line = element_line(colour = "black",size=1),
                                                                                                   axis.text= element_text(size=9*cex.plots), axis.title=element_text(size=11*cex.plots,face="bold"),
                                                                                                   plot.title = element_text(hjust = 0.5,size=12*cex.plots),
                                                                                                   axis.ticks.length = unit(-0.05, "in"),
                                                                                                   axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
                                                                                                   axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
                                                                                                   axis.ticks.x = element_blank(),
                                                                                                   aspect.ratio = 1,
                                                                                                   legend.background = element_rect(color = "black", fill = "white"),
                                                                                                   strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                                                                                   strip.text = element_text(face="bold")) + scale_color_manual(values=c("blue","red"),labels=c("padj>=0.05","padj<0.05")) 
    
    #p1=p1+aes(show.lege)
    p1=p1+theme(legend.title=element_blank())
    #}   
    
    return(p1) # + xlab(x.axis.label)
  })
  
  #save(res1,file="res1.Rda")
  g1=ggarrange(res1[[1]],res1[[2]],res1[[3]],ggpubr::get_legend(res2[[1]]),nrow=2,ncol=2)
  return(g1)
}

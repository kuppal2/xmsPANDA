get_bubbleplot <-
function(Xmat, GroupBy=FALSE,xlab.name="X",ylab.name="Y",cex.plots=0.8,newdevice=TRUE,
                         filename="bubbleplot",file.format="pdf",
                         plot.width=8,plot.height=8,plot.res=300,color.scale=c("blue","red"),statistic.type="pvalue",
                         reverse.yaxis=FALSE){
  options(warn=-1)
  
  if(GroupBy==TRUE){
    
    colnames(Xmat)<-c("X","Y","Statistic","GroupBy")
  }else{
    colnames(Xmat)<-c("X","Y","Statistic")
    
  }
  if(statistic.type=="pvalue"){
  p=ggplot(Xmat,
           aes(x = X, y = Y)) +  suppressWarnings(geom_point(aes(size = abs(Statistic)), pch = 21, show.legend= TRUE))
  }else{
    
    if(statistic.type=="correlation"){
    p=ggplot(Xmat,
             aes(x = X, y = Y)) +  suppressWarnings(geom_point(aes(size = abs(Statistic)), pch = 21, show.legend= TRUE))
    }else{
      
      p=ggplot(Xmat,
               aes(x = X, y = Y)) +  suppressWarnings(geom_point(aes(size = abs(Statistic)), pch = 21, show.legend= TRUE))
    }
  }
  if(GroupBy==TRUE){
    p=p+  facet_wrap(~ GroupBy, scale="free_x")
    
  }
  p=p+labs(x=xlab.name,y=ylab.name)
  p=p +aes(fill = Statistic) + theme_bw()
  p=p+scale_fill_gradient2(low=color.scale[1],mid="white",high=color.scale[2])+scale_size(range=c(0.5,12))
  
  #scale_size(range=c(floor(min(abs(Xmat$Statistic),na.rm=TRUE)),ceiling(max(abs(Xmat$Statistic),na.rm=TRUE))))
  
  #scale_size(range = c(0.5, 12)
             
  p=p+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing=unit(1,"lines"),
            axis.line = element_line(colour = "black",size=1),
            axis.text= element_text(size=14*cex.plots), axis.title=element_text(size=18*cex.plots,face="bold"),
            plot.title = element_text(hjust = 0.5,size=18*cex.plots),
            axis.ticks.length = unit(-0.05, "in"),
            axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
            axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
            axis.ticks.x = element_blank(),
            # aspect.ratio = 1,
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #p=p+scale_fill_manual(values = c("blue","red"))
  
  if(reverse.yaxis==TRUE){
      p=p+scale_y_discrete(limits=rev)
  }
  p=p+theme(axis.title = element_blank())
  if(newdevice==TRUE){
    
    if(file.format=="pdf"){
      fname=paste(filename,".pdf",sep="") 
       pdf(fname,width=plot.width,height=plot.height)
    }else{
      fname=paste(filename,".png",sep="") 
      png(fname,width=plot.width,height=plot.height,res=plots.res,type="cairo",units="in")
    }
    print(p)
    dev.off()
  }else{
    #print(p)
  }
  
  options(warn=0)
  return(p)
}

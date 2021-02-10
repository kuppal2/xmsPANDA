get_manhattanplots <-
function(xvec,yvec,up_or_down,maintext="",ythresh=0.05,y2thresh=NA,ylab,xlab,colorvec=c("darkblue","red3"),
                             col_seq=c("brown","chocolate3","orange3","coral","pink","skyblue","blue","darkblue","purple","violet"),xincrement=100,yincrement=1,
                             pchvec=c(21,21),background.points.col="black",bad.feature.index=NA){
  
  d4<-xvec
  min_val<-min(c(0,d4),na.rm=TRUE)[1]
  max_val<-max(d4,na.rm=TRUE)[1]
  
  windowsize=xincrement
  
  d4<-as.vector(d4)
  pvalues<-as.vector(yvec)
  
  logp<-as.vector(yvec)
  
  if(is.na(up_or_down)==TRUE){
    up_or_down<-rep(1,length(yvec))
  }
  
  max_yval<-max(yvec,na.rm=TRUE)[1]+1.96*(sd(yvec,na.rm=TRUE)/(sqrt(length(yvec))))
  
  
  plot(d4,logp,xaxt="n",ylab=ylab,xlab=xlab,xaxt="n",yaxt="n",cex=0.4,cex.main=0.7,main=maintext,ylim=range(pretty(c(0,max(logp)))))
  
  axis(1, at=seq(min_val , max_val, by=xincrement) , las=2)
  axis(2, at=seq(0 , (max(logp)+2), by=yincrement) , las=2)
  
  if(length(col_seq)>1){
    
    s1<-seq(windowsize,max_val,windowsize)
    points(d4[which(d4>=0 & d4<=windowsize)],logp[which(d4>=0 & d4<=windowsize)],col=col_seq[1],cex=0.4,pch=21,bg=background.points.col)
    for(i in 1:(length(s1)-1))
    {
      points(d4[which(d4>s1[i] & d4<=s1[i+1])],logp[which(d4>s1[i] & d4<=s1[i+1])],col=col_seq[i+1],cex=0.4,pch=21,bg=background.points.col)
    }
  }else{
    
    #points(d4[which(d4>=0 & d4<=windowsize)],logp[which(d4>=0 & d4<=windowsize)],col="black",bg=background.points.col,cex=0.4,pch=21)
    
    points(d4,logp,col=background.points.col,bg=background.points.col,cex=0.4,pch=21)
  }
  
  if(is.na(y2thresh)==TRUE){
    
    goodip<-which(yvec>ythresh)
  }else{
    
    goodip<-which(yvec>y2thresh)
    
  }
  
  if(length(bad.feature.index)>0){
    if(is.na(bad.feature.index)==FALSE){
      if(length(goodip)>0){
        
        
        
        check_bad_feat_index<-which(goodip%in%bad.feature.index)
        
        if(length(check_bad_feat_index)>0){
          goodip<-goodip[-check_bad_feat_index]
        }
        
      }
    }
  }
  
  for(i in goodip){
    if(up_or_down[i]>0){
      points(d4[i],logp[i],col=colorvec[1],cex=0.8,pch=pchvec[1],bg=colorvec[1]); points(d4[i],logp[i],col=colorvec[1],cex=0.2,bg=colorvec[1])
    }else{
      
      points(d4[i],logp[i],col=colorvec[2],cex=0.8,pch=pchvec[2],bg=colorvec[2]); points(d4[i],logp[i],col=colorvec[2],cex=0.2,bg=colorvec[2])
    }
  }
  
  if(length(bad.feature.index)>0){
    if(is.na(bad.feature.index)==FALSE){
      for(i in bad.feature.index){
        
        points(d4[i],logp[i],col=background.points.col,cex=0.4,pch=pchvec[1],bg=background.points.col); #points(d4[i],logp[i],col=colorvec[1],cex=0.2,bg=)
      }
    }
  }
  
  if(length(goodip)>0){
    #hfdrfdrthresh<-logp[which(logp==min(logp[which(yvec>ythresh)],na.rm=TRUE))]
    #abline(h=hfdrfdrthresh,col="gray8",lty=2,lwd=2)
    
    abline(h=ythresh,col="gray8",lty=2,lwd=0.8)
    if(is.na(y2thresh)==FALSE){
      abline(h=y2thresh,col="gray8",lty=2,lwd=0.8)
    }
  }
  
}

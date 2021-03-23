get_volcanoplots <-
function(xvec,yvec,up_or_down,maintext="",ythresh=0.05,y2thresh=NA,
                           ylab,xlab,colorvec=c("darkblue","red3"),
                           col_seq=c("brown","chocolate3","orange3","coral",
                                     "pink","skyblue","blue","darkblue","purple","violet"),
                           xincrement=1,yincrement=1,xthresh=1,pchvec=c(21,21)
                           ,background.points.col="gray50",bad.feature.index=NA,xlim.arg=NA,ylim.arg=NA){
  
  # ###savelist=ls(),file="volcano.Rda")
  d4<-xvec
  min_val<-round(min(d4,na.rm=TRUE)+0.5)
  max_val<-round(max(d4,na.rm=TRUE)+0.5)
  
  windowsize=xincrement
  
  d4<-as.vector(d4)
  
  logp<-as.vector(yvec)
  
  if(is.na(up_or_down)==TRUE){
    up_or_down<-rep(1,length(yvec))
  }
  
  if(is.na(xlim.arg)==TRUE & is.na(ylim.arg)==TRUE){
  plot(d4,logp,xaxt="n",ylab=ylab,xlab=xlab,xaxt="n",yaxt="n",cex=0.4,cex.main=0.7,main=maintext)
    axis(1, at=seq(min_val , max_val, by=xincrement) , las=2)
    axis(2, at=seq(0 , (max(logp)+2), by=yincrement) , las=2)
    
    
   }else{
    plot(d4,logp,xaxt="n",ylab=ylab,xlab=xlab,xaxt="n",yaxt="n",cex=0.4,cex.main=0.7,main=maintext,xlim=xlim.arg,ylim=ylim.arg)
     axis(1, at=seq(xlim.arg[1] , xlim.arg[2], by=xincrement) , las=2)
     axis(2, at=seq(ylim.arg[1],ylim.arg[2], by=yincrement) , las=2)
     
  }
 
  
  points(d4,logp,col=background.points.col,cex=0.4,bg=background.points.col,pch=21)
  points(d4,logp,col=background.points.col,cex=0.4,bg=background.points.col,pch=21)
  
  
 # goodip<-which(yvec>y2thresh & abs(xvec)>xthresh)
  
  if(is.na(y2thresh)==TRUE){
    
    goodip<-which(yvec>ythresh & abs(xvec)>xthresh)
  }else{
    
    goodip<-which(yvec>y2thresh & abs(xvec)>xthresh)
    
  }
  
  
  if(length(bad.feature.index)>0){
    
    if(is.na(bad.feature.index)==FALSE){
      if(length(goodip)>0){
        
        
        check_bad_feat_index<-which(goodip%in%bad.feature.index)
        
        if(length(check_bad_feat_index)>0){
          goodip<-goodip[-check_bad_feat_index]
        }
        #goodip<-goodip[-which(goodip%in%bad.feature.index)] #goodip[-bad.feature.index]
      }
    }
  }
  
  for(i in goodip){
    if(up_or_down[i]>0){
      points(d4[i],logp[i],col=colorvec[1],cex=0.8,pch=pchvec[1],bg=colorvec[1]); points(d4[i],logp[i],col=colorvec[1],cex=0.4,bg=colorvec[1])
    }else{
      
      points(d4[i],logp[i],col=colorvec[2],cex=0.8,pch=pchvec[2],bg=colorvec[2]); points(d4[i],logp[i],col=colorvec[2],cex=0.4,bg=colorvec[2])
    }
  }
  if(length(bad.feature.index)>0){
    
    if(is.na(bad.feature.index)==FALSE){
      for(i in bad.feature.index){
        points(d4[i],logp[i],col=background.points.col,cex=0.4,bg=background.points.col,pch=21)
        
      }
    }
  }
  
  if(length(goodip)>0){
    
    
    abline(v=(-1)*xthresh,col="gray8",lty=2,lwd=0.8)
    abline(v=xthresh,col="gray8",lty=2,lwd=0.8)
    
    abline(h=ythresh,col="gray8",lty=2,lwd=0.8)
    
    if(is.na(y2thresh)==FALSE){
      abline(h=y2thresh,col="gray8",lty=2,lwd=0.8)
      
    }
    
  }
  
  
}

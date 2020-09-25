get_fcs_child_v1.0.8.53 <-
function(metab_data,reference_sets,itrs=1000,fcs.min.hits=2,numnodes=2){


 
 # colnames(reference_sets)<-c("SetID","Name","XID")
   count.unique.formula.overlapsize=TRUE
   dup.feature.check=TRUE
   
   
   cnames1<-colnames(metab_data)
   
   cnames1[1:2]<-c("XID","Statistic")
   
   colnames(metab_data)<-cnames1
   
   scale_stat<-scale(metab_data$Statistic)
   
   metab_data<-cbind(metab_data,scale_stat)
   colnames(metab_data)<-c(cnames1,"Scaled.Statistic")
   
 
   set_size=table(reference_sets$SetID)
   set_size<-cbind(names(set_size),set_size)
   set_size<-as.data.frame(set_size)
   colnames(set_size)<-c("SetID","Total.Size")
   
   
    metab_data_sets<-merge(metab_data,reference_sets,by="XID")
    
    # save(metab_data_sets,metab_data,reference_sets,set_size,file="metab_data_sets0.Rda")
   
   metab_data_sets<-metab_data_sets[order(metab_data_sets$SetID,metab_data_sets$Statistic,decreasing=TRUE),]
   
   if(count.unique.formula.overlapsize==TRUE){
   
   if(length(grep(colnames(metab_data_sets),pattern="Formula"))>0){
            xid_formula<-paste(metab_data_sets$SetID,"_",metab_data_sets$Formula,sep="")
            
            dup_check<-which(duplicated(xid_formula)==TRUE)
            if(length(dup_check)>0){
                metab_data_sets<-metab_data_sets[-dup_check,]
            }
   }
   }
   #  write.table(metab_data_sets,file="fcs_sets_with_statistic.txt",sep="\t",row.names=FALSE)
   metab_data_sets<-metab_data_sets[,c(1:5)]
    
  
    
    if(nrow(metab_data_sets)<1){
        
        metabset_res_mat<-{}
        
    }else{
            sid_list<-metab_data_sets$SetID;
            sid_list<-unique(sid_list)
           
            metab_data_sets$Statistic<-as.numeric(as.character(metab_data_sets$Statistic))
            
            metab_data_sets$Scaled.Statistic<-as.numeric(as.character(metab_data_sets$Scaled.Statistic))
            
              metab_data_sets<-merge(metab_data_sets,set_size,by="SetID")
              
               cl<-makeCluster(max(numnodes,detectCores()*0.5))
              
              agg_res<-parLapply(cl,1:length(sid_list),function(c,metab_data_sets,sid_list){
                  
                    
                        x<-metab_data_sets$Statistic[which(metab_data_sets$SetID==sid_list[c])]
                        
                        x<-as.numeric(as.character(x))
                        
                        scaledx<-metab_data_sets$Scaled.Statistic[which(metab_data_sets$SetID==sid_list[c])]
                        
                        scaledx<-abs(as.numeric(as.character(scaledx)))
                        
                        set_size=as.numeric(as.character(metab_data_sets$Total.Size[which(metab_data_sets$SetID==sid_list[c])]))
                        
                        agg.stat<-sum((x))/set_size
                        
                        agg.scaled.stat<-sum(scaledx)/set_size
                        
                        z.score=sqrt(set_size)*agg.scaled.stat
                        
                        max_mean<-max(c(mean(x[which(x>0)],na.rm=TRUE),mean(x[which(x<0)],na.rm=TRUE)),na.rm=TRUE)/set_size
                        
                        resv=cbind(as.character(sid_list[c]),round(agg.stat,3),round(agg.scaled.stat,3),round(z.score,3),round(max_mean,3),length(x),set_size,round(sum(abs(x)),3))
                        resv=as.data.frame(resv)
                        resv=unique(resv)
                        # print(resv)
                        return(resv)
                  
              },metab_data_sets=metab_data_sets,sid_list=sid_list)
              
                                                 
                                                 set_statistic1<-ldply(agg_res,rbind) #do.call(data.frame,agg_res)
                                                    colnames(set_statistic1)=c("SetID","Agg.Statistic","Agg.Scaled.Statistic","Z.score","MaxMean","Num.Hits","Total.Size","Sum.Statistic")
                                                    
                                                    
                                                    set_statistic1$MaxMean<-as.numeric(as.character(set_statistic1$MaxMean))
                                                    set_statistic1$Agg.Statistic<-as.numeric(as.character(set_statistic1$Agg.Statistic))
                                                    set_statistic1$Z.score<-as.numeric(as.character(set_statistic1$Z.score))
                                                    
                                                    set_statistic1$Sum.Statistic<-as.numeric(as.character(set_statistic1$Sum.Statistic))
                                                    
                                                    set_statistic1$Z.score<-scale(set_statistic1$MaxMean)
                                                   
                                                   # save(metab_data_sets,sid_list,set_size,set_statistic1,file="sid1.rda")
                                                    
                                                    
                                                    #metab_data_sets<-metab_data_sets[,-c(4)]
    #metabset_res<-parLapply(cl,1:length(sid_list),function(c,dup.feature.check,fcs.method,fcs.permutation.type){
        z.test1 = function(x=NA,mu=NA,std=NA,z.score=NA){

                                             one.tail.p <- NULL

                                             if(is.na(z.score)==TRUE){
                                                  z.score <- round((mean(x)-mu)/(std),3)
                                             }
                                             one.tail.p <- round(pnorm(abs(z.score),lower.tail = FALSE),3)

                                             
                                             return(one.tail.p)
                                     }
                                        if(dup.feature.check==TRUE){

                                            if(length(grep(colnames(metab_data_sets),"mz")[1])>0){ if(length(duplicated(paste(metab_data_sets$SetID,"_",metab_data_sets$mz,"_",metab_data_sets$time,sep=""))==TRUE)>0){
                                              metab_data_sets<-metab_data_sets[-which(duplicated(paste(metab_data_sets$SetID,"_",metab_data_sets$mz,"_",metab_data_sets$time,sep=""))==TRUE),]
                                            }
                                            }
                                           
                                       }
                                        
                   
                   metab_data_sets$SetID<-as.character(metab_data_sets$SetID)
                   
                   #  print("Start here")
                   #save(metab_data_sets,sid_list,itrs,file="sid2.rda")
                                 cl<-makeCluster(max(numnodes,detectCores()*0.5))
                                 
                                  clusterEvalQ(cl,library(plyr))
                                 #   a=Sys.time()
                                 #parLapply(cl,
                                 randmetabset_res<-parLapply(cl,1:itrs,function(i,metab_data_sets,sid_list){
          
                                                set.seed(i)
                                                                            all <- sample(1:nrow(metab_data_sets),nrow(metab_data_sets))
                                                                            randmetab_data_sets<-metab_data_sets
                                                                            randmetab_data_sets$XID<-metab_data_sets$XID[all]
                                                                            
                                                                            randmetab_data_sets$Statistic<-metab_data_sets$Statistic[all]
                                                                            
                                                                            randmetab_data_sets$Scaled.Statistic<-metab_data_sets$Scaled.Statistic[all]
                                                                              
                                                                            agg_res<-lapply(1:length(sid_list),function(c){
                                                                             
                                                                             x<-randmetab_data_sets$Statistic[which(metab_data_sets$SetID==sid_list[c])]
                                                                                                   
                                                                                                   x<-as.numeric(as.character(x))
                                                                                                   
                                                                                                  
                                                                                                   set_size=as.numeric(as.character(metab_data_sets$Total.Size[which(metab_data_sets$SetID==sid_list[c])]))
                                                                                                   
                                                                                                   agg.stat<-sum((x))/set_size
                                                                                                   
                                                                                                   agg.scaled.stat<-agg.stat
                                                                                                   
                                                                                                   z.score=sqrt(set_size)*agg.scaled.stat
                                                                                                   
                                                                                                   max_mean<-max(mean(x[which(x>0)],na.rm=TRUE),mean(x[which(x<0)],na.rm=TRUE),na.rm=TRUE)/set_size
                                                                                                   
                                                                                                   resv=cbind(as.character(sid_list[c]),round(agg.stat,3),round(agg.scaled.stat,3),round(z.score,3),round(max_mean,3),length(x),set_size,round(sum(abs(x)),3))
                                                                                                   resv=as.data.frame(resv)
                                                                                                   
                                                                                                   resv=unique(resv)
                                                    
                                                                    return(resv)
                                                                            })
                                                                            
                                                                            agg_res<-ldply(agg_res,rbind)
                                                                            return(agg_res)
                    
                                },metab_data_sets=metab_data_sets,sid_list=sid_list)
                                
                                # b=Sys.time()
                                #print(b-a)
                                stopCluster(cl)
                                
                                
                                
                                #  print("End here")
        
        randmetabset_res<-ldply(randmetabset_res,rbind) #do.call(data.frame,agg_res)
        colnames(randmetabset_res)=c("SetID","Agg.Statistic","Agg.Scaled.Statistic","Z.score","MaxMean","Num.Hits","Total.Size","Sum.Statistic")
        
        # save(randmetabset_res,file="r1.rda")
 
 randmetabset_res$MaxMean<-as.numeric(as.character(randmetabset_res$MaxMean))
 
 randmetabset_res[,-c(1)]<-apply(randmetabset_res[,-c(1)],2,function(x){as.numeric(as.character(x))})
 
 randmetabset_res_maxmean_avg<-mean(randmetabset_res$MaxMean,na.rm=TRUE)
 randmetabset_res_maxmean_std<-sd(randmetabset_res$MaxMean,na.rm=TRUE)
       
        #checkhere
                       
                              
                                    meanp.agg<-function (p)
                                    {
                                        keep <- (p >= 0) & (p <= 1)
                                        invalid <- sum(1L * keep) < 2
                                        if (invalid) {
                                            warning("Must have at least four valid p values")
                                            res <- list(z = NA_real_, p = NA_real_, validp = p[keep])
                                        }
                                        else {
                                            pi <- mean(p[keep])
                                            k <- length(p[keep])
                                            z <- (0.5 - pi) * sqrt(12 * k)
                                            if (k != length(p)) {
                                                warning("Some studies omitted")
                                            }
                                            res <- list(z = z, p = pnorm(z, lower.tail = FALSE),
                                                validp = p[keep])
                                        }
                                        #class(res) <- c("meanp", "metap")
                                        res
                                    }

                                    set_statistic1<-as.data.frame(set_statistic1)
                                    
                                    # save(set_statistic1,randmetabset_res,sid_list,file="r2.rda")
                                    
                                   
                                            metabset_res1<-lapply(1:length(sid_list),function(c){
            
         
                                            cur_sid_rand<-randmetabset_res[which(as.character(randmetabset_res$SetID)==sid_list[c]),]
                                            cur_sid_orig<-set_statistic1[which(as.character(set_statistic1$SetID)==sid_list[c]),]
                        
                        
                        if(nrow(cur_sid_orig)>0){
                                                    pval_agg.stat=length(which(as.numeric(as.character(cur_sid_rand[,2]))>as.numeric(as.character(cur_sid_orig[1,2]))))/itrs
                                                    
                                                    
                                                    
                                                    #maxmean z-transformation
                                                    z.score=(as.numeric(as.character(cur_sid_orig[1,5]))-mean(as.numeric(as.character(cur_sid_rand[,5])),na.rm=TRUE))/sd(as.numeric(as.character(cur_sid_rand[,5])),na.rm=TRUE)
                                                    
                                                    #randmetabset_res_maxmean_avg)/randmetabset_res_maxmean_std
                                                    
                                                    pval_zscore=round(pnorm(abs(z.score),lower.tail = FALSE),3) #length(which(as.numeric(as.character(cur_sid_rand[,3]))>as.numeric(as.character(cur_sid_orig[3]))))/itrs
                                                    pval_stdmaxmean=length(which(as.numeric(as.character(cur_sid_rand[,5]))>as.numeric(as.character(cur_sid_orig[1,5]))))/itrs
                                                    
                                                    # meanp_val=mean(c(pval_agg.stat,pval_zscore,pval_stdmaxmean),na.rm=TRUE)
                                                    #zval=(0.5-meanp_val)*sqrt(12*3)
                                                    
                                                    
                                                    # pval_meta <- pval_stdmaxmean
                                                    
                                                    pval_meta <-meanp.agg(c(pval_agg.stat,pval_stdmaxmean))$p
                                                    
                                                    #save(pval_meta,pval_agg.stat,pval_zscore,pval_stdmaxmean,file="pdebug.Rda")
                                                    #round(pnorm((zval),lower.tail = FALSE),3)
                                                    
                                                    pval_mat<-c(as.character(sid_list[c]),round(pval_agg.stat,3),round(pval_zscore,3),round(pval_stdmaxmean,3),round(pval_meta,3))
                                                    #   return(pval_mat)
                                                     }
                                           })
                                           
                                           metabset_res_pvalmat<-ldply(metabset_res1,rbind)
                                          # save(metabset_res_pvalmat,file="metabset_res_pvalmat.Rda")
                                           metabset_res_qvalmat<-apply(metabset_res_pvalmat[,-c(1)],2,function(x){p.adjust(as.numeric(as.character(x)),method="BH")})
                                           
                                           metabset_res_pvalmat<-cbind(metabset_res_pvalmat,round(metabset_res_qvalmat,3))
                                           
                                           colnames(metabset_res_pvalmat)<-c("SetID","pval.Agg.Statistic","pval.Z.score","pval.MaxMean","pval.meta","FDR.Agg.Statistic","FDR.Z.score","FDR.MaxMean","FDR.meta")
                                           
                                           res<-merge(set_statistic1,metabset_res_pvalmat,by="SetID")
                                           
                                           res<-merge(res,set_size,by="SetID")
                                           
                                       save(res,file="reschild2.Rda")
                                           
                                           pdf("QQplot.pdf")
                                          my.pvalues=as.numeric(as.character(res$pval.MaxMean))
                                          exp.pvalues<-(rank(my.pvalues, ties.method="first")+.5)/(length(my.pvalues)+1)
                                          plot(-log10(exp.pvalues), -log10(my.pvalues), asp=1,main="QQplot-maxmean")
                                          abline(0,1)
                                          dev.off()
                                          
                                          if(FALSE){
                                          par(mfrow=c(2,2))
                                           my.pvalues=as.numeric(as.character(res$pval.Agg.Statistic))
                                           exp.pvalues<-(rank(my.pvalues, ties.method="first")+.5)/(length(my.pvalues)+1)
                                           plot(-log10(exp.pvalues), -log10(my.pvalues), asp=1,main="QQplot-meanStatistic")
                                           abline(0,1)
                                           
                                           my.pvalues=as.numeric(as.character(res$pval.Z.score))
                                                                                exp.pvalues<-(rank(my.pvalues, ties.method="first")+.5)/(length(my.pvalues)+1)
                                                                                plot(-log10(exp.pvalues), -log10(my.pvalues), asp=1,main="QQplot-Zscore")
                                                                                abline(0,1)
                                                                                
                                        
                                        
                                        my.pvalues=as.numeric(as.character(res$pval.meta))
                                        exp.pvalues<-(rank(my.pvalues, ties.method="first")+.5)/(length(my.pvalues)+1)
                                        plot(-log10(exp.pvalues), -log10(my.pvalues), asp=1,main="QQplot-meta")
                                        abline(0,1)
                                        
                                           dev.off()
                                          }
                                           
    return(res)

}
}

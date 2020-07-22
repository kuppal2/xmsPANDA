get_fcs <-
function(target.data,target.data.annot=NA,kegg_species_code="hsa",database="pathway",reference_set=NA,type.statistic="pvalue",fcs.min.hits=2,itrs=100, numnodes=2){
    
    #metab_data: ID, Statistic
    #reference set: SetID, Name
    suppressMessages(library('KEGGREST'))
    
    count.unique.formula.overlapsize=TRUE
    dup.feature.check=TRUE
    
    if(is.na(target.data.annot)==FALSE){
        
        if(target.data.annot=="none"){
            
            target.data.annot=NA
        }
    }
    if(is.na(reference_set)==TRUE){
        
        if(database=="pathway"){
                    #homo sapiens: humans
                    if(kegg_species_code=="hsa"){
                        data(kegg_hsa)
                        g1=kegg_hsa
                        rm(kegg_hsa)
                    }else{
                        #Mus musculus: mouse
                        if(kegg_species_code=="mmu"){
                            data(kegg_mmu)
                            g1=kegg_mmu
                            rm(kegg_mmu)
                        }else{
                            
                            #Pan troglodytes: Chimpanzee
                            if(kegg_species_code=="ptr"){
                                data(kegg_ptr)
                                g1=kegg_ptr
                                rm(kegg_ptr)
                            }else{
                                
                                #Macaca mulatta: Rhesus monkey
                                if(kegg_species_code=="mcc"){
                                    data(kegg_mcc)
                                    g1=kegg_mcc
                                    rm(kegg_mcc)
                                }else{
                                    #Bos taurus: cow
                                    if(kegg_species_code=="bta"){
                                        data(kegg_bta)
                                        g1=kegg_bta
                                        rm(kegg_bta)
                                    }else{
                                        #Rattus norvegicus: rat
                                        if(kegg_species_code=="rno"){
                                            data(kegg_rno)
                                            g1=kegg_rno
                                            rm(kegg_rno)
                                        }else{
                                            
                                            #Danio rerio: Zebrafish
                                            if(kegg_species_code=="dre"){
                                                data(kegg_dre)
                                                g1=kegg_dre
                                                rm(kegg_dre)
                                            }else{
                                                
                                                #C. elegans: nematode
                                                if(kegg_species_code=="cel"){
                                                    data(kegg_cel)
                                                    g1=kegg_cel
                                                    rm(kegg_cel)
                                                }else{
                                                    
                                                    #Drosophila melanogaster: fruit fly
                                                    if(kegg_species_code=="dme"){
                                                        data(kegg_dme)
                                                        g1=kegg_dme
                                                        rm(kegg_dme)
                                                    }
                                                    
                                                }
                                                
                                            }
                                            
                                        }
                                        
                                    }
                                    
                                }
                            }
                            
                        }
                        
                    }
                   
        }else{
            
            if(database=="module"){
                
                data(kegg_modules)
                g1=kegg_modules
                rm(kegg_modules)
            }else{
                if(database=="brite"){
                    
                    data(kegg_brite)
                    g1=kegg_brite
                    rm(kegg_brite)
                }else{
                    if(database=="lipidmaps_mainclass"){
                        
                        data(lipidmaps_mainclass)
                        g1=lipidmaps_mainclass
                        rm(lipidmaps_mainclass)
                    }else{
                        
                        if(database=="lipidmaps_subclass"){
                            
                            data(lipidmaps_subclass)
                            g1=lipidmaps_subclass
                            rm(lipidmaps_subclass)
                        }else{
                            
                            if(database=="refmet_superclass"){
                                
                                data(refmet_superclass)
                                g1=refmet_superclass
                                rm(refmet_superclass)
                            }else{
                                if(database=="refmet_mainclass"){
                                    
                                    data(refmet_mainclass)
                                    g1=refmet_mainclass
                                    rm(refmet_mainclass)
                                }else{
                                    
                                    if(database=="refmet_subclass"){
                                        
                                        data(refmet_subclass)
                                        g1=refmet_subclass
                                        rm(refmet_subclass)
                                    }
                                }
                                
                            }
                        }
                    }
                    
                }
                
                
            }
            
            
        }
    }else{
        g1<-reference_set
        
        if(ncol(g1)<5){
            
            exact_mass<-seq(1,nrow(g1))
            chemformula<-paste("A",seq(1,nrow(g1)),sep="")
            
            
            g1<-cbind(g1,exact_mass,chemformula)
           
        
        }
        
        
    }
    g1<-as.data.frame(g1)
    #colnames(hsa_module_comp)<-c("MetabSetID","KEGGID")
    
    
    if(is.na(target.data.annot)==FALSE){
        
        cnames1<-colnames(target.data)
        
        if(grep("mz",cnames1)>0 && grep("time",cnames1)>0){
            
            #    res_1<-getVenn(dataA=target.data,name_a="target.data", dataB=target.data.annot,name_b="target.data.annot",mz.thresh=1,time.thresh=1,alignment.tool=NA, xMSanalyzer.outloc=getwd(),use.unique.mz=FALSE,plotvenn=FALSE,use.best.match=FALSE)
            
            target.data_1<-merge(target.data,target.data.annot,by=c("mz","time"))
            target.data_1<-target.data_1[,-c(1:2)]
            
        }else{
            
            
            if(grep("mz",cnames1)>0){
                
                target.data_1<-merge(target.data,target.data.annot,by=c("mz"))
                target.data_1<-target.data_1[,-c(1)]
                
            }else{
                
                if(grep("Name",cnames1)>0){
                    
                    target.data_1<-merge(target.data,target.data.annot,by=c("Name"))
                    target.data_1<-target.data_1[,-c(1)]
                }else{
                    
                    warning("Skipping merge with target annotation file as no columns are labeled as mz or Name in the annotation file.")
                }
            }
        }
        metab_data<-target.data_1
        rm(target.data_1)
        rm(target.data)
        
        #   colnames(metab_data)<-c("Statistic","mz","time")
        
        #colnames(metab_annot)<-c("XID","mz","time")
        
        #metab_data_1<-merge(metab_data,metab_annot,by="mz")
        metab_data<-metab_data[,c("XID","Statistic")]
        
    }else{
        
        colnames(target.data)<-c("XID","Statistic")
        metab_data<-target.data
        
        rm(target.data)
        
       
        metab_data<-metab_data[,c("XID","Statistic")]
    }
    
    if(type.statistic=="pvalue" || type.statistic=="p-value" || type.statistic=="p.value"){
        
       
          metab_data$Statistic=(-1)*log10(metab_data$Statistic)
           
    }
    
   colnames(g1)<-c("XID","SetID","Name","ExactMass","Formula")
   res<-get_fcs_child(metab_data=metab_data,reference_sets=g1,fcs.min.hits=fcs.min.hits,itrs=itrs,numnodes=numnodes)
   
   
   if(length(which(duplicated(g1$SetID)==TRUE))>0){
    path_names_ids<-g1[-which(duplicated(g1$SetID)==TRUE),]
   }else{
       path_names_ids<-g1
   }

   #   print(head(res))
    
    if(length(res)<1){
        res<-{}
        
        print("Functional class scoring returned no rsults. Please make sure that the features match between the feature table and the annotation files.")
    }else{
        res<-merge(res,path_names_ids[,c(1:3)],by.x="SetID",by.y="SetID")
    
    cnames1<-colnames(res)
    #print(head(res))
    cnames1[length(cnames1)]<-c("SetName")
    
    colnames(res)<-cnames1
    
    res<-res[order(as.numeric(as.character(res$FDR.meta)),decreasing=FALSE),]
    
   
    if(is.na(fcs.min.hits)==FALSE){
        res<-res[which(res$Num.Hits>=fcs.min.hits),]
     
    }
    
   
   if(length(res)>0){
       res<-res[,c("SetID","Agg.Statistic","Z.score","MaxMean","Total.Size","Num.Hits","pval.meta","FDR.meta","SetName")]
       
        res_temp<-merge(res,g1,by="SetID")
        res_temp<-merge(res_temp,metab_data,by="XID")
        
        res1<-aggregate(res_temp$XID,by=list(res_temp$SetID),function(x){paste(x,sep="",collapse=";")})
        res1<-do.call(data.frame,res1)
        colnames(res1)<-c("SetID","XID")
        
        res<-merge(res,res1,by="SetID")
         res<-res[order(as.numeric(as.character(res$FDR.meta)),decreasing=FALSE),]
   }else{
       
       res<-{}
   }
    }
    
    #print(head(res))
    #save(res,file="res.Rda")
    return(res)
    
}

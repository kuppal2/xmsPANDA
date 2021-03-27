get_fcs <-
function(target.data,target.data.annot=NA,kegg_species_code="hsa",database="pathway",
                  reference_set=NA,type.statistic=c("pvalue","t-statistic","fold.change","VIP"),fcs.min.hits=2,itrs=100, numnodes=2){
  
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  
  #metab_data: ID, Statistic
  #reference set: SetID, Name
  suppressMessages(library('KEGGREST'))
  
  colnames(target.data)<-c("XID","Statistic")
  
 # print(fcs.min.hits)
  #print(itrs)
 # print(numnodes)
  #print(head(target.data))
  #print(head(reference_set))
  #print(type.statistic)
  
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
                  }else{
                    if(database=="reactome_compound"){
                      
                      data(reactome_atlas)
                      g1<-reactome_atlas[which(reactome_atlas$IDtype=="compound"),]
                      rm(reactome_atlas)
                    }else{
                      if(database=="reactome_atlas"){
                        
                        data(reactome_atlas)
                        g1<-reactome_atlas
                        rm(reactome_atlas)
                      }else{
                        if(database=="kegg_atlas"){
                          #load("/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/xmsPANDA_v1.1.65/data/kegg_atlas.rda")
                           data(kegg_atlas)
                          g1<-kegg_atlas
                          rm(kegg_atlas)
                        }
                        
                      }
                      
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
    
    
    colnames(reference_set)<-c("XID","SetID","SetName")
    
    g1<-reference_set
    
    if(FALSE){
    if(ncol(g1)<5){
      
      
      exact_mass<-seq(1,nrow(g1))
      chemformula<-paste("A",seq(1,nrow(g1)),sep="")
      
      
      g1<-cbind(g1,exact_mass,chemformula)
      
      
    }
    
  }
    
  }
  g1<-as.data.frame(g1)
  #colnames(hsa_module_comp)<-c("MetabSetID","KEGGID")
  
  
  if(is.na(target.data.annot)==FALSE){
    
    cnames1<-colnames(target.data)
    
    if(grep("mz",cnames1)>0 && grep("time",cnames1)>0){
      
      
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
  
  
  
  cnames_g1<-colnames(g1)
  
  check_colnames<-length(which(cnames_g1%in%c("XID","SetID","SetName")))
  
  if(check_colnames<3){
    
    check_colnames<-length(which(cnames_g1%in%c("XID","SetID","Name")))
    
    if(check_colnames<3){
    
        stop("The columns in the reference database file should be: XID, SetID, and SetName")
    }
  }
  
  
  cnames_g1[1:3]<-c("XID","SetID","SetName")
  
  
  colnames(g1)<-cnames_g1 #c("XID","SetID","Name","ExactMass","Formula")
  res<-get_fcs_child(metab_data=metab_data,reference_sets=g1,fcs.min.hits=fcs.min.hits,itrs=itrs,numnodes=numnodes)
  
  save(res,g1,metab_data,file="r5.rda")
  
  g2=g1[which(g1$XID%in%metab_data$XID),]
  
  g1agg<-aggregate(g2$XID,by=list(g2$SetID),function(x){paste(x,sep="",collapse=";")})
  colnames(g1agg)<-c("SetID","XID")
  
  if(length(which(duplicated(g1$SetID)==TRUE))>0){
    path_names_ids<-g1[-which(duplicated(g1$SetID)==TRUE),c("SetID","SetName")]
  }else{
    path_names_ids<-g1[,c("SetID","SetName")]
  }
  
  #   print(head(res))
  res[,-c(1)]<-apply(res[,-c(1)],2,function(x){as.numeric(as.character(x))})
  
  
  if(length(res)<1){
    res<-{}
    
    print("Functional class scoring returned no rsults. Please make sure that the features match between the feature table and the annotation files.")
  }else{
    res<-merge(res,path_names_ids,by.x="SetID",by.y="SetID")
    
    cnames1<-colnames(res)
    #print(head(res))
    cnames1[length(cnames1)]<-c("SetName")
    
    colnames(res)<-cnames1
    
    res<-res[order(as.numeric(as.character(res$FDR.meta)),decreasing=FALSE),]
    
    
    if(is.na(fcs.min.hits)==FALSE){
      res<-res[which(res$Num.Hits>=fcs.min.hits),]
      
    }
    
    
    
    if(length(res)>0){
      
      res<-merge(res,g1agg,by="SetID")
      res<-res[order(as.numeric(as.character(res$pval.MaxMean)),decreasing=FALSE),]
      #  #save(res,file="res.Rda")
      #res<-res[,c("SetID","Agg.Statistic","Z.score","MaxMean","Total.Size","Num.Hits","pval.Agg.Statistic","pval.Z.score","pval.MaxMean","FDR.MaxMean","pval.meta","FDR.meta","SetName","XID")]
      # res<-res[,c("SetID","MaxMean","Total.Size.x","Num.Hits","pval.MaxMean","FDR.MaxMean","SetName","XID")]
      
      #   res<-res[which(abs(res$Z.score)>0),]
      
      # res<-res[,c("SetID","Agg.Statistic","Z.score","MaxMean","Total.Size","Num.Hits","pval.Agg.Statistic","pval.Z.score","pval.MaxMean","pval.meta","FDR.meta","SetName")]
      res<-res[,c("SetID","Total.Size","Num.Hits","Agg.Statistic","pval.Z.score","pval.MaxMean","FDR.Z.score","FDR.MaxMean","SetName","XID")]
      
      colnames(res)<-c("SetID","Total.Size","Num.Hits","Agg.Statistic","p.Z.stat","p.MaxMean","FDR.Z.stat","FDR.MaxMean","SetName","XID")
      
      #  res<-res[,c("SetID","MaxMean","MaxMean.Std","Total.Size","Num.Hits","pval.MaxMean","FDR.MaxMean","SetName","XID")]
      
      # colnames(res)<-res[,c("SetID","Agg.MaxMean.Statistic","Total.Size","Num.Hits","pval.Agg.Statistic","pval.Z.score","pval.MaxMean","pval.meta","FDR.meta","SetName","XID")]
      
      
    }else{
      
      res<-{}
    }
  }
  
  #print(head(res))
  ##save(res,file="res.Rda")
  return(res)
  
}

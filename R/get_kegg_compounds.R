get_kegg_compounds <-
function(kegg_species_code="hsa",database="pathway"){
  
  options(warn=-1)
  
  
  
  if(database=="pathway"){
    path_list<-keggList(kegg_species_code,database=database)
    
    kegg_pathway_ids<-names(path_list)
    kegg_pathway_ids<-gsub(kegg_pathway_ids,pattern="path:",replacement="")
  }else{
    
    if(database=="module"){
      
      path_list<-keggList(database=database)
      
      kegg_pathway_ids<-names(path_list)
      kegg_pathway_ids<-gsub(kegg_pathway_ids,pattern="md:",replacement="")
    }else{
      if(database=="brite"){
        path_list<-keggList(database="cpd")
        kegg_cpd_ids<-names(path_list) #[1:2]
        
        cl<-makeCluster(2)
        clusterExport(cl,"get_kegg_compounds_brite")
        clusterEvalQ(cl,library(KEGGREST))
        clusterEvalQ(cl,library(stringr))
        res<-parLapply(cl,1:length(kegg_cpd_ids),function(x,kegg_cpd_ids){
          
          Sys.sleep(0.5)
          keggid<-kegg_cpd_ids[x]
          restemp<-try(get_kegg_compounds_brite(keggid),silent=TRUE)
          Sys.sleep(0.5)
          if(is(restemp,"try-error")){
            # print(restemp)
            #print(keggid)
            #fname1<-paste("bad",keggid,".Rda",sep="")
            ##save(keggid,file=fname1)
          }else{
            return(restemp)
          }
        },kegg_cpd_ids=kegg_cpd_ids)
        
        res<-do.call(rbind,res)
        res<-as.data.frame(res)
        colnames(res)<-c("XID","SetID","Name","ExactMass","Formula")
        
        kegg_cpd_mat<-cbind(names(path_list),path_list)
        colnames(kegg_cpd_mat)<-c("KEGGID","Name")
        #  #save(kegg_cpd_mat,file="kegg_cpd_mat.Rda")
        stopCluster(cl)
        
        return(res)
        #kegg_pathway_ids<-gsub(kegg_pathway_ids,pattern="br:",replacement="")
      }else{
        
        if(database=="reaction"){
          path_list<-keggList(database=database)
          kegg_pathway_ids<-names(path_list)
          
          kegg_pathway_ids<-gsub(kegg_pathway_ids,pattern="rn:",replacement="")
        }else{
          
          if(database=="disease"){
            path_list<-keggList(database=database)
            kegg_pathway_ids<-names(path_list)
            
            kegg_pathway_ids<-gsub(kegg_pathway_ids,pattern="ds:",replacement="")
          }
        }
      }
      
    }
  }
  
  
  
  kegg_pathway_ids<-kegg_pathway_ids
  kegg_comp_list<-{}
  map_res<-{}
  kegg_module_list<-{}
  
  path_name_id_mapping<-cbind(kegg_pathway_ids,path_list)
  path_name_id_mapping<-as.data.frame(path_name_id_mapping)
  
  
  colnames(path_name_id_mapping)<-c("SetID","Name")
  
  
  
  path_comp_mat<-{}
  cl<-makeCluster(detectCores()*0.25)
  clusterExport(cl,"keggGet",envir = .GlobalEnv)
  clusterEvalQ(cl,library(KEGGREST))
  #
  path_comp_mat<-parLapply(cl,1:length(kegg_pathway_ids),function(p,kegg_pathway_ids)
  {
    kegg_pathway_id=kegg_pathway_ids[p]
    Sys.sleep(0.5)
    k1<-keggGet(dbentries=kegg_pathway_id)
    
    kegg_comp_list<-c(kegg_comp_list,k1[[1]]$COMPOUND)
    
    if(length(kegg_comp_list)<1){
      print("bad")
      print(kegg_pathway_id)
      kegg_module_list<-c(kegg_module_list,k1[[1]]$MODULE)
      
    }else{
      path_comp_mat_temp<-cbind(kegg_pathway_id,names(k1[[1]]$COMPOUND))
      
      if(length(path_comp_mat_temp)>0){
        
        if(ncol(path_comp_mat_temp)==2){
          
          return(path_comp_mat_temp)
        }
      }
    }
    Sys.sleep(0.05)
  },kegg_pathway_ids=kegg_pathway_ids)
  stopCluster(cl)
  path_comp_mat<-do.call(rbind,path_comp_mat)
  
  
  colnames(path_comp_mat)<-c("SetID","KEGGID")
  
  path_comp_mat<-merge(path_name_id_mapping,path_comp_mat,by.x="SetID",by.y="SetID")
  
  ##save(path_comp_mat,file="path_comp_mat.Rda")
  
  t1=table(path_comp_mat$SetID)
  
  
  # load("/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/xmsPANDA_v1.1.46/data/kegg_rno.rda")
  #path_comp_mat<-kegg_rno
  cl<-makeCluster(max(2,detectCores()*0.1))
  clusterExport(cl,"keggGet",envir = .GlobalEnv)
  clusterEvalQ(cl,library(KEGGREST))
  
  unique_keggIDs<-unique(path_comp_mat$KEGGID)
  
  #length(unique_keggIDs)
  kegg_comp_mat<-parLapply(cl,1:length(unique_keggIDs),function(p,unique_keggIDs)
  {
    Sys.sleep(0.5)
    k1<-keggGet(dbentries=as.character(unique_keggIDs[p]))
    
    res<-c(as.character(unique_keggIDs[p]),k1[[1]]$EXACT_MASS,k1[[1]]$FORMULA)
    Sys.sleep(0.5)
    return(res)
  },unique_keggIDs=unique_keggIDs)
  stopCluster(cl)
  
  
  kegg_comp_mat<-do.call(rbind,kegg_comp_mat)
  # #save(kegg_comp_mat,file="kegg_comp_mat.Rda")
  
  kegg_comp_mat<-as.data.frame(kegg_comp_mat)
  colnames(kegg_comp_mat)<-c("XID","ExactMass","Formula")
  path_comp_mat=merge(path_comp_mat,kegg_comp_mat,by="KEGGID")
  
  colnames(path_comp_mat)<-c("XID","SetID","Name","ExactMass","Formula")
  #kegg_rno<-path_comp_mat
  #    #save(kegg_rno,file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/xmsPANDA_v1.1.46/data/new/kegg_rno.Rda")
  
  #c("KEGGID","SetID","Name","ExactMass","Formula")
  
  
  options(warn=0)
  return(path_comp_mat)
}

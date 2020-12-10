get_kegg_compounds_brite <-
function(keggid) {
<<<<<<< HEAD
  Sys.sleep(0.5)
  b1 <- keggGet(keggid) #paste("cpd:", keggid, sep = ""))
  
  
  if(length(b1[[1]]$BRITE)>0){
=======
    Sys.sleep(0.5)
    b1 <- keggGet(keggid) #paste("cpd:", keggid, sep = ""))
    
    
    if(length(b1[[1]]$BRITE)>0){
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
    
    brite_inf<-b1[[1]]$BRITE
    
    
    brite_inf <- brite_inf[grep(brite_inf,pattern="BR:br")]
    
    brite_ids<-str_extract(brite_inf,pattern="BR:br[0-9]*")
    brite_name<-gsub(brite_inf,pattern=" [BR:br[0-9]*]",replacement="")
    
    exact_mass<-b1[[1]]$EXACT_MASS
    cformula<-b1[[1]]$FORMULA
    
    
    #otherdb_inf <- paste(b1[[1]]$DBLINKS, collapse = ";")
    r1 <- cbind(as.character(keggid),as.character(brite_ids), as.character(brite_name), exact_mass,cformula)
    colnames(r1)<-c("XID","SetID","Name","ExactMass","Formula")
    
    return(r1)
<<<<<<< HEAD
  }
  
  
=======
    }
    
    
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
}

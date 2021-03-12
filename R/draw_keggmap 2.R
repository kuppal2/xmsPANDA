draw_keggmap <-
function(targeted_table,outloc,foldchange.thresh=2,minhit=3,highcolor='red',lowcolor='blue',percent_node=0.2){
  
  suppressMessages(library(pathview))
  
  if(length(grep('keggid',colnames(targeted_table),ignore.case=TRUE))==0){
    stop("No KEGG ID was detected in reference list. To draw kegg map, please add one more column named 'KEGGID' to reference list for KEGG ID.")
  }
  
  iddata=targeted_table[,"KEGGID"]
  
  targeted_table[,'map'] <- unlist(mclapply(1:length(iddata), function(j) {
    API=paste("http://rest.kegg.jp/get/",as.character(iddata[j]),sep="")
    content<-try(readLines(API),silent=TRUE)
    if (class(content) == "try-error") {
      return(NA)
    }else{
      return(paste(gsub('^.*(map[0-9]*).*$','\\1',content[grep('map[0-9]',content)],perl=TRUE),collapse=","))
    }
  },mc.cores=detectCores()))
  
  totalmap=unique(unlist(strsplit(paste(targeted_table[,"map"][grep("map",targeted_table[,"map"])],collapse=","), ",")))
  mapdata=data.frame(mapID=totalmap)
  
  for(ii in 1:length(totalmap)){
    mapdata[ii,"ref_mz_time"]=paste(paste(round(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"ref_mz"],5),round(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"ref_time"],1),sep="_"),collapse="/")
    mapdata[ii,"mz_time"]=paste(paste(round(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"mz"],5),round(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"time"],1),sep="_"),collapse="/")
    mapdata[ii,"keggid"]=paste(as.character(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"KEGGID"]),collapse="/")
    mapdata[ii,"foldchange"]=paste(targeted_table[grep(totalmap[ii],targeted_table[,"map"]),"foldchange"],collapse="/")
  }
  
  write.table(mapdata, paste(outloc,"keggmapdata.txt",sep="/"), sep= "\t", row.names=FALSE)
  
  suppressWarnings(dir.create(paste(outloc,"KEGGmaps",sep="/"),,showWarnings = FALSE))
  
  print_png <- function(i,mapdata,foldchange.thresh,lowcolor,highcolor,minhit,outloc){
    
    #load the required packages
    library(pathview)
    
    if(as.character(mapdata[i,1])%in%c('map01100','map01110','map01130','map01120','map00121','map01060')){
      
    }else{
      cpd.data=as.numeric(unlist(strsplit(as.character(mapdata[i,"foldchange"]), "/")))
      names(cpd.data)=unlist(strsplit(as.character(mapdata[i,"keggid"]), "/"))
      cpd.data=cpd.data[!duplicated(cpd.data)]
      if(length(cpd.data)<minhit){
        
      }else{
        log.cpd.data <- log(cpd.data,base = foldchange.thresh)
        pathwayid=gsub("map","",as.character(mapdata[i,"mapID"]))
        suffix=paste("highlighted","_",length(log.cpd.data),'hits',sep="")
        setwd(paste(outloc,"KEGGmaps",sep="/"))
        try(pathview(gene.data =NULL, cpd.data=log.cpd.data, pathway.id=pathwayid, cpd.idtype='kegg', same.layer=TRUE, plot.col.key =TRUE,new.signature=FALSE,low = list(cpd = lowcolor), mid =list(cpd = "grey"),high=list(cpd = highcolor),out.suffix=suffix,species = "ko"),silent=TRUE)
        unlink(paste(paste("ko",pathwayid,sep=""),"png",sep="."), recursive = FALSE)
        unlink(paste(paste("ko",pathwayid,sep=""),"xml",sep="."), recursive = FALSE)
        print(getwd())
      }
    }
  }
  
  nu_cores <- detectCores()
  cl <- makeCluster(ceiling(nu_cores*as.numeric(percent_node)))
  parLapply(cl,1:nrow(mapdata),print_png,mapdata=mapdata,foldchange.thresh=foldchange.thresh,lowcolor=lowcolor,highcolor=highcolor,minhit=minhit,outloc=outloc)
  stopCluster(cl)
  
}

E.W.k <-
function(x, d.power=1,B=100) {
  
 
  rng.x1 <- apply(x, 2L, range)
  logWks <- matrix(0, B, 1)
 
  xs <- scale(x, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], 
                                               max = M[2]), nn = n)
  
  cl <- parallel::makeCluster(getOption("cl.cores", 0.5*detectCores()))
  clusterEvalQ(cl,library(WGCNA))
  clusterEvalQ(cl,library(fastcluster))
  clusterEvalQ(cl,library(dynamicTreeCut))
  clusterExport(cl,"W.k",envir = .GlobalEnv)
  res<-parLapply(cl,1:B,function(b){
    z <-x[sample(1:nrow(x)),sample(1:ncol(x))]#z1 + m.x
    
   
    dist2<-as.dist(1-WGCNA::cor(z))
    f2=fastcluster::hclust(d=dist2,method = "complete")
    
    mycl_metabs2 <-suppressWarnings(cutreeDynamic(f2,distM=as.matrix((dist2)),cutHeight = 0.95,verbose = FALSE,
                                                  ,deepSplit = 4,minClusterSize = 2,pamStage = FALSE,pamRespectsDendro = FALSE))
    names(mycl_metabs2)<-f2$labels
    
    wk<-W.k(x=t(z),clus=mycl_metabs2)
    #if(wk>0){
    cval1<- wk
    #}else{
      
     # cval1<- log(0.001)
    #}
    return(cval1)
  })
  stopCluster(cl)
  res<-unlist(res)
  print(summary(res,na.rm=TRUE))
  E.logW<-(mean(res,na.rm=TRUE))
  return(E.logW)
}

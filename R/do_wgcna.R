do_wgcna <-
function(feature_table_file=NA,class_labels_file=NA,X=NA,Y=NA,sigfeats=NA){
  
  if(typeof(X)=="logical"){
    data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
    
  }else{
    
    data_matrix<-X
  }
  if(typeof(Y)=="logical"){
    classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
    
  }else{
    
    classlabels<-Y
  }
  
  if(dim(classlabels)[2]>2){
    classgroup<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
  }else{
    
    classgroup<-classlabels[,2]
  }
  classlabels<-as.data.frame(classlabels)
  
  class_labels_levels<-levels(as.factor(classgroup))
  
  rnames<-paste(data_matrix$mz,data_matrix$time,sep="_") #paste(sprintf("%.4f",data_matrix$mz),data_matrix$time,sep="_")
  
  data_matrix_orig<-data_matrix
  data_matrix<-data_matrix[,-c(1:2)]
  data_matrix<-na.omit(data_matrix)
  
  rnamesAB<-gsub(pattern="NA_NA",replacement=NA,x=rnames)
  rnamesAB<-na.omit(rnamesAB)
  
  nSets = length(class_labels_levels);
  multiExpr = vector(mode = "list", length = nSets)
  data_matrix_list<-new("list")
  num_samps_groups<-new("list")
  for(i in 1:length(class_labels_levels)){
    
    data_matrix_list[[i]]<-t(data_matrix[,which(classgroup==class_labels_levels[i])])
    num_samps_groups[[i]]<-dim(data_matrix_list[[i]])[1]
    #print(dim(data_matrix_list[[i]]))
    multiExpr[[i]]<-list(data = as.data.frame(data_matrix_list[[i]]));
    rownames(multiExpr[[i]]$data)=c(paste(rep(class_labels_levels[i],num_samps_groups[[i]]),seq(1,num_samps_groups[[i]]),sep=""))
  }
  
  data_matrix_all<-t(data_matrix)
  
  
  # We work with two sets:
  # For easier labeling of plots, create a vector holding descriptive names of the two sets.
  setLabels = as.character(class_labels_levels) #c("Slow", "Rapid")
  shortLabels = as.character(class_labels_levels) #c("Slow", "Rapid")
  
  ####saveclass_labels_levels,file="class_labels_levels.Rda")
  ####savemultiExpr,file="multiExpr.Rda")
  
  exprSize = checkSets(multiExpr)
  
  
  sampleTrees = list()
  for (set in 1:nSets)
  {
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
  }
  
  
  #pdf(file = "SampleClustering.pdf", width = 12, height = 12);
  
  
  for (set in 1:nSets)
  {
    # Find clusters cut by the line
    #labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
    # Keep the largest one (labeled by the number 1)
    #keep = (labels==1)
    multiExpr[[set]]$data = multiExpr[[set]]$data
  }
  collectGarbage();
  
  
  # Check the size of the leftover data
  exprSize = checkSets(multiExpr)
  exprSize;
  
  traitData<-as.data.frame(classlabels)
  
  dim(traitData)
  names(traitData)
  
  # See how big the traits are and what are the trait and sample names
  # Form a multi-set structure that will hold the clinical traits.
  Traits = vector(mode="list", length = nSets);
  for (set in 1:nSets)
  {
    
    Traits[[set]] = list(data = data_matrix_list[[i]] );
    rownames(Traits[[1]]$data) = rownames(data_matrix_list[[i]]);
    #Traits[[2]] = list(data = dataexprB );
    #rownames(Traits[[2]]$data) = rownames(dataexprB);
  }
  collectGarbage();
  # Define data set dimensions
  nGenes = exprSize$nGenes;
  nSamples = exprSize$nSamples;
  
  ####savemultiExpr, Traits, nGenes, nSamples, file = "Consensus-dataInput.RData");
  
  
  
  
  
  
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(data=data_matrix_all, dataIsExpr=TRUE,powerVector = powers, verbose = 0)
  power_val=sft$powerEstimate
  
  if(is.na(power_val)==TRUE){
    power_val=6
  }
  
  netclassA = suppressWarnings(blockwiseConsensusModules(multiExpr, power = power_val, minModuleSize = 10, deepSplit = 2,
                                        pamRespectsDendro = FALSE,
                                        mergeCutHeight = 0.2, numericLabels = TRUE,saveTOMs = FALSE, verbose = 0,saveIndividualTOMs=FALSE,useDiskCache=FALSE))
  
  
  # open a graphics window
  sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(netclassA$colors)
  if(FALSE){
    # Plot the dendrogram and the module colors underneath
    plotDendroAndColors(netclassA$dendrograms[[1]], mergedColors[netclassA$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  }
  
  classAmoduleLabels = netclassA$colors
  classAmoduleColors = labels2colors(netclassA$colors)
  classAMEs = netclassA$MEs;
  classAgeneTree = netclassA$dendrograms[[1]];
  
  
  
  
  # Recalculate consMEs to give them color names
  consMEsC = multiSetMEs(multiExpr, universalColors = classAmoduleColors);
  
  MET = consensusOrderMEs(consMEsC);
  suppressWarnings(dir.create("Figures",showWarnings = FALSE))
  suppressWarnings(dir.create("Tables",showWarnings = FALSE))
  #setwd("NetworkAnalysis/")
  save(MET,setLabels,file="WGCNA_modules.Rda")
#  save(setLabels,file="setLabels.Rda")
  pdf("Figures/Module_preservation_analysis.pdf",width=10,height=8)
  #tiff("module_preservation.tiff",width=2000,height=2000)
  #sizeGrWindow(8,10);
  #par(cex = 1)
  #plotEigengeneNetworks_custom(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),zlimPreservation = c(0.5, 1), plotPreservation = "standard",plotDendrograms=FALSE)
  #dev.off()
  ##save(MET,setLabels,multiExpr,classAmoduleColors,consMEsC,file="wgcnamp.Rda")
  
  #mp = modulePreservation(MET, setLabels,dataIsExpr = TRUE)
  ##save(mp,file="mp.Rda")
  
  #sizeGrWindow(8,10);
  #tiff("module_preservation.tiff",width=2000,height=2000)
  #par(cex = 0.9)
  plotEigengeneNetworks_custom(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),zlimPreservation = c(0.5, 1), plotPreservation = "standard",
                               plotDendrograms=FALSE)
  try(dev.off(),silent=TRUE)
  graphics.off()
  graphics.off()
  
  
  colors_met<-cbind(data_matrix_orig[,c(1:2)],MET[[1]]$validColors)
  
  colors_met<-as.data.frame(colors_met)
  
  colnames(colors_met)<-c("mz","time","Module")
  
  #print(table(MET[[1]]$validColors))
  #print(table(classAmoduleColors))
  write.table(colors_met,file="Tables/Allfeature_modules.txt",sep="\t",row.names=FALSE)
  
  if(is.na(sigfeats)==FALSE){
    sub_colors_met<-colors_met[which(data_matrix_orig$mz%in%sigfeats$mz),]
    write.table(sub_colors_met,file="Tables/Sigfeature_modules.txt",sep="\t",row.names=FALSE)
  }
  
  ####saveMET,file="MET.Rda")
  
  
  
  return(list(MET=MET,setLabels=setLabels))
}

plotEigengeneNetworks_custom <-
function (multiME, setLabels, letterSubPlots = FALSE, Letters = NULL,
                                        excludeGrey = TRUE, greyLabel = "grey", plotDendrograms = TRUE,
                                        plotHeatmaps = TRUE, setMargins = TRUE, marDendro = NULL,
                                        marHeatmap = NULL, colorLabels = TRUE, signed = TRUE, heatmapColors = NULL,
                                        plotAdjacency = TRUE, printAdjacency = FALSE, cex.adjacency = 0.8,
                                        coloredBarplot = TRUE, barplotMeans = TRUE, barplotErrors = FALSE,
                                        plotPreservation = "standard", zlimPreservation = c(0, 1),
                                        printPreservation = FALSE, cex.preservation = 0.8)
{
  size = checkSets(multiME, checkStructure = TRUE)
  if (!size$structureOK) {
    multiME = fixDataStructure(multiME)
  }
  if (is.null(Letters))
    Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  if (is.null(heatmapColors))
    if (signed) {
     # heatmapColors = greenWhiteRed(50)
    }
  else {
    #heatmapColors = topo.colors(30) #heat.colors(30)
  }
  heatmapColors = topo.colors(30)
  
  heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
  heatmap_cols<-rev(heatmap_cols)
  
  nSets = length(multiME)
  
  cex.before <- par("cex")
  cex = par("cex")
  mar = par("mar")
  nPlotCols = 1 #nSets
  nPlotRows = 1 #as.numeric(plotDendrograms) + nSets * as.numeric(plotHeatmaps)
  if (nPlotRows == 0)
    stop("Nothing to plot: neither dendrograms not heatmaps requested.")
 par(mfrow = c(nPlotRows, nPlotCols))
  par(cex = cex)
  if (excludeGrey)
    for (set in 1:nSets) multiME[[set]]$data = multiME[[set]]$data[,
                                                                   substring(names(multiME[[set]]$data), 3) != greyLabel]
  plotPresTypes = c("standard", "hyperbolic", "both","change","differences")
  ipp = pmatch(plotPreservation, plotPresTypes)
  if (is.na(ipp))
    stop(paste("Invalid 'plotPreservation'. Available choices are",
               paste(plotPresTypes, sep = ", ")))
  letter.ind = 1

#  save(multiME,nSets,setLabels,colorLabels,heatmap_cols,file="debugwgcna.rda")
  if (plotHeatmaps)
    for (i.row in (1:nSets)) for (i.col in (i.row:nSets)) {
      letter.ind = i.row * nSets + i.col
      if (letterSubPlots) {
        letter = paste(substring(Letters, first = letter.ind,
                                 last = letter.ind), ".  ", sep = "")
      }
      else {
        letter = NULL
      }
      par(cex = cex)
      if (setMargins) {
        if (is.null(marHeatmap)) {
          if (colorLabels) {
            par(mar = c(1, 2, 3, 4) + 0.2)
          }
          else {
            par(mar = c(6, 7, 3, 5) + 0.2)
          }
        }
        else {
          par(mar = marHeatmap)
        }
      }
      nModules = dim(multiME[[i.col]]$data)[2]
      textMat = NULL
      
     
      par(cex = 0.7)
      
    
      if (i.row == i.col) {
        corME = WGCNA::cor(multiME[[i.col]]$data, use = "p")
        pME = corPvalueFisher(corME, nrow(multiME[[i.col]]$data))
        if (printAdjacency) {
          textMat = paste(signif(corME, 2), "\n", signif(pME,1))
          dim(textMat) = dim(corME)
        }
        
        corrplot::corrplot(corME,col=heatmap_cols,cl.cex=0.8,tl.cex=0.45,title=setLabels[i.col],cl.offset=2,tl.offset=0.2,mar=c(1, 2, 4, 2) + 0.01)
      }
      else {
        corME1 = WGCNA::cor(multiME[[i.col]]$data, use = "p")
        corME2 = WGCNA::cor(multiME[[i.row]]$data, use = "p")
        cor.dif = (corME1 - corME2)/2
        
        corrplot::corrplot(cor.dif,col=heatmap_cols,cl.cex=0.9,tl.cex=0.45,title=paste("Delta correlation: ",setLabels[i.row], " vs ", setLabels[i.col],sep=""),
                           cl.offset=2,tl.offset=0.2,mar=c(1, 2, 4, 2) + 0.01)
        
     #   save(corME1,corME2,i.row,i.col,multiME,setLabels,file="debugbarplot.Rda")
        
        d = tanh((corME1 - corME2)/(abs(corME1) + abs(corME2))^2)
       
              dp = 1 - abs(cor.dif)
              method = "Preservation"
          
          diag(dp) = 0
          #write.table(dp,file="preservation_matrix.txt",sep="\t")
         
            sum_dp = mean(dp[upper.tri(dp)])
            means = apply(dp, 2, sum)/(ncol(dp) - 1)
            if (barplotErrors) {
              errors = sqrt((apply(dp^2, 2, sum)/(ncol(dp) -1) - means^2)/(ncol(dp) - 2))
            }
            else {
              errors = NULL
            }
            Dmatrix=cbind(names(multiME[[i.col]]$data),means)
            Dmatrix=as.data.frame(Dmatrix)
            colnames(Dmatrix)<-c("Module","meanPreservationScore")
            
            fname<-paste("Tables/mean_preservation_matrix_",i.row,"vs",i.col,".txt",sep="")
            write.table(Dmatrix,file=fname,sep="\t",row.names=FALSE)
            labeledBarplot(means, names(multiME[[i.col]]$data),
                           main = paste(method, " score for ",i.row," vs ", i.col, ":",signif(sum_dp,
                                                                                                      2)), ylim = c(0, 1), 
                           stdErrors = errors,cex.lab=cex.preservation,cex.main=cex.preservation)
            
                           #colorLabels = colorLabels,
               #            colored = coloredBarplot, setStdMargins = FALSE,
                          
          # barplot(as.numeric(as.character(Dmatrix$meanPreservationScore)),Dmatrix$Module,main = paste("Preservation score for ",setLabesl[i.row]," vs ", setLabels[i.col], ":",signif(sum_dp,2),sep=""),cex.main=0.7,ylim = c(0, 1))
       
        }
      
    }
  
  par(cex = cex.before)
}

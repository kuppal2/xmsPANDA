do_vst_norm_count <-
function(countdata){
  library(DESeq)
  condition <- factor(rep("Class", ncol(countdata)))
  countdata <- newCountDataSet(countdata,condition )
  countdata <- estimateSizeFactors( countdata )
  cdsBlind <- DESeq::estimateDispersions( countdata, method="blind")
  vstdata <- varianceStabilizingTransformation( cdsBlind )
  return(exprs(vstdata))
}

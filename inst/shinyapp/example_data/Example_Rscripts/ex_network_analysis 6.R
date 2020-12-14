

library(xmsPANDA)



###Change file locations#######
feature_table_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/xmsPANDA_v1.1.57/examples_and_manual/Example_feature_table_and_classlabels/H1N1/exh1n1_network_analysis_feature_table.txt"
parentoutput_dir="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1out/"
class_labels_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_classlabels.txt"

sig_feat_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/xmsPANDA_v1.1.57/examples_and_manual/Example_feature_table_and_classlabels/H1N1/exh1n1_network_analysis_selected_feature_table.txt"

################################

net_res<-metabnet(feature_table_file=feature_table_file,
target.metab.file=NA, sig.metab.file=sig_feat_file,
parentoutput_dir=parentoutput_dir,
class_labels=NA,num_replicates=1,cor.method="spearman",abs.cor.thresh=0.4,cor.fdrthresh=0.2,target.mzmatch.diff=10,
target.rtmatch.diff=NA,max.cor.num=100,feat.filt.thresh=NA,summarize.replicates=FALSE,group.missing.thresh=0.7,
log2transform=TRUE,medcenter=FALSE,znormtransform=TRUE,quantile_norm=TRUE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0, 
networktype="complete", summary.na.replacement="zeros",samplermindex=NA,net_node_colors=c("green","red"), net_legend=FALSE)


names(net_res)

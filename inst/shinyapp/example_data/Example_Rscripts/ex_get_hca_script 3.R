
library(xmsPANDA)


feature_table_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_metabolome.txt"
parentoutput_dir="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1out/"
class_labels_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_classlabels.txt"



get_hca(feature_table_file,parentoutput_dir,class_labels_file,heatmap.col.opt="RdBu",cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
sample.col.opt="journal")

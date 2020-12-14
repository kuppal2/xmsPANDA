

library(xmsPANDA)

avg_data<-data_preprocess(feature_table_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_metabolome.txt",
parentoutput_dir="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1out/",
class_labels_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_classlabels.txt",
num_replicates=1,
rep.max.missing.thresh=0.3,
summarize.replicates=TRUE,summary.method="median",
 all.missing.thresh=0.5,group.missing.thresh=0.8,
normalization.method = "log2transform",
missing.val=0,
summary.na.replacement="knn",
pairedanalysis = FALSE,
input.intensity.scale = "raw",
create.new.folder = TRUE
)

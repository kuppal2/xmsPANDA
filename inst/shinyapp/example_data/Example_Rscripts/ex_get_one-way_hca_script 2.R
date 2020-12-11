
library(xmsPANDA)


feature_table_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/lcms_table.txt"
parentoutput_dir="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/testlive_hca2/"
class_labels_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/classlabels_gender.txt"

heatmap.col.opt="RdBu"

cor.method="spearman"
plots.width=2000
plots.height=2000
plots.res=300
is.data.znorm=FALSE
analysismode="classification"
sample.col.opt="rainbow"


get_hca(feature_table_file,parentoutput_dir,class_labels_file,heatmap.col.opt="RdBu",cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300,alphacol=0.3,hca_type="one-way")


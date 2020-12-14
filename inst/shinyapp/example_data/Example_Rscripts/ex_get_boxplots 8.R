
library(xmsPANDA)


feature_table_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_metabolome.txt"
parentoutput_dir="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1out/"
class_labels_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_classlabels.txt"


boxplot.col.opt="journal"
boxplot.type="ggplot" #Other options: "simple"
filename="H1N1boxplots"

get_boxplots(feature_table_file=feature_table_file,X=NA,Y=NA,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,boxplot.col.opt=boxplot.col.opt,boxplot.type=boxplot.type,filename=filename,newdevice=TRUE)

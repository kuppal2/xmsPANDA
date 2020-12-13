
library(xmsPANDA)


##training set files###
feature_table_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_metabolome.txt"
class_labels_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_classlabels.txt"
###################

##test set files### Same files are used in this example
testfeature_table_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_metabolome.txt"
testclass_labels_file="/Users/karanuppal/Downloads/examples_and_manual 2/Example_feature_table_and_classlabels/H1N1/exh1n1_classlabels.txt"
###############

##Read files####
featuretable<-read.table(feature_table_file,sep="\t",header=TRUE)
testfeaturetable<-read.table(testfeature_table_file,sep="\t",header=TRUE)

classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
testclasslabels<-read.table(testclass_labels_file,sep="\t",header=TRUE)

classlabels<-classlabels[,2] #column with class information
testclasslabels<-testclasslabels[,2] #column with class information

kfold=10 #dim(featuretable[,-c(1:2)])[2] #LOOCV recommended for only small studies N<30; or can be set to 5 or 10 for k-fold CV

errortype="BAR"; #other options: "AUC", "total"


#pdf("testacc_carnitines.pdf")

res<-get_classification.accuracy(kfold=kfold,featuretable=featuretable,classlabels=classlabels,classifier="svm",testfeaturetable=testfeaturetable,testclasslabels=testclasslabels,errortype="BAR",kernelname="radial")

print(res)

#dev.off()



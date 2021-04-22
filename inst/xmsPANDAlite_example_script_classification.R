#Example script to run diffexp.lite() to compare categorical groups "classification" mode
#load xmsPANDA
suppressPackageStartupMessages(library(xmsPANDA))

#change the input and output locations
feature_table_file<-"~/exh1n1_metabolome.txt"
class_labels_file<-"~/exh1n1_classlabels.txt"
outloc<-"~/xmsPANDAliteout/"


Xmat<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
Ymat<-read.table(class_labels_file,sep="\t",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)

#limma
demetabs_reslite1<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="limma",normalization.method = "log2transform",
                                pvalue.thresh = 0.05,fdrthresh=0.1,fdrmethod="BH",foldchangethresh = 0)

#pls
demetabs_reslite2<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="pls",normalization.method = "log2transform",pls_vip_thresh=2,foldchangethresh = 0)


#find common features selected by limma and pls
demetabs_reslite_1and2<-merge(demetabs_reslite1$diffexp_metabs[,c("Name","P.value","adjusted.P.value")],demetabs_reslite2$diffexp_metabs,by=c("Name"))

#pamr
demetabs_reslite3<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="pamr",normalization.method = "log2transform",fdrthresh=0.05)

#rf
demetabs_reslite4<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="rf",normalization.method = "log2transform")

#lmreg
demetabs_reslite5<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="lmreg",normalization.method = "log2transform",fdrthresh=0.05,fdrmethod="none")

#lm1wayanova
demetabs_reslite6<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="lm1wayanova",normalization.method = "log2transform",fdrthresh=0.05,fdrmethod="none")

#wilcox
demetabs_reslite7<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="wilcox",normalization.method = "log2transform",fdrthresh=0.05,fdrmethod="none")


#TTEST
demetabs_reslite8<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="ttest",normalization.method = "log2transform",fdrthresh=0.05,fdrmethod="none")



#####################################################


####################################################################
#Options for featselmethod:
#"limma": for one-way ANOVA using LIMMA (mode=classification)
#"limma2way": for two-way ANOVA using LIMMA (mode=classification)
#"limma1wayrepeat": for one-way ANOVA repeated measures using LIMMA (mode=classification)
#"limma2wayrepeat": for two-way ANOVA repeated measures using LIMMA (mode=classification)
#"lm1wayanova": for one-way ANOVA using linear model (mode=classification)
#"lm2wayanova": for two-way ANOVA using linear model (mode=classification)
#"lm1wayanovarepeat": for one-way ANOVA repeated measures using linear model (mode=classification)
#"lm2wayanovarepeat": for two-way ANOVA repeated measures using linear model (mode=classification)
#"lmreg": variable selection based on p-values calculated using a linear regression model; 
#allows adjustment for covariates (mode= regression or classification)
#"logitreg": variable selection based on p-values calculated using a logistic regression model; 
# allows adjustment for covariates (mode= classification)
#"rfesvm": uses recursive feature elimination SVM algorithm for variable selection; 
#(mode=classification)
#"wilcox": uses Wilcoxon tests for variable selection; 
#(mode=classification)
#"RF": for random forest based feature selection (mode= regression or classification)
#"RFconditional": for conditional random forest based feature selection (mode= regression or classification)
#"pamr": for prediction analysis for microarrays algorithm based on the nearest shrunked centroid method (mode=classification)
#"MARS": for multiple adaptive regression splines (MARS) based feature selection
#(mode= regression or classification)
#"pls": for partial least squares (PLS) based feature selection
#(mode= regression or classification)
#"spls": for sparse partial least squares (PLS) based feature selection
#(mode= regression or classification)
#"o1pls": for orthogonal partial least squares (OPLS) based feature selection
#(mode= regression or classification)
####################################################################

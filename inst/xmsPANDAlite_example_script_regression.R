#Example script to run diffexp.lite() to compare categorical groups "classification" mode
#load xmsPANDA
suppressPackageStartupMessages(library(xmsPANDA))

#change the input and output locations
feature_table_file<-"~/exh1n1_metabolome.txt"
class_labels_file<-"~/exh1n1_regression.txt"
outloc<-"~/xmsPANDAliteoutreg/"

Xmat<-read.table(feature_table_file,sep="\t",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
Ymat<-read.table(class_labels_file,sep="\t",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)

#lmreg
demetabs_reslite1<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="lmreg",normalization.method = "log2transform",
                                pvalue.thresh = 0.05,fdrthresh=0.1,fdrmethod="BH",analysismode="regression")

#pls
demetabs_reslite2<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="pls",normalization.method = "log2transform",
                                vipthresh=2,analysismode="regression")

#rf
demetabs_reslite3<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="rf",normalization.method = "log2transform",analysismode="regression")

#lmregrobust-HC3
demetabs_reslite4<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="lmregrobust",
                                normalization.method = "log2transform",analysismode="regression",vcovHC.type="HC3")

#lmregrobust
demetabs_reslite5<-diffexp.lite(Xmat=Xmat,Ymat=Ymat,outloc=outloc,featselmethod="lmregrobust",
                                normalization.method = "log2transform",analysismode="regression",vcovHC.type="HC0")


#Venn diagram showing overlap between different methods
input.i<-list(lmreg=demetabs_reslite1$diffexp_metabs$Name,
              pls=demetabs_reslite2$diffexp_metabs$Name,
              rf=demetabs_reslite3$diffexp_metabs$Name,
              lmregrobustHC3=demetabs_reslite4$diffexp_metabs$Name,
              lmregrobustHC0=demetabs_reslite5$diffexp_metabs$Name)
venn(input.i)


#####################################################


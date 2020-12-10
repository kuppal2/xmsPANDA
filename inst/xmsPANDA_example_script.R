
<<<<<<< HEAD
#load xmsPANDA
library(xmsPANDA)

#change the input and output locations
feature_table_file<-"/Users/karanuppal/Desktop/H1N1/exh1n1_metabolome.txt"
=======
#load xmsPANDA-v1.1
library(xmsPANDA)

source("/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/xmsPANDA_v1.1.67/inst/shinyapp/R/source_codes/xmsPANDA_v1.0.8.52.R")

#change the input and output locations

feature_table_file<-"/Users/karanuppal/Desktop/H1N1/exh1n1_metabolome.txt"

>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
class_labels_file<-"/Users/karanuppal/Desktop/H1N1/exh1n1_classlabels.txt"
outloc<-"/Users/karanuppal/Desktop/H1N1/testlog2knn4/"


#start: see manual for additional arguments and description
demetabs_res<-diffexp(
        #1) arguments for input files
        feature_table_file=feature_table_file,
        parentoutput_dir=outloc,
        class_labels_file=class_labels_file,
        input.intensity.scale="raw",

        ##2) data preprocessing order: 1) summarization, 2) filtering by missing values, 3) imputation; 4) transformation and normalization
        num_replicates = 1,
        summarize.replicates =TRUE, summary.method="median",summary.na.replacement="knn",
        rep.max.missing.thresh=0.3,
        all.missing.thresh=0.1, group.missing.thresh=0.8, missing.val=0,
        log2transform = FALSE, medcenter=FALSE, znormtransform = FALSE,
        quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE,
        TIC_norm=FALSE,
        rsd.filt.list = c(1),normalization.method="autoscaling",

        ##3) options for feature seletion: "limma","ttest","wilcox","lm1wayanova","lmreg","pls",
        #"pamr","spls","pls","o1pls","MARS","RF","rfesvm","logitreg", "poissonreg",
        #"ttestrepeat","wilcoxrepeat", "lm1wayanovarepeat","limma1wayrepeat","spls1wayrepeat"
        #"lm2wayanova","lm2wayanovarepeat","limma2way","limma2wayrepeat","spls2wayrepeat"
        pairedanalysis = FALSE, featselmethod=c("limma"),
        pvalue.thresh=0.05,
        fdrthresh = 0.05, fdrmethod="BH",
        kfold=5,networktype="complete",
        samplermindex=NA,numtrees=5000,analysismode="classification",pls_vip_thresh = 2, num_nodes = 3,
        max_varsel = 100, pls_ncomp = 5,pred.eval.method="BER",rocfeatlist=seq(2,10,1),
        rocfeatincrement=TRUE,
        rocclassifier="svm",foldchangethresh=1,
        optselect=TRUE,max_comp_sel=1,saveRda=FALSE,pls.permut.count=NA,
        pca.ellipse=TRUE,ellipse.conf.level=0.95,svm.acc.tolerance=5,pamr.threshold.select.max=FALSE,
        aggregation.method="none",mars.gcv.thresh=50,

        #4) arguments for centrality analysis, WGCNA and global clustering analysis (HCA and EM clustering)
        degree_rank_method=NA,differential.network.analysis=TRUE, wgcnarsdthresh=1,WGCNAmodules=TRUE,globalclustering=TRUE,

        #5) arguments for correlation and network analysis using the selected features
        cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
        globalcor=TRUE,target.metab.file=NA,
        target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,

        #6) arguments for graphical options: see manual for additional arguments
        output.device.type="png",pca.cex.val=4,legendlocation="bottomleft",
        net_node_colors=c("green","red"),
        net_legend=FALSE,aggregation.max.iter=100,
        heatmap.col.opt="redblue",manhattanplot.col.opt=c("darkblue","red3"),
        boxplot.col.opt=c("journal"),barplot.col.opt=c("journal"),individualsampleplot.col.opt="journal",
        lineplot.col.opt="journal",hca_type="two-way",cex.plots=0.6,
        lineplot.lty.option=c("dotted", "solid", "dashed", "dotdash", "longdash", "twodash"),
<<<<<<< HEAD
        timeseries.lineplots=FALSE,lme.modeltype="RI",ylab_text="Intensity",boxplot.type="simple",
        multiple.figures.perpanel = FALSE,fill.plots=TRUE
=======
        timeseries.lineplots=FALSE,lme.modeltype="RI",ylab_text="Normalized intensity",boxplot.type="simple"

>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3
)
sink(file=NULL)
#end
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

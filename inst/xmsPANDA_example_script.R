
#load xmsPANDA
suppressPackageStartupMessages(library(xmsPANDA))

#change the input and output locations
feature_table_file<-"~/exh1n1_metabolome.txt"
class_labels_file<-"~/exh1n1_classlabels.txt"
outloc<-"~/test/"


#start: see manual for additional arguments and description
demetabs_res<-diffexp(
        #1) arguments for input files
        feature_table_file=feature_table_file,
        parentoutput_dir=outloc,
        class_labels_file=class_labels_file,
        input.intensity.scale="raw",

        ##2) data preprocessing order: 1) summarization, 2) filtering by missing values, 3) imputation; 4) transformation and normalization
        #options for normalization methods: log2quantilenorm, log2transform, znormtransform, lowess_norm, quantile_norm, rangescaling, paretoscaling, mstus, eigenms_norm,
        #vsn_norm, sva_norm, none, tic_norm, cubicspline_norm, mad_norm
        num_replicates = 1,
        summarize.replicates =TRUE, summary.method="median",summary.na.replacement="knn",
        rep.max.missing.thresh=0.3,
        all.missing.thresh=0.1, group.missing.thresh=0.8, missing.val=0,
        rsd.filt.list = c(1),
        normalization.method="log2transform",

        ##3) options for feature seletion: "limma","ttest","wilcox","lm1wayanova","lmreg","pls",
        #"pamr","spls","pls","o1pls","MARS","RF","rfesvm","logitreg", "poissonreg",
        #"ttestrepeat","wilcoxrepeat", "lm1wayanovarepeat","limma1wayrepeat","spls1wayrepeat"
        #"lm2wayanova","lm2wayanovarepeat","limma2way","limma2wayrepeat","spls2wayrepeat"
        pairedanalysis = FALSE, featselmethod=c("limma"),
        pvalue.thresh=0.05,
        fdrthresh = 0.05, fdrmethod="BH",
        kfold=5,networktype="complete",
        analysismode="classification",pls_vip_thresh = 2,
        num_nodes = 3,
        foldchangethresh=0,
        optselect=TRUE,max_comp_sel=3,saveRda=FALSE,pls.permut.count=NA,
        pca.ellipse=TRUE,
        aggregation.method="none",

        #4) arguments for centrality analysis, WGCNA and global clustering analysis (HCA and EM clustering)
        differential.network.analysis=FALSE, wgcnarsdthresh=1,WGCNAmodules=FALSE,globalclustering=FALSE,

        #5) arguments for correlation and network analysis using the selected features
        cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
    
        #6) arguments for graphical options: see manual for additional arguments
        output.device.type="png",pca.cex.val=4,legendlocation="bottomleft",
        net_node_colors=c("green","red"),
        net_legend=FALSE,aggregation.max.iter=100,
        heatmap.col.opt="redblue",manhattanplot.col.opt=c("darkblue","red3"),
        color.palette=c("journal"),hca_type="two-way",cex.plots=0.6,
        lineplot.lty.option=c("dotted", "solid", "dashed", "dotdash", "longdash", "twodash"),
        timeseries.lineplots=FALSE,lme.modeltype="RI",ylab_text="Intensity",boxplot.type="ggplot",
        multiple.figures.perpanel = FALSE,add.jitter=FALSE,add.pvalues=FALSE,ggplot.type1=TRUE,
        hca.labRow.value = TRUE, hca.labCol.value = FALSE,hca.cex.legend=0.5,
        limma.contrasts.type=c("contr.sum"),
        plot.boxplots.raw=FALSE,vcovHC.type="HC3"
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

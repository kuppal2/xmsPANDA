\name{diffexp}
\alias{diffexp}
\title{diffexp}
\description{
This is a wrapper function for the data analysis workflow. The function performs data preprocessing (filtering, normalization, transformation), PCA, biomarker discovery, correlation-based network analysis, clustering analysis, and functional class scoring.  
The "featselmethod" allows users to select the method for selecting features for classification, regression, ANOVA, and time-series designs. The function evaluates the k-fold cross-validation accuracy using Support Vector Machine, performs hierarchical clustering analysis, PCA analysis (R2/Q2 diagnostics), and generates 
boxplots, bar plots, and line plots for discriminatory features.  
}
\usage{
function(Xmat=NA,Ymat=NA,feature_table_file,parentoutput_dir,
class_labels_file,num_replicates=3,summarize.replicates=TRUE,
summary.method="mean",summary.na.replacement="zeros",missing.val=0,
rep.max.missing.thresh=0.3,
 all.missing.thresh=0.1,group.missing.thresh=0.7,
input.intensity.scale="raw",
log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,
quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,TIC_norm=FALSE,
rangescaling=FALSE,mstus=FALSE,paretoscaling=FALSE,sva_norm=FALSE,
eigenms_norm=FALSE,vsn_norm=FALSE,
normalization.method=c("log2transform"),rsd.filt.list=1,
pairedanalysis=FALSE,featselmethod=c("limma","pls"),fdrthresh=0.05,
fdrmethod="BH",cor.method="spearman",networktype="complete",
abs.cor.thresh=0.4,cor.fdrthresh=0.05,kfold=10,
pred.eval.method="BER",globalcor=TRUE,
target.metab.file=NA,target.mzmatch.diff=10,target.rtmatch.diff=NA,
max.cor.num=100,
numtrees=20000,analysismode="classification",net_node_colors=c("green","red"), net_legend=TRUE,
svm_kernel="radial",heatmap.col.opt="brewer.RdBu",
manhattanplot.col.opt=c("darkblue","red3"),boxplot.col.opt=c("journal"),
barplot.col.opt=c("journal"),sample.col.opt="journal",
lineplot.col.opt="journal",hca_type="two-way",pls_vip_thresh=2,
num_nodes=2,max_varsel=100,
pls_ncomp=5,pca.stage2.eval=FALSE,scoreplot_legend=TRUE,pca.global.eval=TRUE,
rocfeatlist=seq(2,6,1),rocfeatincrement=TRUE,rocclassifier="svm",
foldchangethresh=1,wgcnarsdthresh=20,WGCNAmodules=FALSE,
optselect=TRUE,max_comp_sel=1,saveRda=FALSE,legendlocation="topleft",
pcacenter=TRUE,pcascale=TRUE,pca.cex.val=6,
pca.ellipse=FALSE,ellipse.conf.level=0.95,pls.permut.count=NA,
svm.acc.tolerance=5,limmadecideTests=TRUE,pls.vip.selection="max",
globalclustering=FALSE,plots.res=600,plots.width=8,plots.height=8,
plots.type="cairo",
output.device.type="pdf",pvalue.thresh=0.05,individualsampleplot.col.opt="journal",
pamr.threshold.select.max=FALSE,aggregation.method="RankAggreg",
aggregation.max.iter=1000,mars.gcv.thresh=1,error.bar=TRUE,cex.plots=1,
lme.modeltype="RI",
barplot.xaxis="Factor1",lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),degree_rank_method=NA,match_class_dist=TRUE,timeseries.lineplots=FALSE,
alphabetical.order=FALSE,
kegg_species_code="hsa",database="pathway",reference_set=NA,metab_annot=NA,
add.pvalues=FALSE,add.jitter=FALSE,fcs.permutation.type=1,fcs.method="zscore",
fcs.min.hits=2,names_with_mz_time=NA,ylab_text="Abundance",xlab_text=NA,
boxplot.type="ggplot",samplermindex=NA,
differential.network.analysis = FALSE, degree.centrality.method = "hybrid.DEC", log2.transform.constant = 1,
balance.classes = FALSE, balance.classes.sizefactor = 10, balance.classes.seed = 1, cv.perm.count = NA, 
multiple.figures.perpanel = FALSE, hca.labRow.value = FALSE, hca.labCol.value = FALSE, alpha.col = 1,
similarity.matrix = "correlation", outlier.method = c("pcout", "sumtukey", "pcatukey", "MDchisq"),
removeRda = TRUE, color.palette = c("journal", "npg", "nejm", "jco", "lancet", "custom1", "brewer.RdYlBu", 
"brewer.RdBu", "brewer.PuOr", "brewer.PRGn", "brewer.PiYG", "brewer.BrBG", "brewer.Set2", "brewer.Paired",
"brewer.Dark2", "brewer.YlGnBu", "brewer.YlGn", "brewer.YlOrRd", "brewer.YlOrBr", "brewer.PuBuGn",
"brewer.PuRd", "brewer.PuBu", "brewer.OrRd", "brewer.GnBu", "brewer.BuPu", "brewer.BuGn", "brewer.blues", 
"black", "grey65", "terrain", "rainbow", "heat", "topo"), plot_DiNa_graph = FALSE, 
limma.contrasts.type = c("contr.sum", "contr.treatment"), hca.cex.legend = 0.7, plot.boxplots.raw = FALSE, 
vcovHC.type = "HC3", ggplot.type1 = TRUE, facet.nrow = 1, pairwise.correlation.analysis = FALSE, 
generate.boxplots = TRUE, pvalue.dist.plot = TRUE, ...)
}
\arguments{
 \item{Xmat}{
 R object for feature table. If this is given, then feature table can be set to NA.
 }
  \item{Ymat}{
  R object for response/class labels matrix. If this is given, then class can be set to NA.
  }
  \item{feature_table_file}{
Feature table that includes the mz, retention time, and measured intensity in each sample 
for each analyte. The first 2 columns should be the mz and time. The remaining columns
should correspond to the samples in the class labels file with each column including the intensity profile
of a sample.
Full path required. Eg: C:/My Documents/test.txt
The feature table should be in a tab-delimited format. An example of the input file is provided under the
"example" folder.
}
  \item{parentoutput_dir}{
Provide full path of the folder where you want the results to be written.
Eg: C:/My Documents/ProjectA/results/
}
  \item{class_labels_file}{
File with class labels information for each sample. Samples should be in the same order
as in the feature table. Please use the same format as in the example folder. If you
want to adjust for covariates in "lmreg" option, then you can add additional columns,
one per covariate. Categorical variables should be strings (eg: "male", "female").
Please see "classlabels_gender.txt" file as an example.
}
  \item{num_replicates}{
Number of technical replicates
}
  \item{feat.filt.thresh}{
Percent Intensity Difference or Coefficient of variation threshold; feature filtering
Use NA to skip this step. 
}
  \item{summarize.replicates}{
Do the technical replicates per sample need to be averaged or median summarized? 
}
 \item{summary.method}{
Method for summarizing the replicates. Options: "mean" or "median"
}
 \item{summary.na.replacement}{
How should the missing values be represented? 
Options: "zeros", "halffeaturemin", "halfsamplemin","halfdatamin", "none"
"zeros": replaces missing values by 0
"halfsamplemin": replaces missing value by one-half of the lowest signal intensity in the
corresponding sample
"halfdatamin": replaces missing value by one-half of the lowest signal intensity in the
complete dataset
"halffeaturemin": replaces missing value by one-half of the lowest signal intensity for the
current feature
"none": keeps missing values as NAs

Users are recommended to perform imputation prior to performing biomarker discovery.
}
\item{missing.val}{
How are the missing values represented in the input data? Options: "0" or "NA"
}
 \item{rep.max.missing.thresh}{
What propotion of replicates are allowed to have missing values during the averaging or 
median summarization step of each biological sample? If the number of replicates with
missing values is greater than the defined threshold, then the summarized value is 
represented by the "missing.val" parameter. If the number of replicates with missing values
is less than or equal to the defined threshold, then the summarized value is equal to the 
mean or the median of the non-missing values. Default: 0.5
}
  \item{all.missing.thresh}{
What propotion of total number of samples should have an intensity?
Default: 0.5
}
 \item{input.intensity.scale}{
Are the intensities in the input feature table at raw scale or log2 scale?
eg: "raw" or "log2"
Default: "raw"
}
  \item{group.missing.thresh}{
What propotion of samples in either of the two groups should have an intensity?
If at least x% of the samples in either group have a signal, then the feature is retained
for further analysis. Default: 0.7
}
 \item{normalization.method}{

  Data transformation and normalization method. Options:
  
  "log2transform": log2 transformation
  "log2quantilenorm": log2 transformation and quantile normalization
  "znormtransform": auto scaling; each variable will have a mean of 0 and unit variance
  "quantile_norm": Performs quantile normalization
  "lowess_norm": Performs lowess normalization
  "rangescaling": Performs range scaling; scale by the min and max range
  "paretoscaling": Performs Pareto scaling; scale by the square root of the standard deviation
  "mstus": MS Total Useful Signal (MSTUS) normalization
  "sva": Surrogate Variable Analysis (SVA) normalization
  "eigenms_norm": EigenMS normalization
  "vsn_norm": variance stabilizing normalization
  "tic_norm": totial intensity normalization
  "cubicspline_norm": Cubic spline normalization
  "mad_norm": Median absolute deviation normalization

}
  \item{rsd.filt.list}{
This parameter allows to perform feature filtering based on overall variance (across all samples)
prior to performing hypothesis testing. Eg: seq(0,30,5). 
}
  \item{pairedanalysis}{
 Is this a paired-study design? TRUE or FALSE
 If samples are paired, then the feature table and the class labels file should be organized so that the paired samples 
 are arranged in the same order in each group. For example, the first sample in group A and the first sample in 
 group B should be paired.  
}
\item{featselmethod}{
Options:
"limma": for one-way ANOVA using LIMMA (mode=classification)
"limma2way": for two-way ANOVA using LIMMA (mode=classification)
"limma1wayrepeat": for one-way ANOVA repeated measures using LIMMA (mode=classification)
"limma2wayrepeat": for two-way ANOVA repeated measures using LIMMA (mode=classification)
"lm1wayanova": for one-way ANOVA using linear model (mode=classification)
"lm2wayanova": for two-way ANOVA using linear model (mode=classification)
"lm1wayanovarepeat": for one-way ANOVA repeated measures using linear model (mode=classification)
"lm2wayanovarepeat": for two-way ANOVA repeated measures using linear model (mode=classification)
"lmreg": variable selection based on p-values calculated using a linear regression model; 
allows adjustment for covariates (mode= regression or classification)
"logitreg": variable selection based on p-values calculated using a logistic regression model; 
 allows adjustment for covariates (mode= classification)
"rfesvm": uses recursive feature elimination SVM algorithm for variable selection; 
(mode=classification)
"RF": for random forest based feature selection (mode= regression or classification)
"RFconditional": for conditional random forest based feature selection (mode= regression or classification)
"pamr": for prediction analysis for microarrays algorithm based on the nearest shrunken centroid
method (mode= classification)
"MARS": for multiple adaptive regression splines (MARS) based feature selection
(mode= regression or classification)
"pls": for partial least squares (PLS) based feature selection
(mode= regression or classification)
"spls": for sparse partial least squares (PLS) based feature selection
(mode= regression or classification)
"spls1wayrepeat": for sparse partial least squares (PLS) based feature selection
for one-way repeated measures (mode= regression or classification)
"spls2wayrepeat": for sparse partial least squares (PLS) based feature selection
for two-way repeated measures (mode= regression or classification)
"o1pls": for orthogonal partial least squares (OPLS) based feature selection
(mode= regression or classification)
}
\item{pvalue.thresh}{
p-value threshold. Eg: 0.05^M
}
\item{fdrthresh}{
False discovery rate threshold. Eg: 0.05
}
 \item{fdrmethod}{
Options: "BH", "ST", "Strimmer", "none"
"BH": Benjamini-Hochberg (1995) (Default: more conservative than "ST" and "Strimmer")
"ST": Storey & Tibshirani (Storey 2001, PNAS) algorithm implemented in the qvalue package
"Strimmer": (Strimmer 2008, Bioinformatics) algorithm implemented in the fdrtool package
"none": No FDR correction will be performed. fdrthresh will be treated as raw p-value 
cutoff
}
\item{cor.method}{
Correlation method. Options: "pearson" or "spearman". Default: "spearman"
}
 \item{networktype}{
Options: "complete" or "GGM"
"complete": performs network analysis using ordinary Pearson or Spearman correlation
statistic
"GGM": generates network based on partial correlation analysis using the 
GeneNet package
}
  \item{abs.cor.thresh}{
Absolute Pearson or Spearman correlation coefficient for network analysis. Default: 0.4
}
  \item{cor.fdrthresh}{
False discovery rate threshold for correlation analysis. Default: 0.05
}
\item{kfold}{
	Number of subsets in which the data should be divided 
	for cross-validation. If kfold=10, then the data set will
	be divided into 10 subsets of size (N/10), where N is 
	the total number of samples. 9 subsets are used for
	training and the remaining 1 is used for testing. This
	process is repeated 10 times and the CV-accuracy
	would be the mean of the classification accuracy
	from the 10 iterations. The same will be true for any
	other value of k. 
	WARNING: The kfold value should be
	less than or equal to the total number of samples.
}
 \item{pred.eval.method}{
	Criteria for evaluating the performance of the
	model. CV: Overall Cross-validation classification accuracy, 
	balanced error rate (BER): (sum of accuracy in each class)/(number of classes)
	area under the curve (AUC)
	Eg: "CV", "BER", or "AUC". Default: "BER"
 }
\item{globalcor}{
Perform correlation analysis between selected features and all other features?

Options: TRUE or FALSE
}
  \item{target.metab.file}{
File that includes the mz and/or retention time of the targeted metabolites.
See example.
}
  \item{target.mzmatch.diff}{
+/- ppm mass tolerance for searching the target m/z in the current
	  feature table
}
  \item{target.rtmatch.diff}{
+/- retention time tolerance for searching the target m/z in the current
	  feature table
}
  \item{max.cor.num}{
Maximum number of correlated metabolites to be included in the network
figure. Default: 100
}
\item{pcacenter}{
Data centering for PCA. Options: "TRUE" or "FALSE". Default=TRUE
}
\item{pcascale}{
Data scaling for PCA. Options: "TRUE" or "FALSE". Default=TRUE
}
\item{samplermindex}{
Column index of any additional or irrelevant columns to be deleted.
Options: "NA" or list of column numbers. eg: c(1,3,4) Default=NA
}
\item{numtrees}{
Number of trees to be used for random forest method. Default=500
}
\item{analysismode}{
"classification" for group-wise comparison (case vs control) or
"regression" for continuous response variables. Default: "classification"
}
\item{net_node_colors}{
 Colors of nodes in the correlation networks. Eg: c("pink", "skyblue"), 
or ("red","green") 
}
\item{net_legend}{
Should the network be displayed for the correlation network? eg:
TRUE or FALSE
}
\item{max_var}{
Max number of variables to be used for sPLS, rfesvm, and Random Forest? 
eg:150
}
\item{svm_kernel}{
SVM kernel eg: "radial" or "linear"
}
\item{rf_selmethod}{
Random forest VIP based selection method. 
If rankbased option is selected, variables are ranked 
based on the Variable Importance Measure. 
Only the top "max_varsel" variables are selected. If
absVIMthresh is selected, then all features with VIM greater 
than the absolute value of the lowest VIM are selected.
eg: "absVIMthresh" or "rankbased"
}
\item{pls_vip_thresh}{
Threshold for VIP score from PLS/O1PLS.
eg: 1
}
\item{max_varsel}{
Maximum number of variables to keep if "rankbased"
RF or spls is used. eg: 100
}
\item{pls_ncomp}{
Maximum number of components to be considered 
during the PLS optimal number of components 
selection step.
eg: 2
}
\item{pca.stage2.eval}{
Should PCA diagnostics be performed in stage 2?
eg: TRUE or FALSE
}
\item{scoreplot_legend}{
Should legends be included in score plots?
eg: TRUE or FALSE
}
\item{pca.global.eval}{
Should global PCA evaluation be performed? Default:TRUE
eg: TRUE or FALSE
}
\item{rocfeatlist}{
Vector indicating number of features to be used for ROC
evaluation: eg: c(2,4,6) will generate ROC for top 2, top 4,
and top 6 feautres. Default: seq(2,10,1)
}
\item{rocclassifier}{
Classifier to be used for ROC evaluation. Options: "svm" or
"logitreg". Default: "svm"
}
\item{foldchangethresh}{
Secondary feature selection criteria based on fold change
threshold. This is performed after
statistical significance or importance evaluation.
}
\item{wgcnarsdthresh}{
Relative standard deviation or coefficient of variation
(across all samples) based filtering threshold before performing
WGCNA module preservation analysis. Default: 20
}
\item{WGCNAmodules}{
Perform WGCNA module preservation analysis. TRUE or FALSE
Default: TRUE
}
\item{optselect}{
Determine optimal number of PLS components. Default: TRUE
}
\item{max_comp_sel}{
Number of PLS components to use for VIP or sparse loading
selection (sPLS). Default=1
}
\item{saveRda}{
Should the results be saved in a binary R object. Default: TRUE
}
\item{legendlocation}{
Legend location for PLS or PCA plots
}
 \item{pca.cex.val}{
Size of points on PCA plots. eg: 4
}
\item{pca.ellipse}{
Should ellipse be plotted on PCA plots? 
eg: TRUE or FALSE
}
\item{ellipse.conf.level}{
Confidence interval for PCA ellipses
eg: 0.95
}
\item{pls.permut.count}{
Number of permutations for calculating 
p-values for PLS or sPLS models.
eg: 1000
}
\item{svm.acc.tolerance}{
Stopping criteria for forward feature selection
using "rfeSVM" method. If the difference between
best accuracy and current accuracy based on the newly
added feature drops below the tolerance level, the
forward selection process is terminated.
eg: 5
}
\item{pamr.threshold.select.max}{
If two or more thresholds for shrinking the d statistic
in the PAM algorithm (Tibshirani et al. Statistical
Science 2003) have equal accuracy,
should the maximum value (lowest number of features) be used? 
Default: FALSE
}
\item{aggregation.method}{
Method for combining the results from mutliple feature selection methods
Options: Consensus: will only keep features that are selected in all methods
RankAggreg: will use the cross entropy algorithm with Spearman footrule
distance as the distance measure (RankAggreg; Pihur et al. BMC Bioinformatics 2009)
}
\item{aggregation.max.iter}{
	Maximum number of iterations used in the rank aggregation step.
	Default: 1000
}
\item{mars.gcv.thresh}{
	Minimum generalized cross-validation value (range: 0 to 100) for
	a feature to be selected. Default: 10
}
\item{limmadecideTests}{
	Perform decide tests for LIMMA to perform multiple testing and assign up, 
	down, or not significant. TRUE or FALSE. 
}
\item{pls.vip.selection}{
	How to summarize VIP across multiple PLS components? Options: "max" to take
	the maximum VIP across all selected components or "mean" to take the average
	VIP across all selected components.
	Default: "max"
}
\item{globalclustering}{
	Perform global clustering using all features based on EM and hierachical
	clustering analysis. TRUE or FALSE.
	Default:FALSE
}
\item{plots.res}{
	Resolution of PNG files.
	Default: 600
}
\item{plots.width}{
	Width dimension for PNG files.
	Default: 8
}
\item{plots.height}{
	Height dimension for PNG files.
	Default: 8
}
\item{output.device.type}{
	Output device: "png" or "pdf"
	Default:"pdf"
}
\item{heatmap.col.opt}{
      Color scheme for heatmaps:

      Default: "brewer.RdBu"
}
\item{error.bar}{
	Plot error bars. TRUE or FALSE. Default: TRUE
}
\item{cex.plots}{
	Relative factor to change font size of text in plots
	e.g.: 0.8 or 2
	Default: 1
}
\item{lme.modeltype}{
	Options for mixed-effects models:
	RI:Random intercept
	RIRS: random intercept and random slope models
	Default: "RI"
}
\item{barplot.xaxis}{
	Label for x-axis in barplots
	Default: "Factor1"
}
\item{lineplot.lty.option}{

	Default: c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
}
\item{timeseries.lineplots}{
	Generate lineplots showing longitudinal pattern: TRUE or FALSE
	Default: FALSE
}	
\item{alphabetical.order}{
	Arrange class labels in alphabetical order versus arranging them based on
	which class appears first in the class labels file. TRUE or FALSE
	Default: TRUE
}
\item{boxplot.type}{
	Type of boxplots: "simple" using the boxplot() function in R 
	or "ggplot" for ggboxplot and geom_boxplot functions
}
\item{add.pvalues}{
	Add p-values in boxplots: TRUE or FALSE
	Default: FALSE
}
\item{add.jitter}{
	Add jitter in boxplots: TRUE or FALSE
	Default: FALSE
}
\item{ylab_text}{
	Y-axis label in barplots, boxplots, and lineplots
	Default: "Abundance"
}
  \item{differential.network.analysis}{
  Compare network centrality of each variable between the groups (classes).
  e.g. TRUE or FALSE
}
  \item{degree.centrality.method}{
  Method used for calculating degree centraliy. e.g. "eigenvector", "betweenness",
  "hybrid.DEC","hybrid.DBC"
  
  hybrid.DEC uses local connectivity and eigenvector centrality
  
  hybrid.DBC uses local connectivity and betweenness centrality
  
}
  \item{log2.transform.constant}{
  Small constant to add to values for log2 transformation to avoid
  getting undefined values when calculating log2 of 0
  e.g. 1
}
  \item{balance.classes}{
  Balance classes using the ROSE package for cases when the 
  classes are imbalanced. Only useful for binary or two-class comparisons.
  e.g. FALSE or TRUE
}
  \item{balance.classes.sizefactor}{
  Multiplicative factor (balance.classes.sizefactor x original N) for defining the 
  number of samples in the resulting dataset. e.g. 10
}
  \item{balance.classes.seed}{
  Random number for the ROSE method e.g. 999
}
  \item{cv.perm.count}{
  Number of iterations for calculating the permuted k-fold CV accuracy
  e.g. 100 or 1000
  
  Set to NA to turn off this step
}
  \item{multiple.figures.perpanel}{
  Generate multiple figures per panel in box plots
  e.g. FALSE or TRUE
}
  \item{hca.labRow.value}{
  Show variable (row) names in hierarchical clustering
  analysis heatmaps e.g. TRUE or FALSE
}
  \item{hca.labCol.value}{
  Show sample (column) names in hierarchical clustering
  analysis heatmaps e.g. TRUE or FALSE
}
  \item{alpha.col}{
  Numeric valute to adjust color transparency e.g. 1 or 0.8
}
  \item{similarity.matrix}{
  Similarity matrix to use in hierarchical clustering
  analysis e.g. "correlation" for Pearson or Spearman
  or "TOM" for topological overlap matrix as used in WGCNA
}
  \item{outlier.method}{
  Method for outlier detection for samples
  e.g. "pcout" for the pcout function in the mvoutlier package
  "MDchisq" for using Mahalanobis distance with a Chi-square test
  "pcatukey" for using PCA scores with Tukey's +/- 1.5 interquartile
  range (IQR) rule
  "sumtukey" for applying Tukey's +/- 1.5 IQR rule on summed abundance
  or intensities across all variables for each sample
}
  \item{removeRda}{
  Remove intermediate .Rda files 
  e.g. TRUE or FALSE
}
  \item{color.palette}{
  Color theme for plots. default: "journal"
  Options:
        1. "journal": color-blind friendly palette
        2. built-in R color palettes: "rainbow", "terrain", "heat","topo"
        3. RColorBrewer pallettes: "brewer.YlOrRd", "brewer.Purples","brewer.YlGn",
        "brewer.BuPu","brewer.BuGn","brewer.GnBu", "brewer.YlGnBu", "brewer.RdBu",
        "brewer.RdYlBu","brewer.PuOr","brewer.PRGn" (color codes: Yl-yellow; Rd-red,
        Bu-blue, Or-orange, Gn-green, PR-purple)
        4. Generate a custom palette by providing colors (e.g. c("orange","blue","green"))
        
  Please the color_palettes_xmsPANDA.pdf
  file on the GitHub page under xmsPANDA/inst
}
  \item{plot_DiNa_graph}{
  Show network plots generated using differential network 
  analysis e.g. TRUE or FALSE
}
  \item{limma.contrasts.type}{
  Contrasts method for LIMMA
  e.g. "contr.sum" for ANOVA like sum contrasts method
  "contr.treatment" to treat the first group as the reference group
  and all other groups are compared to the reference
}
  \item{hca.cex.legend}{
  Numeric value indicating the amount by which plotting text and
  symbols should be scaled relative to the default. e.g 0.7
  
  Set to NA to hide the HCA legend
}
  \item{plot.boxplots.raw}{
  Should the box plots using untransformed (raw)
  data be generated? e.g. TRUE or FALSE
}
  \item{vcovHC.type}{
  Heteroscedasticity-consistent covariance matrix estimation type
  Please see the type argument in the vcovHC
  function in R package sandwich for more details
  e.g. "HC3", "HC0"
}
  \item{ggplot.type1}{
  Should the boxplots be sub-grouped (faceted) by Factor2?
  Applicable for two-factor study designs. e.g. TRUE or FALSE
  
}
  \item{facet.nrow}{
  Number of rows in box plots if using ggplots for two-factor
  study designs. e.g. 1
}
  \item{pairwise.correlation.analysis}{
  Should the pairwise correlation analysis between selected variables
  be performed? e.g. TRUE or FALSE
  
  Note that this is only for the variables meeting the feature selection
  criteria.
}
  \item{generate.boxplots}{
  Should the boxplots be generated? e.g. TRUE or FALSE
}
  \item{pvalue.dist.plot}{
  Should the p-value distribution plots be generated? e.g. TRUE or FALSE
}
}
\details{
This function performs data transformation, normalization, FDR analysis using LIMMA, 
variable selection using random forests, evaluates the predictive accuracy of
the FDR significant features using k-fold cross-validation with a Support Vector
Machine classifier, performs two-way hierarchical clustering analysis, and principal
component analysis. An optimizaiton scheme is used to select the best set of features
from different log2 fold change filtering thresholds. Finally, metabolome-wide and
targeted correlation based network analysis of the FDR significant features is performed.
}
\value{
  \item{diffexp_metabs}{Best set of discriminatory features.}
  \item{all_metabs}{Results for all features.}
  \item{mw.an.fdr}{Metabolome-wide significant correlation network of 
  differentially expressed metabolites.}
  \item{targeted.an.fdr}{Correlation network of  differentially expressed
  metabolites with targeted metabolites.}
  
  Following files are generated in the parent output location:
Manhattan plots: showing metabolome wide p-values; 
Heatmap from Two-way hierarchical clustering analysis;
Pairwise score plots from Principal Component Analysis;
PCA score distribution plots;
ROC plots;
List of differentially expressed metabolites;
Boxplots of differentially expressed metabolites;
Correlation network figure and matrix;
Pairwise correlation matrix CIRCOS format ready to be uploaded to:
  http://mkweb.bcgsc.ca/tableviewer/visualize/
Or uploaded to Cytoscape gml format

}
\author{Karan Uppal <kuppal3gt@gmail.com>}

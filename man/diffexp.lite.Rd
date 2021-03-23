\name{diffexp.lite}
\alias{diffexp.lite}
\title{diffexp.lite}
\description{
A "lite" version of the diffexp() function with key features for data preprocessing, feature selection, PCA/HCA, and correlation analysis
}
\usage{
diffexp.lite(Xmat = NA, Ymat = NA, outloc = NA, summary.na.replacement = "halffeaturemin", missing.val = 0, all.missing.thresh = 0.1, group.missing.thresh = 0.7, input.intensity.scale = "raw", normalization.method = c("log2transform", "znormtransform", "lowess_norm", "log2quantilenorm", "quantile_norm", "rangescaling", "paretoscaling", "mstus", "eigenms_norm", "vsn_norm", "sva_norm", "tic_norm", "cubicspline_norm", "mad_norm", "none"), rsd.filt.list = 1, pca.global.eval = TRUE, pairedanalysis = FALSE, featselmethod = c("limma"), pvalue.thresh = 0.05, fdrthresh = 0.05, fdrmethod = "BH", foldchangethresh = 0, analysismode = "classification", pls_vip_thresh = 2, num_nodes = 2, optselect = TRUE, max_comp_sel = 1, hca_type = "two-way", output.device.type = "png", timeseries.lineplots = FALSE, alphabetical.order = FALSE, ylab_text = "Abundance", boxplot.type = "ggplot", color.palette = c("journal", "npg", "nejm", "jco", "lancet", "custom1", "brewer.RdYlBu", "brewer.RdBu", "brewer.PuOr", "brewer.PRGn", "brewer.PiYG", "brewer.BrBG", "brewer.Set2", "brewer.Paired", "brewer.Dark2", "brewer.YlGnBu", "brewer.YlGn", "brewer.YlOrRd", "brewer.YlOrBr", "brewer.PuBuGn", "brewer.PuRd", "brewer.PuBu", "brewer.OrRd", "brewer.GnBu", "brewer.BuPu", "brewer.BuGn", "brewer.blues", "black", "grey65", "terrain", "rainbow", "heat", "topo"), generate.boxplots = TRUE, hca.cex.legend = 0.7, lme.modeltype = "lme.RI", globalcor = FALSE, abs.cor.thresh = 0.4, cor.fdrthresh = NA, net_legend = TRUE, rocclassifier = "svm", limma.contrasts.type = c("contr.sum", "contr.treatment"), limmadecideTests = FALSE, cex.plots = 0.9,hca.labRow.value = TRUE,hca.labCol.value = TRUE,...)
}
\arguments{
  \item{Xmat}{
 R object for feature table. If this is given, then feature table can be set to NA.
 }
  \item{Ymat}{
  R object for response/class labels matrix. If this is given, then class can be set to NA.
  }
  \item{outloc}{
Provide full path of the folder where you want the results to be written.
Eg: C:/My Documents/ProjectA/results/
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
  \item{all.missing.thresh}{
What propotion of total number of samples should have an intensity?
Default: 0.5
}
  \item{group.missing.thresh}{
What propotion of samples in either of the two groups should have an intensity?
If at least x% of the samples in either group have a signal, then the feature is retained
for further analysis. Default: 0.7
}
 \item{input.intensity.scale}{
Are the intensities in the input feature table at raw scale or log2 scale?
eg: "raw" or "log2"
Default: "raw"
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
 \item{pairedanalysis}{
 Is this a paired-study design? TRUE or FALSE
 If samples are paired, then the feature table and the class labels file should be organized so that the paired samples 
 are arranged in the same order in each group. For example, the first sample in group A and the first sample in 
 group B should be paired.  
}
 \item{pca.global.eval}{
  Perform PCA using all variables? TRUE or FALSE
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
Options: "BH", "ST", "Strimmer", "BY","none"
"BH": Benjamini-Hochberg (1995) (Default: more conservative than "ST" and "Strimmer")
"ST": Storey & Tibshirani (Storey 2001, PNAS) algorithm implemented in the qvalue package
"Strimmer": (Strimmer 2008, Bioinformatics) algorithm implemented in the fdrtool package
"none": No FDR correction will be performed. fdrthresh will be treated as raw p-value 
cutoff
}
\item{foldchangethresh}{
Secondary feature selection criteria based on fold change
threshold. This is performed after
statistical significance or importance evaluation.
}
\item{analysismode}{
"classification" for group-wise comparison (case vs control) or
"regression" for continuous response variables. Default: "classification"
}
\item{pls_vip_thresh}{
Threshold for VIP score from PLS/O1PLS.
eg: 1
}
  \item{num_nodes}{
 Number of CPU cores to use e.g. 2
}
\item{optselect}{
Determine optimal number of PLS components. Default: TRUE
}
\item{max_comp_sel}{
Number of PLS components to use for VIP or sparse loading
selection (sPLS). Default=1
}
  \item{hca_type}{
  "one-way" or "two-way" HCA
}
  \item{output.device.type}{
  "pdf" or "png"
}
 \item{timeseries.lineplots}{
	Generate lineplots showing longitudinal pattern: TRUE or FALSE
	Default: FALSE
}	
 \item{alphabetical.order}{
	Arrange class labels in alphabetical order versus arranging them based on
	which class appears first in the class labels file. TRUE or FALSE

}
 \item{ylab_text}{
	Y-axis label in barplots, boxplots, and lineplots
	Default: "Abundance"
}
  \item{boxplot.type}{
	Type of boxplots: "simple" using the boxplot() function in R 
	or "ggplot" for ggboxplot and geom_boxplot functions
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
  \item{generate.boxplots}{
  Should the boxplots be generated? e.g. TRUE or FALSE
}
  \item{hca.cex.legend}{
  Numeric value indicating the amount by which plotting text and
  symbols should be scaled relative to the default. e.g 0.7
  
  Set to NA to hide the HCA legend
}
  \item{lme.modeltype}{
	Options for mixed-effects models:
	RI:Random intercept
	RIRS: random intercept and random slope models
	Default: "RI"
}
 \item{globalcor}{
Perform correlation analysis between selected features and all other features?

Options: TRUE or FALSE
}
 \item{abs.cor.thresh}{
Absolute Pearson correlation coefficient for network analysis. Default: 0.4
}
 \item{cor.fdrthresh}{
False discovery rate threshold for correlation analysis. Default: 0.05
}
\item{net_legend}{
Should the network be displayed for the correlation network? eg:
TRUE or FALSE
}
  \item{rocclassifier}{
 Set to NA to turn off k-fold CV classification accuracy and ROC analysis
 Default: "svm"
}
 \item{limma.contrasts.type}{
  Contrasts method for LIMMA
  e.g. "contr.sum" for ANOVA like sum contrasts method
  "contr.treatment" to treat the first group as the reference group
  and all other groups are compared to the reference
}
 \item{limmadecideTests}{
	Perform decide tests for LIMMA to perform multiple testing and assign up, 
	down, or not significant. TRUE or FALSE. 
}
\item{cex.plots}{
	Relative factor to change font size of text in plots
	e.g.: 0.8 or 2
	Default: 1
}
\item{hca.labRow.value}{
  Show variable (row) names in hierarchical clustering
  analysis heatmaps e.g. TRUE or FALSE
}
  \item{hca.labCol.value}{
  Show sample (column) names in hierarchical clustering
  analysis heatmaps e.g. TRUE or FALSE
}
  
}
\details{
The "lite" version requires fewer computational resources. The function performs 
data transformation, normalization, feature selection, evaluates the predictive accuracy of
the FDR significant features using k-fold cross-validation with a Support Vector
Machine classifier, performs hierarchical clustering analysis, correlation analysis,
and principal component analysis. 
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

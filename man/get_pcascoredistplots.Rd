\name{get_pcascoredistplots}
\alias{get_pcascoredistplots}
\title{get_pcascoredistplots}
\description{PCA score distribution (25th percentile, median, 75th percentile) plots}
\usage{
get_pcascoredistplots(X = NA, Y = NA, feature_table_file, parentoutput_dir,
class_labels_file, sample.col.opt = "rainbow", plots.width = 2000,
plots.height = 2000, plots.res = 300, alphacol = 0.3, col_vec,
pairedanalysis = FALSE, pca.cex.val = 3, legendlocation = "topright",
pca.ellipse = TRUE, ellipse.conf.level = 0.5, filename = "all")
}
\arguments{
  \item{X}{
 R object for feature table. If this is given, then feature table can be set to NA.
}
  \item{Y}{
 R object for response/class labels matrix. If this is given, then class labels file can be set to NA.
}
  \item{feature_table_file}{
Path and name of feature table that includes the mz, retention time, and measured intensity
in each sample for each analyte. The first 2 columns should be the mz and time. 
The remaining columns should correspond to the samples in the class labels file 
with each column including the intensity profile of a sample.
Full path required. Eg: C:/My Documents/test.txt
The feature table should be in a tab-delimited format. An example of the input file is
provided under the "example" folder.
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
  \item{sample.col.opt}{
Color scheme for PCA and heatmap sample axis 
eg: "rainbow", "heat" or "topo"
}
  \item{alphacol}{
Color scaling parameter
eg:0.3
}
  \item{col_vec}{
Vector of colors for each sample.
}
  \item{pairedanalysis}{
 Is this a paired-study design? TRUE or FALSE
 If samples are paired, then the feature table and the class labels file should be organized so that the paired samples 
 are arranged in the same order in each group. For example, the first sample in group A and the first sample in 
 group B should be paired.  
}
  \item{pca.cex.val}{
Size of points on PCA plots. eg: 4
}
  \item{legendlocation}{
Legend location on PCA plots
}
  \item{pca.ellipse}{
Should ellipse be plotted on PCA plots? 
eg: TRUE or FALSE
}
  \item{ellipse.conf.level}{
Confidence interval for PCA ellipses
eg: 0.95
}
  \item{filename}{
Name of output PDF file
}
}
\details{
This function performs PCA and generates pariwise score plots
as well as score distribution plots (per group). It uses the Y vector
and classlabels for color coding the samples in the pairwise
score plots.
}
\value{
The output includes:
Pairwise PCA score plots,
PCA score distribution plots, 
PCA scores and loadings text files.
}
\author{
Karan Uppal
}
\note{
The plots can be sent to an external device
by running the following commands:
pdf("get_pcascoredistplots.pdf")
get_pcascoredistplots(...)
dev.off()
}
\keyword{ PCA }

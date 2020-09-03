\name{metabnet}
\alias{metabnet}
\title{metabnet}
\description{
Function for correlation (complete or partial) based metabolome-wide network analysis. 
Additionally, users have the option to provide a matrix of m/z features corresponding to 
chemicals of interest such as (phenylalanine, choline, etc) and/or a matrix of m/z
features corresponding to discriminatory metabolites.
}
\usage{
metabnet(feature_table_file, target.metab.file, sig.metab.file,
 class_labels_file = NA, parentoutput_dir, num_replicates = 1,
 cor.method = "spearman", abs.cor.thresh = 0.4, cor.fdrthresh = 0.05,
 target.mzmatch.diff = 10, target.rtmatch.diff = NA, max.cor.num = 100, 
feat.filt.thresh = NA, summarize.replicates = TRUE, summary.method = "mean", all.missing.thresh = 0.5, group.missing.thresh = 0.7, log2transform = TRUE,
 medcenter = TRUE, znormtransform = FALSE, quantile_norm = TRUE,
 lowess_norm = FALSE, madscaling = FALSE, missing.val = 0, 
networktype = "complete", samplermindex = NA, rep.max.missing.thresh = 0.3, summary.na.replacement = "zeros", net_node_colors = c("pink", "skyblue"), 
net_legend = FALSE, netrandseed = 555, TIC_norm = FALSE, 
normalization.method = c("none"), input.intensity.scale = "raw")
}
\arguments{
 \item{feature_table_file}{
Feature table that includes the mz, retention time, and measured intensity in each sample 
for each analyte. The first 2 columns should be the mz and time. The remaining columns
should correspond to the samples in the class labels file with each column including the intensity profile
of a sample.
Full path required. Eg: C:/My Documents/test.txt
The feature table should be in a tab-delimited format. An example of the input file is provided under the
"example" folder.
}
 \item{target.metab.file}{
File that includes the mz and/or retention time of the targeted metabolites corresponding
to pathways or chemicals of interest.
See example.
}
  \item{sig.metab.file}{
File that includes the mz and/or retention time of the discriminatory metabolites.
See example.
}
 \item{class_labels_file}{
File with class labels information for each sample. Samples should be in the same order
as in the feature table. Please use the same format as in the example folder.
}
  \item{parentoutput_dir}{
Provide full path of the folder where you want the results to be written.
Eg: C:/My Documents/ProjectA/results/
}
  \item{num_replicates}{
Number of technical replicates
}
\item{cor.method}{
Correlation method. Options: "pearson" or "spearman". Default: "spearman"
}
  \item{abs.cor.thresh}{
Absolute Pearson correlation coefficient for network analysis. Eg: 0.5
}
  \item{cor.fdrthresh}{
False discovery rate threshold for correlation analysis. Eg: 0.05
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
Options: "zeros", "halfsamplemin","halfdatamin", "none"
"zeros": replaces missing values by 0
"halfsamplemin": replaces missing value by one-half of the lowest signal intensity in the
corresponding sample
"halfdatamin": replaces missing value by one-half of the lowest signal intensity in the
complete dataset
"none": keeps missing values as NAs

Users are recommended to perform imputation prior to performing biomarker discovery.
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
  \item{log2transform}{
Data transformation: Please refer to http://www.biomedcentral.com/1471-2164/7/142
Try different combinations; such as log2transform=TRUE, znormtransfrom=FALSE
or log2transform=FALSE, znormtransfrom=TRUE
}
  \item{medcenter}{
Median centering of metabolites
}
  \item{znormtransform}{
Auto scaling; each metabolite will have a mean of 0 and unit variance
}
  \item{quantile_norm}{
Performs quantile normalization. Normalization options: Please set only one of the options to be TRUE
}
  \item{lowess_norm}{
Performs lowess normalization. Normalization options: Please set only one of the options to be TRUE
}
  \item{madscaling}{
Performs median adjusted scale normalization. Normalization options: Please set only one of the options to be TRUE
}
\item{missing.val}{
How are the missing values represented in the input data? Options: "0" or "NA"
}
 \item{networktype}{
Options: "complete" or "GGM"
"complete": performs network analysis using ordinary Pearson or Spearman correlation
statistic
"GGM": generates network based on partial correlation analysis using the 
GeneNet package
}
\item{samplermindex}{
Column index of any additional or irrelevant columns to be deleted.
Options: "NA" or list of column numbers. eg: c(1,3,4) Default=NA
}
 \item{rep.max.missing.thresh}{
What propotion of replicates are allowed to have missing values during the averaging or 
median summarization step of each biological sample? If the number of replicates with
missing values is greater than the defined threshold, then the summarized value is 
represented by the "missing.val" parameter. If the number of replicates with missing values
is less than or equal to the defined threshold, then the summarized value is equal to the 
mean or the median of the non-missing values. Default: 0.5
}
\item{net_node_colors}{
 Colors of nodes in the correlation networks. Eg: c("pink", "skyblue"), 
or ("red","green") 
}
\item{net_legend}{
Should the network be displayed for the correlation network? eg:
TRUE or FALSE
}
}
\details{
Function for metabolomic network analysis
}
\value{
Correlation matrix and network of metabolites.
}
\references{

}
\author{Karan Uppal <kuppal2@emory.edu>}

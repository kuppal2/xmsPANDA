\name{data_preprocess}
\alias{data_preprocess}
\title{data_preprocess}
\description{
This function performs data transformation, normalization
}
\usage{
function(Xmat = NA, Ymat = NA, feature_table_file = NA,
                 parentoutput_dir, class_labels_file = NA,
                 num_replicates = 3, feat.filt.thresh = NA,
                 summarize.replicates = TRUE, summary.method = "mean",
                 all.missing.thresh = 0.5, group.missing.thresh = 0.7,
                 log2transform = TRUE, medcenter = FALSE,
                 znormtransform = FALSE, quantile_norm = TRUE,
                 lowess_norm = FALSE, rangescaling = FALSE,
                 paretoscaling = FALSE, mstus = FALSE, sva_norm =
                 FALSE, TIC_norm = FALSE, eigenms_norm = FALSE,
                 madscaling = FALSE, vsn_norm = FALSE, cubicspline_norm
                 = FALSE, missing.val = 0, samplermindex = NA,
                 rep.max.missing.thresh = 0.5, summary.na.replacement =
                 "zeros", featselmethod = NA, pairedanalysis = FALSE,
                 normalization.method = "none", input.intensity.scale =
                 "raw", create.new.folder = TRUE)
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
as in the feature table. Please use the same format as in the example folder.
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
}
\value{
Pre-processed data matrix.
}
\author{Karan Uppal <kuppal3gt@gmail.com>}

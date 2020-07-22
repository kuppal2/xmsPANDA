\name{do_wgcna}
\alias{do_wgcna}
\title{do_wgcna}
\description{This function performs module preservation analysis 
using WGCNA.}
\usage{
do_wgcna(feature_table_file = NA, class_labels_file = NA, X = NA, 
Y = NA, sigfeats = NA)
}
\arguments{
  \item{feature_table_file}{
Path and name of feature table that includes the mz, retention time, and measured intensity
in each sample for each analyte. The first 2 columns should be the mz and time. 
The remaining columns should correspond to the samples in the class labels file 
with each column including the intensity profile of a sample.
Full path required. Eg: C:/My Documents/test.txt
The feature table should be in a tab-delimited format. An example of the input file is
provided under the "example" folder.
}
  \item{class_labels_file}{
File with class labels information for each sample. Samples should be in the same order
as in the feature table. Please use the same format as in the example folder. If you
want to adjust for covariates in "lmreg" option, then you can add additional columns,
one per covariate. Categorical variables should be strings (eg: "male", "female").
Please see "classlabels_gender.txt" file as an example.
}
  \item{X}{
R object for feature table. If this is given, then feature table can be set to NA.
}
  \item{Y}{
R object for response/class labels matrix. If this is given, then class labels file
can be set to NA.
}
  \item{sigfeats}{
List of differentially expressed features. Default: NA
}
}
\details{
This function calls WGCNA to perform module preservation
analysis between different classes or groups.
}
\value{
PDF plots for module preservation from WGCNA and
preservation matrix
}
\references{
WGCNA (Horvath 2007)
}
\author{
Karan Uppal
}
\keyword{ ~WGCNA }

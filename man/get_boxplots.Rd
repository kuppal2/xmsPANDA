\name{get_boxplots}
\alias{get_boxplots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get_boxplots}
\description{
This function generates boxplots for m/z features. The input intensity matrix could be transformed
or non-transformed intensities. Sample labels in the class labels file should be
in the same order as the intensity matrix or feature table.
}
\usage{
get_boxplots(feature_table_file, parentoutput_dir, class_labels_file, sample.col.opt = "rainbow", plots.width = 2000, plots.height = 2000, plots.res = 300, alphacol = 0.3)
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
  \item{parentoutput_dir}{
Provide full path of the folder where you want the results to be written.
Eg: C:/My Documents/ProjectA/results/
}
  \item{class_labels_file}{
File with class labels information for each sample. Samples should be in the same order
as in the feature table. Please use the same format as in the example folder.
}
\item{sample.col.opt}{
Color scheme for PCA and heatmap sample axis 
eg: "heat" or "topo"
}
\item{alphacol=0.3}{
Color scaling parameter
eg:0.3
}
}
\value{
Creates a PDF with boxplots for each m/z feature.
}
\author{Karan Uppal <kuppal2@emory.edu>}
\name{get_hca}
\alias{get_hca}
\title{get_hca}
\description{
This function performs two-way hierarchical clustering analysis and generates a heatmap showing
the clustering results. The input intensity matrix could be transformed
or non-transformed intensities. Sample labels in the class labels file should be
in the same order as the intensity matrix or feature table.
}
\usage{
get_hca(feature_table_file, parentoutput_dir, class_labels_file, heatmap.col.opt = "RdBu", cor.method = "spearman",
is.data.znorm = FALSE, analysismode = "classification", sample.col.opt = "rainbow", plots.width = 2000,
plots.height = 2000, plots.res = 300, alphacol = 0.3,hca_type="two-way")
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
\item{heatmap.col.opt}{
Color scheme for HCA hetmap eg: "RdBu", "topo", "heat", or "terrain"
}
\item{cor.method}{
Correlation method. Options: "person" or "spearman". Default: "spearman"
}
  \item{is.data.znorm}{
Is the data already auto-scaled or z-scaled? eg: TRUE or FALSE
}
\item{analysismode}{
"classification" for group-wise comparison (case vs control) or
"regression" for continuous response variables. Default: "classification"
}
\item{sample.col.opt}{
Color scheme for PCA and heatmap sample axis 
eg: "heat" or "topo" or "rainbow"
}
  \item{plots.width}{
Width of the tiff file. eg: 2000
}
  \item{plots.height}{
Height of the tiff file. eg: 2000
}
  \item{plots.res}{
Resolution of the tiff file. eg: 300
}
\item{alphacol}{
Color scaling parameter
eg:0.3
}
\item{hcatype}{
Color scaling parameter
eg:"two-way" or "one-way"
}
}
\value{
Heatmap from Two-way hierarchical clustering analysis;
Intensity matrix in the same order as the dendrograms in
heatmap; 
Sample cluster labels
}
\author{Karan Uppal <kuppal2@emory.edu>}
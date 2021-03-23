\name{get_pca}
\alias{get_pca}
\title{
Perfors PCA analysis
}
\description{
This function uses the pca function implemented in 
the mixOmics package for PCA analysis
}
\usage{
get_pca(X, samplelabels, legendlocation = "topright", filename = NA, 
ncomp = 5, center = TRUE, scale = TRUE, legendcex = 0.5, 
outloc = getwd(), col_vec = NA, sample.col.opt = "default", 
alphacol = 0.3, class_levels = NA, pca.cex.val = 3,
pca.ellipse = TRUE, ellipse.conf.level = 0.5, samplenames = FALSE)
}
\arguments{
  \item{X}{
Data matrix without m/z and time.
}
  \item{samplelabels}{
Vector with class label for each sample.
}
  \item{legendlocation}{
Location of the legend on PCA score plots
}
  \item{filename}{
eg: "all", "signficantfeats"
}
  \item{ncomp}{
Number of components; please use
?pca for more information
}
  \item{center}{
Should the data be centered?; please use
?pca for more information
}
  \item{scale}{
Should the data be scaled?; please use
?pca for more information
}
  \item{legendcex}{
Size of the legend text in the PCA score plots.
e.g.: 0.5 or 0.7
}
  \item{outloc}{
Output folder location
}
  \item{col_vec}{
Provide vector of colors for each
group. eg: NA or  c("red","green")
for cases and controls, respectively.
This argument is ignored if sample.col.opt
is provided
}
  \item{sample.col.opt}{
Select R color palette.
eg: "rainbow", "terrain", "topo".
"heat", "default"
}
  \item{alphacol}{
Semi-transparent colors
eg: 0.2
}
  \item{class_levels}{
Vector with names of different
sample groups. 
eg: c("case", "control")
or NA
}
  \item{pca.cex.val}{
Size of dots in PCA score
plots. eg: 0.4
}
  \item{pca.ellipse}{
Should the score confidence
interval for each group be drawn?
eg: TRUE or FALSE
}
  \item{ellipse.conf.level}{
Confidence interval level
eg: 0.95
}
  \item{samplenames}{
Should the sample names be included in
PCA plots? eg: TRUE or FALSE
}
}
\details{
This function performs PCA analysis. 
The results are saved in a RDA file.

}
\value{
The function returns PCA results
as an object and generates pairwise
score plots for the first three
components
}
\references{
mixOmics
}
\author{
Karan Uppal
}
\keyword{ ~PCA }

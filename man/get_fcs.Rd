\name{get_fcs}
\alias{get_fcs}
\title{get_fcs}
\description{
Functional class scoring
}
\usage{
get_fcs(target.data, target.data.annot = NA, kegg_species_code = "hsa", 
database = c("pathway","module","brite","lipidmaps_mainclass","lipidmaps_subclass",
"refmet_superclass","refmet_mainclass","refmet_subclass","reactome_compound",
"reactome_atlas","kegg_atlas","custom"),
reference_set = NA, type.statistic = "pvalue", fcs.min.hits = 2, itrs = 100, numnodes = 2)
}
\arguments{
  \item{target.data}{
  Data frame with variable ID in column A and the statistic (eg. pvalue, 
  t-statistic, fold change, VIP) in column B
}
  \item{target.data.annot}{
  Optional argument with annotation data for the variables in target.data
}
  \item{kegg_species_code}{
  KEGG species code if using KEGG as the reference database (e.g. hsa)
}
  \item{database}{
  Options include: c("pathway","module","brite","lipidmaps_mainclass",
  "lipidmaps_subclass", "refmet_superclass","refmet_mainclass",
  "refmet_subclass","reactome_compound",
"reactome_atlas","kegg_atlas","custom")

  The reactome_atlas and kegg_atlas options include all genes, compounds,
  and proteins associated with each pathway in the databases.
}
  \item{reference_set}{
  If cutom database option is selected, then a data frame with the 
  reference database should be provided with the functional class/pathway 
  ID in column A, variable ID (e.g. compound name/ID) in column B, 
  and set/pathway name in column C
}
  \item{type.statistic}{
  Type of test statistic (e.g. pvalue, t-statistic,fold change)
}
  \item{fcs.min.hits}{
  Minimum number of hits in a functional class (e.g. 2)
}
  \item{itrs}{
  Number of permutations to generate the null distribution
  (e.g. 100)
}
  \item{numnodes}{
  Number of CPUs to use (e.g. 2)
}
}
\details{
The algorithm uses the z-statistic method proposed by Irizarry et al.(2009)
and the max-mean method proposed by Efron and Tibshirani (2006) for determining
significance of each functional class. The p-values from the two methods are 
aggregated using the Chi-square method.
}
\value{
A data frame with aggregated statistic, p-value 
}
\author{
Karan Uppal
}

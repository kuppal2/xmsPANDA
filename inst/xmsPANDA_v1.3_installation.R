if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
BiocManager::install(version = "3.11")

BiocManager::install(c("mixOmics","CMA","GO.db","impute","limma","qvalue","pcaMethods","KEGGREST","genefilter"),suppressUpdates=TRUE,dependencies=TRUE)

install.packages(c('rgl', 'snow', 'doSNOW', 'foreach', 'e1071', 'WGCNA','reshape2','robust','mvoutlier',
'randomForest', 'party', 'fdrtool', 'GeneNet', 'corpcor', 'earth',
'pROC', 'multcomp', 'pls', 'plsgenomics', 'igraph', 'ROCR',
'flashClust', 'data.table', 'dplyr', 'mclust', 'RankAggreg', 'pamr', "tidyverse","grid",
'sandwich', 'Boruta', 'lsmeans', 'car', 'ggpubr', 'extrafont', 'stepPlr', 'h2o','shinyBS', 'V8', 'shinyWidgets', 
                   'plotly','shinycssloaders', 'colourpicker',"raster","ROSE","devtools"),dependencies=TRUE,
                 repos="https://cran.r-project.org")


Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
library(devtools); install_github("kuppal2/xmsPANDA")
#launch xmsPANDA
library(xmsPANDA)
runApp.xmsPANDA()

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install(c("mixOmics","CMA","GO.db","impute","limma","qvalue","pcaMethods","KEGGREST","genefilter"),suppressUpdates=TRUE,dependencies=TRUE)

install.packages(c('rgl', 'snow', 'doSNOW', 'foreach', 'e1071', 'WGCNA','reshape2',
'randomForest', 'party', 'fdrtool', 'GeneNet', 'corpcor', 'earth',
'pROC', 'multcomp', 'pls', 'plsgenomics', 'igraph', 'ROCR',
'flashClust', 'data.table', 'dplyr', 'mclust', 'RankAggreg', 'pamr',
'sandwich', 'Boruta', 'lsmeans', 'car', 'ggpubr', 'extrafont', 'stepPlr', 'h2o','shinyBS', 'V8', 'shinyWidgets', 
                   'plotly','shinycssloaders', 'colourpicker',"raster","ROSE","devtools"),dependencies=TRUE,
                 repos="https://cran.r-project.org")

library(devtools); install_github("kuppal2/xmsPANDA")

#launch xmsPANDA
library(xmsPANDA)
runApp.xmsPANDA()

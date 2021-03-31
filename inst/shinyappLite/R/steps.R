
data_preprocessing<-fluidRow(
  tags$div(
    id="maindiv1",
    column(width=12,
    mainPanel(
      width=12,
      bsCollapse(open = "Replicate summarization",

                 bsCollapsePanel("Filtering",
                                 column(width=6,numericInput(width="350px","all_missing_thresh", "Minimum non-missing sample ratio:", 0.5, min = 0, max = 1)),
                                 bsTooltip("all_missing_thresh", "What propotion of total number of samples should have an intensity?","bottom"),
                                 column(width=6,numericInput(width="350px","rsd_filt_list", "Minimum overall variance:", 0, min = 0, max = 100)),
                                 bsTooltip("rsd_filt_list", "Minimum relative standard deviation across all samples","bottom"),
                                 column(width=6,numericInput(width="350px","group_missing_thresh", "Minimum non-missing sample ratio for group:", 0.8, min = 0, max = 1)),
                                 bsTooltip("group_missing_thresh", "Minimum propotion of samples in at least one group in which a non-missing signal value should be present","bottom"),
                                 style = "primary"),
                 
                bsCollapsePanel("Imputation, transformation, and normalization",
                 column(width=12,selectInput(width="350px","summary_na_replacement","Choose an imputation method:",c("halffeaturemin","knn","randomforest","zeros","halfsamplemin","halfdatamin","none"))),
                 bsTooltip("summary_na_replacement", "How are the missing values represented?","bottom"),
                column(width=6,selectInput(width="350px","normmethod","Normalization method:",c("No normalization"="none","log2 and quantile normalization"="log2quantilenorm",
                "log2 transformation"="log2transform",
                "auto-scaling"="znormtransform",
                "lowess normalization"="lowess_norm",
                "quantile normalization"="quantile_norm",
                "range scaling"="rangescaling",
                "paretoscaling"="paretoscaling",
                "most useful total signal"="mstus",
                "EigenMS normalizaiton"="eigenms_norm",
                 "Variance stabilizing normalization"="vsn_norm",
                 "Surrogate variable analysis"="sva_norm",
                  "Total ion intensity normalization"="tic_norm",
                  "cubicspline normalization"="cubicspline_norm",
                  "median absolute deviation normalization"="mad_norm"),selected=c("log2transform"))),
#column(width=6,selectInput(width="350px","log2transform","Perform log2 transformation?",c("TRUE","FALSE"))),
#                                column(width=6,selectInput(width="350px","znormtransform","Perform auto-scaling transformation?",c("FALSE","TRUE"))),
#                                column(width=6,selectInput(width="350px","quantile_norm","Perform quantile normalization?",c("TRUE","FALSE"))),
#                                column(width=6,selectInput(width="350px","TIC_norm","Perform total ion count normalization?",c("FALSE","TRUE"))),
                                 style = "primary")
      ))))
)



network_analysis<-fluidRow(
  tags$div(
    id="maindiv3",
    column(12,radioButtons("globalcor", "3.1 Perform correlation-based network analysis of selected features:", inline=TRUE,c(True = "TRUE",False = "FALSE"),selected = "FALSE")),
    conditionalPanel(
      condition = "input.globalcor == 'TRUE'",
      column(width=12, 
             id="inputarea",
             column(width=6,numericInput(width="350px","abs_cor_thresh", "Absolute correlation threshold:", 0.4, min = 0, max = 1)),
             column(width=6,numericInput(width="350px","cor_fdrthresh", "FDR threshold for correlation analysis:", 0.05, min = 0, max = 1)),
             column(width=6,selectInput(width="350px","cor_method","Correlation method:",c("spearman","pearson"))),
             column(width=6,selectInput(width="350px","networktype","Network type:",c("complete","GGM")))
      )
    ),
    column(width=12,radioButtons("WGCNAmodules", "3.2 Perform WGCNA module preservation analysis:", inline=TRUE,c(True = "TRUE",False = "FALSE"),selected = "FALSE")),
    column(width=12,radioButtons("globalclustering", "3.3 Perform global clustering analysis (HCA and EM clustering):", inline=TRUE,c(True = "TRUE",False = "FALSE"),selected = "FALSE"))
  )
)

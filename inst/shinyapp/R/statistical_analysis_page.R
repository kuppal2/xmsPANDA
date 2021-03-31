
source("R/input.R")
source("R/data_preprocessing.R")
source("R/workflow.R")
source("R/network_analysis.R")
source("R/graphical_options.R")
source("R/load_packages.R")

# Define UI for data upload app ----

statistical_analysis_page <- fluidPage(
  
  # Make the Layout of parameter ----
  column(12,style="padding-top:10px;",navlistPanel(
    "Input Files",
    tabPanel("Choose Files (see help and support)",input),
    "Parameter Settings for analysis",
    tabPanel("Data preprocessing", data_preprocessing),
    tabPanel("Statistical, feature selection, and predictive modeling",workflow),
    tabPanel("Network analysis",network_analysis),
    tabPanel("Graphical options",graphical_options),
    
    widths = c(3, 9)
  )),

  column(12,style='padding-left:0px;margin-bottom:20px;',
    mainPanel(
      style='padding-left:0px;',
    verbatimTextOutput("nText2"),
    verbatimTextOutput("nText"),
    bsAlert("alert"),
    
    actionButton("go","Start processing",icon=icon("play-circle")),
    downloadButton("downloadData", label = "Download results"),
#  actionButton("stop", "Stop app",icon=icon("stop-circle")),
actionButton("resetAll", "Refresh app",icon=icon("refresh"))
#  actionButton("stop", "Terminate",icon=icon("stop-circle"))
  )),

  #column(12,style='padding-top:5px;padding-left:0;',tags$div(h4("Output"))),
  uiOutput("output_results")
  
)
  



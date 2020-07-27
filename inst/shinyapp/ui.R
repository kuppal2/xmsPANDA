library(shiny)
library(shinyBS)
library(V8)
library(BiocManager)
options(repos = BiocManager::repositories())
library(xmsPANDA)
source("R/introduction_page.R")
source("R/statistical_analysis_page.R")
source("R/additional_analysis.R")
source("R/help_page.R")
jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

ui<-fluidPage(
  shinyjs::useShinyjs(),
# shinyjs::extendShinyjs(script='js/css.js'),
#  shinyjs::extendShinyjs(text = jsResetCode),

    tags$head(
      tags$meta(charset="utf-8"),
      tags$meta(name="description",content="Free Web tutorials"),
      tags$title("xmsPANDA - v1.1.63b"),
      tags$link(rel = "stylesheet", type = "text/css", href = "mystyle.css")
    ),
    column(12,
        tags$div(id="form",
          h3(tags$img(style="float:left;margin-right:15px;",src="images/xmsPANDA_log.png",height='80px',width="80px"),"xmsPANDA - R package for biomarker discovery, supervised and unsupervised learning, and network analysis (v1.1.63b)")),
        tabsetPanel(
         
            tabPanel("Introduction", introduction_page), 
          tabPanel("Statistical Analysis", statistical_analysis_page),
         
          tabPanel("Additional Analysis", additional_analysis_page), 
          tabPanel("Help and Support", help_page),
         
          type ="tabs"
        )
    ),
    column(style="padding-top:0px;padding-bottom:0px;",12,tags$hr(style="margin-top:0px;margin-bottom:15px;border-top: 0.5px solid #ccccb3;")),
 column(12,  tags$div(style="margin-center",tags$footer(align="center",color="white",style="font-weight:normal;font-size:95%;color:black","Maintained by Chunyu Ma & Karan Uppal (",tags$a(href="mailto:kuppal2@emory.edu","kuppal2@emory.edu"),") at Emory University, Atlanta, GA, USA")))
)

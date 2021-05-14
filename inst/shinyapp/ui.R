
source("R/introduction_page.R")
source("R/statistical_analysis_page.R")
source("R/additional_analysis.R")
source("R/help_page.R")
source("R/load_packages.R")

jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

ui<-fluidPage(
  shinyjs::useShinyjs(),
# shinyjs::extendShinyjs(script='js/css.js'),
#  shinyjs::extendShinyjs(text = jsResetCode),

    tags$head(
      tags$meta(charset="utf-8"),
      tags$meta(name="description",content="Free Web tutorials"),
      tags$title("xmsPANDA - 1.3.2"),
      tags$link(rel = "stylesheet", type = "text/css", href = "mystyle.css")
    ),
    column(12,
        tags$div(id="form",
          h3(tags$img(style="float:left;margin-right:15px;",src="images/xmsPANDA_log.png",height='80px',width="80px"),"xmsPANDA (1.3.2) - R package for biomarker discovery,
             supervised and unsupervised learning, and network analysis")),
        tabsetPanel(
         
          tabPanel("Introduction", introduction_page), 
          
          tabPanel("Statistical Analysis", statistical_analysis_page),
         
          tabPanel("Additional Analysis", additional_analysis_page), 
          tabPanel("Help and Support", help_page),
         
          type ="tabs"
        )
    ),
    column(style="padding-top:0px;padding-bottom:0px;",12,tags$hr(style="margin-top:0px;margin-bottom:15px;border-top: 0.5px solid #ccccb3;")),
 column(12,  tags$div(style="margin-center",tags$footer(align="center",color="white",style="font-weight:normal;font-size:95%;color:black","Please ask questions or report any issues on the ",
                                                        tags$a(href='https://github.com/kuppal2/xmsPANDA/issues',target="_blank","GitHub"), " page"))),
column(12,  tags$div(style="margin-center",tags$footer(align="center",color="white",style="font-weight:normal;font-size:95%;color:black","Release date: 05/14/2021"))),
)

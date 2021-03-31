

introduction_page <- fluidPage(
  column(12),
  column(12,style="",
          tags$div(id="introduction_paragraph2",style="padding-top:5px;",tags$p(align="left",style="font-weight:bold;font-size:120%;color:darkblue","xmsPANDA provides an automated workflow including data preprocessing, statistical analysis, robust feature selection, and network analysis
                  for metabolomics analysis"))
                   #tags$p(align="center",style="font-weight:bold;font-size:80%;color:black","For installing xMWAS locally in R run: ")
                   #tags$p(align="center",style="font-weight:bold;font-size:80%;color:darkred","library(devtools);install_github(\"kuppal2/xMWAS\")")
         ),
  column(12,style="margin-bottom:10px",tags$img(style="margin-left:auto;margin-right:auto;display: block;",src="images/Figure_homepage.png",height='400px',width="70%"))
)

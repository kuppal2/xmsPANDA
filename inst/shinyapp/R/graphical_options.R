library(shiny)

graphical_options<-fluidRow(
  tags$div(
    id="maindiv",
    column(width=12,
           column(width=6,selectInput(width="350px","heatmap_color_scheme","Heatmap color palettes/scheme:",c("redblue","yellowblue","redyellowgreen","yellowwhiteblue","redwhiteblue","topo","heat"))),
           #column(width=6,selectInput(width="350px","output_format","Output format (PDF or PNG for high-res):",c("png","pdf"))),
           column(width=6,numericInput(width="350px","pca_cex_val", "Size of points on PCA plots (1-20 limit):", 4, min = 1, max = 20))
    ),
    column(width=12,
           column(width=6,numericInput(width="350px","ellipse_conf_level", "Confidence interval for PCA ellipses (0-1 limit):", 0.95, min = 0, max = 1)),
           column(width=6,selectInput(width="350px","pca_ellipse","Should ellipse be plotted on PCA plots?",c("TRUE","FALSE")))
    ),
    column(width=12,
              column(width=4,selectInput(width="350px","boxplot.type", "Boxplot type:", c("ggplot","simple"))),
                

column(width=4,selectInput(width="350px","boxplot.jitter", "Add jitter to boxplots:", c("TRUE","FALSE"))),
column(width=4,selectInput(width="350px","boxplot.pvalues", "Add p-values to boxplots:", c("TRUE","FALSE")))
#column(width=6,selectInput(width="350px","plot.color.theme","Color theme for boxplots, barplots, and line plots: ",c("journal","default","topo","heat","rainbow","terrain","black","grey57")))

           
),
column(width=12,
column(width=6,
                  id="inputarea",
           textInput(width="350px","boxplot.color.theme", "Color theme or options for boxplots (e.g. heat, rainbow, red, grey57):","journal",placeholder="Default: journal")


           ),
column(width=6,
                  id="inputarea",
           textInput(width="350px","barplot.color.theme", "Color theme or options for barplots:","journal",placeholder="Default: journal")


           )
),

column(width=12,
         
column(width=6,
       id="inputarea",
textInput(width="350px","sample.color.theme", "Color theme or options for samples:","journal",placeholder="Default: journal")


),
column(width=6,
                  id="inputarea",
           textInput(width="350px","lineplot.color.theme", "Color theme or options for lineplots:","journal",placeholder="Default: journal")


           )

       ),
column(width=12,
         
column(width=6,
       id="inputarea",
selectInput(width="350px","timeseries.lineplots", "Plot time series lineplots (for time-series data):",c("FALSE","TRUE"))

),
column(width=6,
       id="inputarea",
selectInput(width="350px","alphabetical.order", "Plot classes on the x-axis in alphabetical order:",c("TRUE","FALSE"))

),
column(width=6,
       id="inputarea",
textInput(width="350px","ylabel.text", "Label for y-axis in boxplots, barplots, and lineplots:","Abundance",placeholder="Default: Abundance")

),


  )
))

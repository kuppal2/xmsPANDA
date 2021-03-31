

graphical_options<-fluidRow(
  tags$div(
    id="maindiv",
    column(width=12,
           column(width=6,selectInput(width="350px","heatmap_color_scheme","Color theme for up or down status in heatmap/Manhattan plot/Volcano plot:",c("redblue","orangeblue","bluered","blueorange"),selected="redblue")),
           #column(width=6,selectInput(width="350px","output_format","Output format (PDF or PNG for high-res):",c("png","pdf"))),
           #column(width=6,numericInput(width="350px","pca_cex_val", "Size of points on PCA plots (1-20 limit):", 4, min = 1, max = 20))
           column(width=6,selectInput(width="350px","color.palette",
                                      tags$p(style="font-size: 14px;","Color theme for samples in PCA, HCA, boxplots, and line plots. See the palettes ",
                                                                       tags$a(href='https://github.com/kuppal2/xmsPANDA/blob/master/inst/color_palettes_xmsPANDA.pdf',target="_blank","here")),
                                      c("journal","wong","npg","nejm","jco","lancet","custom1","brewer.RdYlBu","brewer.RdBu","brewer.PuOr","brewer.PRGn","brewer.PiYG","brewer.BrBG",
                                        "brewer.Set2","brewer.Paired","brewer.Dark2","brewer.YlGnBu","brewer.YlGn","brewer.YlOrRd","brewer.YlOrBr","brewer.PuBuGn",
                                        "brewer.PuRd","brewer.PuBu",
                                        "brewer.OrRd","brewer.GnBu","brewer.BuPu","brewer.BuGn","brewer.blues","black","grey65",
                                        "topo"),
                                      selected="journal"))
           
    ),
    column(width=12,
           column(width=6,selectInput(width="350px","hca.labRow.value","Print row names in heatmaps:",c("TRUE","FALSE"),
                                                                                                              selected="TRUE")),
           #column(width=6,selectInput(width="350px","output_format","Output format (PDF or PNG for high-res):",c("png","pdf"))),
           #column(width=6,numericInput(width="350px","pca_cex_val", "Size of points on PCA plots (1-20 limit):", 4, min = 1, max = 20))
           column(width=6,selectInput(width="350px","hca.labCol.value","Print column names in heatmaps: ",c("TRUE","FALSE"),
                                      selected="TRUE"))
           
    ),
    column(width=12,
        #   column(width=6,numericInput(width="350px","ellipse_conf_level", "Confidence interval for PCA ellipses (0-1 limit):", 0.95, min = 0, max = 1)),
          column(width=6,selectInput(width="350px","pca_ellipse","Should ellipse be plotted on PCA plots?",c("TRUE","FALSE"),selected="TRUE")),
        column(width=6,numericInput(width="350px","hca.cex.legend", "Size of HCA legend text (set to -1 to turn off):", 0.7, min = 0.001, max = 20)),
         # column(width=6,numericInput(width="350px","pca_cex_val", "Size of points on PCA plots (1-20 limit):", 4, min = 1, max = 20))
    ),
    column(width=12,
              column(width=6,selectInput(width="350px","boxplot.type", "Boxplot type:", c("ggplot","simple"))),
              column(width=6,numericInput(width="350px","cex.plots", "Size of x and y-axis text in the plots (used for all plots):", 0.6, min = 0.001))
),
   column(width=12,
column(width=6,selectInput(width="350px","boxplot.jitter", "Add jitter to boxplots:", c("FALSE","TRUE"))),
column(width=6,selectInput(width="350px","boxplot.pvalues", "Add p-values to boxplots:", c("FALSE","TRUE")))

           
),
column(width=12,
#column(width=6,
 #                 id="inputarea",
           #textInput(width="350px","ggplot.type1", "Color theme or options for boxplots (e.g. heat, rainbow, red, grey57):","journal",placeholder="Default: journal")
column(width=6,selectInput(width="350px","ggplot.type1", "Select how to group samples by factor1 or factor2 in boxplots:", 
                           c("Group/facet by factor 2"="TRUE","Group/facet by factor 1"="FALSE","No facet or sub-grouping"="None"))),
column(width=6,numericInput(width="350px","facet.nrow", "Number of rows in the boxplots (only applicable to factorial designs):", 1, min = 1, max = 10))
),

column(width=12,
         
column(width=6,
       id="inputarea",
selectInput(width="350px","timeseries.lineplots", "Plot time series lineplots (for time-series data):",c("FALSE","TRUE"))

),
column(width=6,
       id="inputarea",
selectInput(width="350px","alphabetical.order", "Plot classes on the x-axis in alphabetical order:",c("Automatically sort by alphabetical order"="TRUE",
                                                                                                      "Use the order in the class labels file"="FALSE"),selected="FALSE")

)),
column(width=12,
      
column(width=6,
       id="inputarea",
textInput(width="350px","ylabel.text", "Label for y-axis in boxplots, barplots, and lineplots:","Abundance",placeholder="Default: Abundance")),

column(width=6,
       id="inputarea",
       selectInput(width="350px","plot.boxplots.raw", "Plot boxplots using raw/non-transformed data:",c("FALSE","TRUE"))
       
)
),
column(width=12,
       #column(width=6,
       #                 id="inputarea",
       #textInput(width="350px","ggplot.type1", "Color theme or options for boxplots (e.g. heat, rainbow, red, grey57):","journal",placeholder="Default: journal")
       column(width=6,numericInput(width="350px","plots.height", "Height of the PDF/PNG files in inches:", 8, min = 1, max = 20)),
       column(width=6,numericInput(width="350px","plots.width", "Width of the PDF/PNG files in inches:", 8, min = 1, max = 20))
)


  
))

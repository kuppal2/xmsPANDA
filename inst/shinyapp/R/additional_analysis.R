source("R/load_packages.R")
# Define UI for data upload app ----

additional_analysis_page <- fluidPage(
  
  # Make the Layout of parameter ----
    tags$div(
      id="maindiv4",
      style='margin-top:12px',
      bsCollapse(open = "",
       

#get_pcascoredistplots
bsCollapsePanel("Principal Component Analysis",column(width=12,
            column(width=6,
                   id="pcainputarea",
                   fileInput("pcainput1", "Select your feature table file ('.csv' or .txt')",
                             multiple = FALSE,
            #  width="50px",
                             accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))),
            column(width=6,id="pcainputarea",
                   fileInput("pcainput2", "Select your class labels file ('.csv' or .txt')",
                   multiple = FALSE,
            # width="50px",
                   accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))),
            
            ),column(width=12,
                     
           # column(width=6,
            #       id="pcainputarea",
            #textInput(width="350px","pca.sample.color.theme", "Color theme or options for samples:","journal",placeholder="Default: journal")
            
            
            #),
           column(width=6,selectInput(width="350px","pca.sample.color.theme",
                                      tags$p(style="font-size: 14px;","Color theme. See the palettes ",
                                             tags$a(href='https://github.com/kuppal2/xmsPANDA/blob/master/inst/color_palettes_xmsPANDA.pdf',target="_blank","here")),
                                      c("journal","wong","npg","nejm","jco","lancet","custom1","brewer.RdYlBu","brewer.RdBu","brewer.PuOr","brewer.PRGn","brewer.PiYG","brewer.BrBG",
                                        "brewer.Set2","brewer.Paired","brewer.Dark2","brewer.YlGnBu","brewer.YlGn","brewer.YlOrRd","brewer.YlOrBr","brewer.PuBuGn",
                                        "brewer.PuRd","brewer.PuBu",
                                        "brewer.OrRd","brewer.GnBu","brewer.BuPu","brewer.BuGn","brewer.blues","black","grey65",
                                        "topo"),
                                      selected="journal")),
            column(width=6,
                              id="pcainputarea",
                     #  textInput(width="350px","pca.lineplot.color.theme", "Color theme or options for lineplots:","journal",placeholder="Default: journal")
                     selectInput(width="350px","pca.ellipse", "Plot PCA ellipse:",c("FALSE","TRUE"))
            
            
                       )
            
                   ),
            column(width=12,
                     
            column(width=6,
                   id="pcainputarea",
            selectInput(width="350px","pca.timeseries.lineplots", "Plot time series lineplots (for time-series data):",c("FALSE","TRUE"))
            
            ),
            column(width=6,
                   id="pcainputarea",
            selectInput(width="350px","pca.alphabetical.order", "Plot classes on the x-axis in alphabetical order:",c("TRUE","FALSE"))
            
            ),
            ),column(width=12,
                     
            column(width=6,
                   id="pcainputarea",
            selectInput(width="350px","pca.analysistype", "Analysis type:",c("multiclass or one-way ANOVA"="oneway","one-way ANOVA with repeat measures"="onewayrepeat","two-way ANOVA with repeat measures"="twowayrepeat","two-way ANOVA"="twoway","continuous outcome or regression"="regression"))
            
            ),
            column(width=6,
                   id="pcainputarea",
            selectInput(width="350px","pca.pairedanalysis", "Paired analysis (time-series):",c("FALSE","TRUE"))
            
            )
            ),
            column(width=12,
            column(width=6,
            id="pcainputarea",actionButton("pcastart","Start processing",icon=icon("play-circle")),
            downloadButton("downloadpcaData", label = "Download results"))),
            column(width=12,verbatimTextOutput("pcaText4")),
            column(width=12,verbatimTextOutput("pcaText5")),
            #textOutput("pcaoutput_results"),
              style='primary'
            ),



#get_boxplotscoredistplots
bsCollapsePanel("Boxplots",column(width=12,
                                                      column(width=6,
                                                             id="boxplotinputarea",
                                                             fileInput("boxplotinput1", "Select your feature table file ('.csv' or .txt')",
                                                                       multiple = FALSE,
                                                                       #  width="50px",
                                                                       accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))),
                                                      column(width=6,id="boxplotinputarea",
                                                             fileInput("boxplotinput2", "Select your class labels file ('.csv' or .txt')",
                                                                       multiple = FALSE,
                                                                       # width="50px",
                                                                       accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))),
                                                      
),column(width=12,
         
        # column(width=6,
         #       id="boxplotinputarea",
          #      textInput(width="350px","boxplot.sample.color.theme", "Color theme or options for samples:","journal",placeholder="Default: journal")
                
                
         #),
         
         column(width=6,selectInput(width="350px","boxplot.sample.color.theme",
                                    tags$p(style="font-size: 14px;","Color theme. See the palettes ",
                                           tags$a(href='https://github.com/kuppal2/xmsPANDA/blob/master/inst/color_palettes_xmsPANDA.pdf',target="_blank","here")),
                                    c("journal","wong","npg","nejm","jco","lancet","custom1","brewer.RdYlBu","brewer.RdBu","brewer.PuOr","brewer.PRGn","brewer.PiYG","brewer.BrBG",
                                      "brewer.Set2","brewer.Paired","brewer.Dark2","brewer.YlGnBu","brewer.YlGn","brewer.YlOrRd","brewer.YlOrBr","brewer.PuBuGn",
                                      "brewer.PuRd","brewer.PuBu",
                                      "brewer.OrRd","brewer.GnBu","brewer.BuPu","brewer.BuGn","brewer.blues","black","grey65",
                                      "topo"),
                                    selected="journal")),
         column(width=6,
                id="boxplotinputarea",
                selectInput(width="350px","boxplot.type", "Select function to generate boxplots:",c("boxplot (simple)"="simple","ggplot"="ggplot"),selected="ggplot")
                
                
         )
         
),
column(width=12,
       
       column(width=6,
              id="boxplotinputarea",
              selectInput(width="350px","boxplot.ggplot.type1", "Group (facet) by factor 2:",c("TRUE","FALSE"))
              
       ),
       column(width=6,
              id="boxplotinputarea",
              selectInput(width="350px","boxplot.alphabetical.order", "Plot classes on the x-axis in alphabetical order:",c("TRUE","FALSE"))
              
       )
),column(width=12,
         
         column(width=6,
                id="boxplotinputarea",
                selectInput(width="350px","boxplot.analysistype", 
                            "Analysis type:",c("multiclass or one-way ANOVA"="oneway","one-way ANOVA with repeat measures"="onewayrepeat",
                                               "two-way ANOVA with repeat measures"="twowayrepeat","two-way ANOVA"="twoway"))
                
         ),
        
         column(width=6,
                id="inputarea",
                textInput(width="350px","boxplot.ylabel.text", "Label for y-axis in boxplots:","Abundance",placeholder="Default: Abundance"))
         
         
),
column(width=12,
       
       column(width=6,
              id="boxplotinputarea",
              selectInput(width="350px","boxplot_jitter_1", "Add jitter:",c("TRUE","FALSE"),selected="TRUE")
              
       ),
       column(width=6,
              id="boxplotinputarea",
              selectInput(width="350px","boxplot_pvalues_1", "Add p-values:",c("TRUE","FALSE"),selected="TRUE")
              
       )
),
column(width=12,
       #column(width=6,
       #                 id="inputarea",
       #textInput(width="350px","ggplot.type1", "Color theme or options for boxplots (e.g. heat, rainbow, red, grey57):","journal",placeholder="Default: journal")
       column(width=6,numericInput(width="350px","boxplot.plots.height", "Height of the PDF file in inches:", 8, min = 1, max = 20)),
       column(width=6,numericInput(width="350px","boxplot.plots.width", "Width of the PDF file in inches:", 8, min = 1, max = 20))
),
column(width=12,
       #column(width=6,
       #                 id="inputarea",
       #textInput(width="350px","ggplot.type1", "Color theme or options for boxplots (e.g. heat, rainbow, red, grey57):","journal",placeholder="Default: journal")
       column(width=6,numericInput(width="350px","boxplot.min.ylim", "Lower limit on y-axis (set to -1 for auto-select):", -1)),
       column(width=6,numericInput(width="350px","boxplot.max.ylim", "Upper limit on y-axis (set to -1 for auto-select):", -1))
),
column(width=12,
       column(width=6,
              id="boxplotinputarea",actionButton("boxplotstart","Start processing",icon=icon("play-circle")),
              downloadButton("downloadboxplotData", label = "Download results"))),
column(width=12,verbatimTextOutput("boxplotText4")),
column(width=12,verbatimTextOutput("boxplotText5")),
#textOutput("boxplotoutput_results"),
style='primary'
),




  
###Functional Class Scoring
bsCollapsePanel("Functional Class Scoring",
  column(width=12,
         column(width=6,
                id="inputarea",
                fileInput("clusterinput", "Upload your target list ('.csv' or .txt')",
                          multiple = FALSE,
                          width="350px",
                          accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))
         ),
         column(width=6,selectInput(width="350px","fcs.database","Select a database:",c("KEGG Pathways (default)"="pathway",
                                                                                        "KEGG Modules"="module",
                                                                                        "KEGG Brite"="brite",
                                                                                        "RefMet:SuperClass"="refmet_superclass",
                                                                                        "RefMet:MainClass"="refmet_mainclass",
                                                                                        "RefMet:SubClass"="refmet_subclass",
                                                                                        "LipidMaps:MainClass"="lipidmaps_mainclass",
                                                                                        "LipidMaps:SubClass"="lipidmaps_subclass",
                                                                                        "Reactome:Atlas"="reactome_atlas",
                                                                                        "KEGG:Atlas"="kegg_atlas",
                                                                                        "Custom (user-defined)"="custom"),selected=c("pathway")))           
            
),
column(width=12,
       column(width=6,
              conditionalPanel(
                condition = "output.checkfcsdatabase",
                fileInput("clusterinput2", "Upload your custom reference database ('.csv' or .txt')",
                          multiple = FALSE,
                          width="350px",
                          accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))
              )
              
       ),
       column(width=6,
                       conditionalPanel(
                            condition = "output.checkannotfile",
                                   fileInput("clusterinput3", "Upload your annotation file ('.csv' or .txt')",
                                             multiple = FALSE,
                                             width="350px",
                                             accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))
                            )
                                 
            )


),

column(width=12,column(width=6,selectInput(width="350px","type.statistic","Statistic type is p-value?",c("TRUE","FALSE")))),
column(width=12,
         column(width=6,
                style='margin-top:24px',
                actionButton("optionbutton2", "More options"),
                bsModal("moreoptions2", "More options for analysis", "optionbutton2", size = "large",
                        tags$div(
                          width=12,
                          style="height:300px",
                          
                          
                          column(width=6,selectInput(width="350px","kegg_species_code","Select a species (only for KEGG pathway):",c("Homo sapiens(default)",
                          "Mus musculus",
                          "Pan troglodytes",
                          "Macaca mulatta",
                          "Bos taurus",
                          "Rattus norvegicus",
                          "Danio rerio",
                          "C. elegans",
                          "Drosophila melanogaster"))),
#   column(width=6,selectInput(width="350px","count.unique.formula.overlapsize","Only count unique formula for overlap size (for redundant matches):", c("TRUE","FALSE"))),
                          column(width=6,numericInput(width="350px","fcs.min.hits","Minimum hits in a pathway or custom set:", value = 2, min = 1, max = 100000)),
                            column(width=6,selectInput(width="350px","fcs.upload.annotation.file","Upload annotation file if the input data is not mapped to database IDs:", c("FALSE","TRUE"))),
                              
# column(width=6,selectInput(width="350px","fcs.method","FCS method:", c("Z.score"="zscore","Permutation test"="permutation"))),
  column(width=6,selectInput(width="350px","path.bg.color","KEGG pathway background color:", c("blue","red","orange","green","brown","purple","pink","grey","black","white"))),
                            column(width=6,numericInput(width="350px","fcs.itrs","Number of iterations for the permutation method:", value = 100, min = 100, max = 10000))
                          
                        )
                ),
                actionButton("start2","Start processing",icon=icon("play-circle"))
                
         )),
column(width=12,
        column(width=6,
        style='margin-top:24px',align="right",actionButton("optionbutton3", "View results")),
        column(width=6,
        style='margin-top:24px',align="right",actionButton("optionbutton3B", "Help"))
),
        
  column(width=12,verbatimTextOutput("nText4A")),
  column(width=12,verbatimTextOutput("nText5A")),
  
  column(width=12,
tags$head(tags$style("#moreoptions3 .modal-dialog{ width:1250px}")),
bsModal("moreoptions3", "View results", "optionbutton3", size = "big",
conditionalPanel(condition = "output.checktable1",

                          uiOutput("downloadbutton"), div(DT::dataTableOutput("pathwaytb"), style = "font-size: 90%; width: 20%") #DTOutput('pathwaytb',width=6)
#,tableOutput('pathwaytb')

))
        
  ),
column(width=12,style="padding-left:25px;padding-right:25px;",

        bsModal("moreoptions3B", "Help", "optionbutton3B", size = "large",
    conditionalPanel(condition = "",
tags$div(tags$h5(style='font-weight:bold',"Description:")),
#tags$p(style='font-weight:normal',"Functional class scoring determines enrichment of functional classes or sets (e.g. pathways, modules, classes) in 
 #      association with a phenotype based on the p-value, VIP, fold change, and other statistics. 
  #     In stage one, statistics of individual entities (features, metabolites) within each functional class are combined to assign an 
   #    aggregated statistic at the functional class level. Only the IDs with unique chemical formula are used for aggregation within each functional class. 
    #   In case of redundancies, the entry with the highest statistic is used. The software applies three methods for aggregation: 
     #  1) mean statistic calculated as sum of the statistic of ; 
      # 2) z-score calculated as sqrt(number of entities in a set)*average value of the statistic (based on Irizarry 2009; PMC3134237);
      # 3) max-mean method by Efron and Tibshiranil (Efron, Bradley; Tibshirani, Robert. On testing the significance of sets of genes. Ann. Appl. Stat. 1 (2007), no. 1, 107--129. doi:10.1214/07-AOAS101. https://projecteuclid.org/euclid.aoas/1183143731). 
      # For each method, a permutation test at the functional class/set level is used to determine p-values. 
      # A meta-pvalue is assigned by aggregating the p-values from permutation tests for each method using the Fisher's p-value aggregation method (Fisher, R. A., 1932. Statistical Methods for Research Workers, 4th Edition. Oliver and Boyd, Edinburgh.)
      # Benjamini-Hochberg method is used for false discovery correction. The input file can be provided in three formats as described below. For the formats where the input file does not contain database identifiers, please see formats 2 and 3 and the corresponding annotation files as a separate input."),
tags$p(style='font-weight:normal',"Functional class scoring determines enrichment of functional classes or sets (e.g. pathways, modules, classes) in 
      association with a phenotype based on the p-value, VIP, fold change, and other statistics. 
       In stage one, statistics of individual entities (features, metabolites) within each functional class are combined to assign an 
       summary statistic at the functional class level using: 1) z-score calculated as sqrt(number of entities in a set)*average value of the statistic (based on Irizarry 2009; PMC3134237); 2) max-mean method developed by Efron and Tibshiranil (Efron, Bradley; Tibshirani, Robert. On testing the significance of sets of genes. Ann. Appl. Stat. 1 (2007), no. 1, 107--129. 
       doi:10.1214/07-AOAS101. https://projecteuclid.org/euclid.aoas/1183143731). Only the IDs with unique chemical formula are used for aggregation within each functional class. 
       In case of redundancies, the entry with the highest statistic is used. The p-values are calculated using the one-sample z-test for the z-score method and using the permutation test by 
       shuffling the functional sets for the max-mean method. Benjamini-Hochberg method is used for false discovery correction. The input file can be provided in three formats as described below. For the formats where the input file does not contain database identifiers, please see formats 2 and 3 and the corresponding annotation files as a separate input."),
tags$div(tags$h5(style='font-weight:bold',"Input file format(s) for target IDs file:")),
                                tags$p(style='font-weight:normal',"Format 1) column A: XID (e.g. KEGG ID); column B: Statistic (p-value,vip,beta coefficient,fold change)"),
                                tags$p(style='font-weight:normal',"Example:"),
                                tags$table(tags$tr(tags$th("XID"),tags$th("Statistic")),
                                           tags$tr(tags$td("C00019"),tags$td("0.2368648")),
                                           tags$tr(tags$td("C00491"),tags$td("0.2319177")),
                                           tags$tr(tags$td("C00300"),tags$td("0.8520213"))
                                ),
                                tags$p(style='font-weight:normal',"Format 2) column A: mz; column B: time; column C: Statistic (p-value,vip,beta coefficient,fold change)"),
                                                               tags$p(style='font-weight:normal',"Example:"),
                                                               tags$table(tags$tr(tags$th("mz"),tags$th("time"),tags$th("Statistic")),
                                                                          tags$tr(tags$td("148.0678"),tags$td("55"),tags$td("0.2368648")),
                                                                          tags$tr(tags$td("150.0580"),tags$td("65"),tags$td("0.2319177")),
                                                                          tags$tr(tags$td("166.0868"),tags$td("80"),tags$td("0.8520213"))
                                                               ),
                                tags$p(style='font-weight:normal',"Format 3) column A: Name; column B: Statistic (p-value,vip,beta coefficient,fold change)"),
                                                                                             tags$p(style='font-weight:normal',"Example:"),
                                                                                             tags$table(tags$tr(tags$th("Name"),tags$th("Statistic")),
                                                                                                        tags$tr(tags$td("Glutamate"),tags$td("0.2368648")),
                                                                                                        tags$tr(tags$td("Methionine"),tags$td("0.2319177")),
                                                                                                        tags$tr(tags$td("Phenylalanine"),tags$td("0.8520213"))
                                                                                             )

    ),
  conditionalPanel(condition = "",
                                       tags$div(tags$h5(style='font-weight:bold',"Input file format for custom reference database:")),
                                       tags$p(style='font-weight:normal',"The input file should have fvee columns: 1. XID (e.g. KEGG ID); 2. set ID (e.g. KEGG pathway map ID); 3. set name (e.g. pathway name); 4. Optional: Exact mass (can assign a value of 1 if mass unknown); 5. Optional: Chemical formula (can assign a unique alphanumeric value if formula is unknown"),
                                       tags$p(style='font-weight:normal',"Example:"),
                                       tags$table(tags$tr(tags$th("XID"),tags$th("SetID"),tags$th("Name"),tags$th("ExactMass (optional)"),tags$th("Formula (optional)")),
                                                  tags$tr(tags$td("C00217"),tags$td("map00270"),tags$td("Cysteine and methionine metabolism")),
                                                  tags$tr(tags$td("C00073"),tags$td("map00270"),tags$td("Cysteine and methionine metabolism")),
                                                  tags$tr(tags$td("C00079"),tags$td("map00330"),tags$td("Arginine and proline metabolism"))
                                       )
   ),
        conditionalPanel(condition = "",
                                            tags$div(tags$h5(style='font-weight:bold',"Input file formats for annotation file (only applicable to formats 2 and 3 above:")),
                                            tags$p(style='font-weight:normal',"Annotation file for format 2) column A: XID (e.g. KEGG ID); column B: mz; column C: time"),
                                            tags$p(style='font-weight:normal',"Example:"),
                                            tags$table(tags$tr(tags$th("XID"),tags$th("mz"),tags$th("time")),
                                                       tags$tr(tags$td("C00217"),tags$td("148.0678"),tags$td("55"),),
                                                       tags$tr(tags$td("C00073"),tags$td("150.0580"),tags$td("65")),
                                                       tags$tr(tags$td("C00079"),tags$td("166.0868"),tags$td("80"))
                                            ),
                                            tags$p(style='font-weight:normal',"Annotation file for format 3) column A: XID (e.g. KEGG ID); column B: Name"),
                                           tags$p(style='font-weight:normal',"Example:"),
                                           tags$table(tags$tr(tags$th("XID"),tags$th("Name")),
                                                      tags$tr(tags$td("C00217"),tags$td("Glutamate")),
                                                      tags$tr(tags$td("C00073"),tags$td("Methionine")),
                                                      tags$tr(tags$td("C00079"),tags$td("Phenylalanine"))
                                           )
        )

)
                  
        ),
  #column(width=12,verbatimTextOutput('hover')),
  style='primary'
),

bsCollapsePanel("Compare normalization methods",
                column(width=12,
                       column(width=6,
                              id="inputarea2",
                              fileInput("norminput1", "Select your feature table file ('.csv' or .txt')",
                                        multiple = FALSE,
                                        accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))),
                       column(width=6,
                              id="inputarea2",
                              fileInput("norminput2", "Select your class labels file ('.csv' or .txt')",
                                        multiple = FALSE,
                                        accept = c("text/csv","text/comma-separated-values,text/plain",".txt")))
                ),
                column(width=6,selectInput(width="350px","study.design","Select study design:",c("multiclass or one-way ANOVA"="multiclass","one-way ANOVA with repeat measures"="onewayrepeat","two-way ANOVA with repeat measures"="twowayanovarepeat","two-way ANOVA"="twowayanova"))),
                column(width=6,selectInput(width="350px","analysistype","Analysis type:",c("Classification"="classification","Regression"="regression"))),
                column(width=6,selectInput(width="350px","normalization.method","Select normalization methods (click in the box to select more options):",multiple=TRUE,c("all(default)"="all",
                                                                                                                                                                          "log2 and quantile normalization"="log2quantilenorm",
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
                                                                                                                                                                          "cubicspline normalization"="cubicspline_norm", "median absolute deviation normalization"="mad_norm"),selected=c("log2transform","quantile_norm","log2quantilenorm","znormtransform","mstus","eigenms_norm"))),
                column(width=6,selectInput(width="350px","summary_na_replacement","Choose an imputation method:",c("halffeaturemin","knn","randomforest","zeros","halfsamplemin","halfdatamin","none"))),
                column(width=6,numericInput(width="350px","log2.transform.constant","Small constant to add for log2 transformation:",1)),
                column(width=6,
                       style='margin-top:24px',
                       actionButton("normoptionbutton2", "More options"),
                       bsModal("normmoreoptions2", "More options for analysis", "normoptionbutton2", size = "large",
                               tags$div(
                                 width=12,
                                 style="height:300px",
                                 column(width=12,mainPanel(
                                   width=12,bsCollapsePanel("Replicate summarization",
                                                            column(width=6,numericInput(width="350px","numreplicate", "Number of technical replicates (1-10 limit):", 1, min = 1, max = 10)),
                                                            column(width=6,selectInput(width="350px","summarization_method","Choose a replicate summarization method:",c("median","mean"))),
                                                            column(width=6,radioButtons("use_summarizion", "Should the replicates be summarized?", inline=TRUE,c(True = "TRUE",False = "FALSE"),selected = "TRUE")),
                                                            column(width=6,numericInput(width="350px","summarization_ratio", "Maximum missing value ratio:", 0.3, min = 0, max = 1)),
                                                            bsTooltip("summarization_ratio", "What propotion of replicates are allowed to have missing values during the averaging or median summarization step of each biological sample?","bottom"),
                                                            style = "primary"),
                                   
                                   bsCollapsePanel("Filtering",
                                                   column(width=6,numericInput(width="350px","all_missing_thresh", "Minimum non-missing sample ratio:", 0.5, min = 0, max = 1)),
                                                   bsTooltip("all_missing_thresh", "What propotion of total number of samples should have an intensity?","bottom"),
                                                   column(width=6,numericInput(width="350px","rsd_filt_list", "Minimum overall variance:", 0, min = 0, max = 100)),
                                                   bsTooltip("rsd_filt_list", "Minimum relative standard deviation across all samples","bottom"),
                                                   column(width=6,numericInput(width="350px","group_missing_thresh", "Minimum non-missing sample ratio for group:", 0.8, min = 0, max = 1)),
                                                   bsTooltip("group_missing_thresh", "Minimum propotion of samples in at least one group in which a non-missing signal value should be present","bottom"),
                                                   style = "primary"),
                                   
                                   bsCollapsePanel("Correlation analysis",
                                                   column(width=6,numericInput(width="350px","abs.cor.thresh", "Absolute correlation threshold:", -1, min = 0, max = 1)),
                                                   bsTooltip("abs.cor.thresh", "Correlation threshold for pairwise correlations between features","bottom. Set to -1 to turn off correlation analysis."),
                                                   column(width=6,numericInput(width="350px","pvalue.thresh", "p-value threshold for correlations:", 0.05, min = 0, max = 1)),
                                                   bsTooltip("pvalue.thresh", "P-value cutoff based on Student's t-test","bottom"),
                                                   column(width=6,numericInput(width="350px","cor.fdrthresh", "FDR threshold for correlations:", 0.05, min = 0, max = 1)),
                                                   
                                                   style = "primary"),
                                 ))
                                 
                                 
                                 
                               )
                               
                       ),
                       actionButton("normstart","Start processing",icon=icon("play-circle")),
                       downloadButton("downloadnormData", label = "Download results")
                       
                ),
                
                column(width=12,verbatimTextOutput("normText5")),
                
                
                style='primary'
),

bsCollapsePanel("Metabolite Quantification Analysis",
                column(width=12,
                       column(width=4,
                              id="inputarea",
                              fileInput("featuretable_file", "Select feature table file ('.csv' or .txt', 100MB limit)",
                                        multiple = FALSE,
                                        width="350px",
                                        accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))
                       ),
                       column(width=4,
                              id="inputarea",
                              fileInput("classlabel_file", "Select class label file ('.csv' or '.txt')",
                                        multiple = FALSE,
                                        width="350px",
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".txt"))
                       ),
                       column(width=4,
                              style='margin-top:24px',
                              actionButton("optionbutton4", "More options"),
                              bsModal("moreoptions4", "More options for analysis", "optionbutton4", size = "large",
                                      tags$div(
                                        width=12,
                                        style= "height:660px",
                                        column(width=12, tags$p(style="font-weight:bold;font-size:16px;","Select the step you need to run:")),
                                        column(width=12, style='margin-bottom:10px',
                                               prettyCheckbox(inputId = "step1",
                                                              label = "Step1: Draw distribution plots",
                                                              value = TRUE,
                                                              thick = TRUE,
                                                              shape = "round",
                                                              animation = "smooth",
                                                              status = "primary",
                                                              inline = TRUE),
                                               prettyCheckbox(inputId = "step2",
                                                              label = "Step2: Quantify concentration",
                                                              value = TRUE,
                                                              thick = TRUE,
                                                              shape = "round",
                                                              animation = "smooth",
                                                              status = "primary",
                                                              inline = TRUE),
                                               prettyCheckbox(inputId = "step3",
                                                              label = "Step3: Pull the KEGG map from KEGG database",
                                                              value = FALSE,
                                                              thick = TRUE,
                                                              shape = "round",
                                                              animation = "smooth",
                                                              status = "primary",
                                                              inline = TRUE)
                                        ),
                                        column(width=12,
                                               column(width=6,selectInput(width="350px","summarize_replicates","Summarize technical replicate?",c("TRUE","FALSE"))),
                                               column(width=6,
                                                      conditionalPanel(condition = "input.summarize_replicates == 'TRUE'",
                                                                       numericInput(width="350px","num_replicate2", "Number of technical replicates  (1-10 limit):", 1, min = 1, max = 10)
                                                      )
                                               )
                                        ),
                                        column(width=12,
                                               column(width=6,
                                                      conditionalPanel(condition = "input.summarize_replicates == 'TRUE'",
                                                                       numericInput(width="350px","rep_max_missing_thresh", "Maximum missing value ratio:", 0.3, min = 0, max = 1)
                                                      )
                                               ),
                                               column(width=6,
                                                      conditionalPanel(condition = "input.summarize_replicates == 'TRUE'",
                                                                       selectInput(width="350px","summary_method","Choose a replicate summarization method:",c("median","mean"))
                                                      )
                                               )
                                        ),
                                        column(width=12,
                                               column(width=6,numericInput(width="350px","mass_error", "Mass-to-charge tolerance(ppm)  (0-100 limit):", 10, min = 0, max = 100)),
                                               column(width=6,numericInput(width="350px","time_error", "Retention time tolerance(second)  (0-1000 limit):", 30, min = 0, max = 1000))
                                        ),
                                        conditionalPanel(condition = "input.step1 || input.step3",
                                                         column(style="padding:0px;",12,tags$hr(style="border-top: 1px solid #000000;"))
                                                         
                                        ),
                                        conditionalPanel(condition = "input.step1",
                                                         column(width=12, tags$p(style="font-weight:bold;font-size:16px;","Parameters for step1:")),
                                                         column(width=12,
                                                                column(width=6,selectInput(width="300px","groupcheck","More than one group in your sample?",c("FALSE","TRUE"))),
                                                                column(width=6,textInput(width="300px","targetID", "Target IDs:","",placeholder="Default: None"))
                                                         ),
                                                         column(width=12,tags$p(style="color:red","For the target ids, you can enter the sample id that you want to highlight in the sample distribution plot. Each
                                                                                sample id should be separated by comma, e.g. sample1,sample2"))
                                                         ),
                                        conditionalPanel(condition = "input.step3",
                                                         column(width=12, tags$p(style="font-weight:bold;font-size:16px;","Parameters for step3:")),
                                                         column(width=12,
                                                                column(width=6,numericInput(width="300px","foldchange_thresh", "Fold Change Threshold (1-100 limit):", 2, min = 1, max = 100)),
                                                                column(width=6,numericInput(width="300px","minhit", "Minimum #metablites hitted in KEGG map:", 3, min = 1, max = 100))
                                                         ),
                                                         column(width=12,
                                                                column(width=6,colourpicker::colourInput("highcolor", "Color for up-regulation in KEGG map:", "red", showColour = "background")),
                                                                column(width=6,colourpicker::colourInput("lowcolor", "Color for down-regulation in KEGG map:", "blue", showColour = "background"))
                                                         )
                                        )
                                        )
                       ),
                       actionButton("start3","Start processing",icon=icon("play-circle"))
                       )
),
column(width=12,
       column(width=4, 
              id="inputarea",
              fileInput("ref_meta_file", "Select standard metabolite library ('.csv' or .txt')",
                        multiple = FALSE,
                        width="350px",
                        accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))
       ),
       conditionalPanel(condition = "input.step3",
                        column(width=4, 
                               id="inputarea",
                               fileInput("foldchange_file", "Select fold change file ('.csv' or .txt')",
                                         multiple = FALSE,
                                         width="350px",
                                         accept = c("text/csv","text/comma-separated-values,text/plain",".txt"))
                        )
       )
),
column(width=12,verbatimTextOutput("nText6")),
column(width=12,verbatimTextOutput("nText7")),
column(width=12,style="padding-left:10px;padding-right:10px;",
       conditionalPanel(condition = "!output.done",
                        column(width=12, style="text-align:left;",tags$a(target="_blank",href="metabolite_quantification_description.html","Introduction & Input and output file descriptions"))
       )
),
conditionalPanel(condition = "output.done",
                 column(width=12, style="text-align:center;",
                        id="inputarea",
                        downloadButton(style="background-color:#417ee0;color:#ffffff;","downloadQdata", label = "Download results")
                 )
),
style='primary'
),multiple=TRUE))
)

  
  


source("R/input.R")
source("R/data_preprocessing.R")
source("R/workflow.R")
source("R/network_analysis.R")
source("R/graphical_options.R")
source("R/load_packages.R")
# Define UI for data upload app ----

statistical_analysis_page_lite <- fluidPage(

  tags$div(
    id="maindiv1",
    column(width=12,
           mainPanel(
             width=12,
             bsCollapse(open = c("Input Files","Feature Selection"),multiple=TRUE,
                        bsCollapsePanel("Input Files",
                                        column(width=12,
                                               h4("Choose your Files:"),
                                               column(width=6, 
                                                      id="inputarea",
                                                      fileInput("featuretable", "Select feature table file ('.csv' or .txt', 100MB limit)",
                                                                multiple = FALSE,
                                                                width="350px",
                                                                accept = c("text/csv",
                                                                           "text/comma-separated-values,text/plain",
                                                                           ".txt"))
                                                      
                                               ),
                                               column(width=6, 
                                                      id="inputarea",
                                                      fileInput("classlabel", "Select class label file ('.csv' or '.txt', 100MB limit)",
                                                                multiple = FALSE,
                                                                width="350px",
                                                                accept = c("text/csv",
                                                                           "text/comma-separated-values,text/plain",
                                                                           ".txt"))
                                               ),
                                               
                                               column(width=12,style="margin-top:10px",
                                                      #column(width=6,style='padding-left:0;',numericInput(width="350px","pvalue_thresh", "P-value threshold:", 0.05, min = 0, max = 1)),
                                                      #column(width=6,numericInput(width="350px","foldchangethresh", "Fold change threshold (0-100 limit):", 0, min = 0, max = 100))
                                                      column(width=6,style='margin-top:25px;padding-left:0;',actionButton("argumentbutton1", "More arguments")),
                                                      #column(width=6,selectInput(width="350px","example_data","Use example data",c("FALSE","TRUE"),selected="FALSE"))
                                                      #radioButtons
                                                      column(width=6,radioButtons(width="350px","example_data","Use example data",inline=TRUE,c("FALSE","TRUE"),selected="FALSE"))
                                                      
                                               ),
                                               
                                               
                                               
                                               bsModal("argument_modal1", "More arguments for input files", "argumentbutton1", size = "large",
                                                       tags$div(
                                                         width=12,
                                                         style="height:275px",
                                                         tags$div(
                                                           column(width=6,style="margin-left:0;margin-right:0;",radioButtons("input_intensity_scale", "input intensity scale:", inline=TRUE,c("raw"= "raw","log2"="log2"),selected = "raw")),
                                                           bsTooltip("input_intensity_scale", "Are the intensities in the input feature table at raw scale or log2 scale?","bottom", options = list(container = "body")),
                                                           column(width=6,radioButtons("missing_val", "missing value notation:", inline=TRUE,c("0"= "0","NA"="NA"),selected = "0")),
                                                           bsTooltip("missing_val", "How are the missing values represented in the input data?","bottom", options = list(container = "body"))
                                                         ),
                                               column(width=6, 
                                                      id="inputarea",
                                                      textInput(width="350px","outloc", "Output folder name:","xmsPANDAout",placeholder="Default: xmsPANDAout")
                                               ))),
                                        ),
                                        style = "primary"),
                        bsCollapsePanel("Filtering",
                                        column(width=6,numericInput(width="350px","all_missing_thresh", "Minimum non-missing sample ratio:", 0.1, min = 0, max = 1)),
                                        bsTooltip("all_missing_thresh", "What propotion of total number of samples should have an intensity?","bottom"),
                                        column(width=6,numericInput(width="350px","rsd_filt_list", "Minimum overall variance:", 0, min = 0, max = 100)),
                                        bsTooltip("rsd_filt_list", "Features with less than minimum relative standard deviation across all samples are filtered","bottom"),
                                        column(width=6,numericInput(width="350px","group_missing_thresh", "Minimum non-missing sample ratio for group:", 0.8, min = 0, max = 1)),
                                        bsTooltip("group_missing_thresh", "Minimum propotion of samples in at least one group in which a non-missing signal value should be present","bottom"),
                                        style = "primary"),
                        
                        bsCollapsePanel("Imputation, transformation, and normalization",
                                        column(width=12,column(width=6,selectInput(width="350px","summary_na_replacement","Choose an imputation method:",c("halffeaturemin",
                                        "knn","randomforest","zeros","halfsamplemin","halfdatamin","none"))),
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
                                                                                                                        "median absolute deviation normalization"="mad_norm"),
                                                                   selected=c("log2transform")))),style = "primary"),
                        bsCollapsePanel("Feature Selection",column(12,textOutput("txt")),
                                        column(width=12, 
                                               id="inputarea",column(width=2,selectInput("analysismode","Select analysis mode:",c("classification","regression"))),
                                               column(width=2,selectInput("pairedanalysis", "Time-series design?", 
                                                                           c(True = "TRUE",False = "FALSE"),selected = "FALSE")), #inline=TRUE,
                                               column(width=8,
                                                      conditionalPanel(
                                                        condition = "input.analysismode == 'classification' && input.pairedanalysis == 'FALSE'",
                                                       # column(width=6,
                                                       # 
                                                              column(width=12,selectInput("featselmethodi", "Select feature selection methods
                                                                                           for classification \n(click in the box to edit or select more options)", 
                                                                                           multiple = TRUE,width="500px",c('limma (one-way ANOVA using LIMMA)'='limma', 'pls (partial least squares using VIP)'='pls',
                                                                                                             'limma2way (two-way ANOVA using LIMMA)'='limma2way',
                                                                                                             'limmarobust (one-way ANOVA using LIMMA robust)'='limmarobust',
                                                                                                             'limma2wayrobust (two-way ANOVA using LIMMA robust)'='limma2wayrobust',
                                                                                                             'lm1wayanova (one-way ANOVA using linear model)'='lm1wayanova',
                                                                                                             'lm2wayanova (two-way ANOVA using linear model)'='lm2wayanova',
                                                                                                             'lmreg (linear regression)'='lmreg',
                                                                                                             'lmregrobust (linear regression using robust sandwich estimator)'='lmregrobust',
                                                                                                             'logitreg (logistic regression)'='logitreg',
                                                                                                             'logitregrobust (logistic regression using robust sandwich estimator)'='logitregrobust',
                                                                                                             'poissonreg (poisson regression)'='poissonreg',
                                                                                                             'poissonregrobust (poisson regression using robust sandwich estimator)'='poissonregrobust',
                                                                                                             'rfesvm (recursive feature elimination SVM)'='rfesvm',
                                                                                                             'ttest (simple t-test)'='ttest',
                                                                                                             'wilcox (wilcoxon test)'='wilcox',
                                                                                                             'RF (random forest using Boruta algorithm)'='RF',
                                                                                                             'MARS (multiple adaptive regression splines)'='MARS',
                                                                                                             'spls (sparse partial least squares)'='spls',
                                                                                                             'spls2way (two-factor analysis using sparse partial least squares)'='spls2way',
                                                                                                             'o1pls (orthogonal partial least squares)'='o1pls',
                                                                                                             'pamr (nearest shrunked centroid method)'='pamr'),
                                                                                           selected=c('limma','pls')
                                                               ))),
                                                      conditionalPanel(
                                                        condition = "input.analysismode == 'classification' && input.pairedanalysis == 'TRUE'",
                                                        column(width=6,style='margin-top:25px;padding-left:0;',
                                                               #  div(style="display: inline-block;vertical-align:top",tags$p(style="margin:0;padding-top:4px;padding-bottom:4px;font-weight:bold","Choose feature selection methods:  ")),
                                                               # div(style="display: inline-block;vertical-align:top",actionButton("methodbuttonii", "View")),
                                                               # bsModal("method_modalii", "Method List", "methodbuttonii", size = "large",
                                                               column(width=6,selectInput("featselmethodii", width="500px","Select feature selection methods for 
                                                                                          classification with repeat measures (click in the box to edit or select more options)", multiple = TRUE,
                                                                                           c('limma1wayrepeat (one-way ANOVA repeated measures using LIMMA)'='limma1wayrepeat',
                                                                                             'limma2wayrepeat (two-way ANOVA repeated measures using LIMMA)'='limma2wayrepeat',
                                                                                             'limma1wayrepeatrobust (one-way ANOVA repeated measures using LIMMA robust)'='limma1wayrepeatrobust',
                                                                                             'limma2wayrepeatrobust (two-way ANOVA repeated measures using LIMMA robust)'='limma2wayrepeatrobust',
                                                                                             'lm1wayanovarepeat (one-way ANOVA repeated measures using linear model)'='lm1wayanovarepeat',
                                                                                             'lm2wayanovarepeat (two-way ANOVA repeated measures using linear model)'='lm2wayanovarepeat',
                                                                                             'spls1wayrepeat (one-way ANOVA repeated measures using sparse partial least squares)'='spls1wayrepeat',
                                                                                             'spls2wayrepeat (two-way ANOVA repeated measures using sparse partial least squares)'='spls2wayrepeat',
                                                                                             'ttestrepeat (simple t-test repeated measures)'='ttestrepeat',
                                                                                             'wilcoxrepeat (wilcoxon test repeated measures)'='wilcoxrepeat'),
                                                                                           selected=c('lm1wayanovarepeat'))))),
                                                      conditionalPanel(
                                                        condition = "input.analysismode == 'regression' && input.pairedanalysis == 'FALSE'",
                                                        column(width=6,style='margin-top:25px;padding-left:0;',
                                                               
                                                               column(width=6,selectInput("featselmethodiii", width="500px","Select feature selection methods for regression", multiple = TRUE,
                                                                                           c('lmreg (linear regression)'='lmreg',
                                                                                             'lmregrobust (linear regression using robust sandwich estimator)'='lmregrobust',
                                                                                             'poissonreg (poisson regression)'='poissonreg',
                                                                                             'poissonregrobust (poisson regression using robust sandwich estimator)'='poissonregrobust',
                                                                                             'RF (random forest using boruta algorithm)'='RF',
                                                                                             'MARS (multiple adaptive regression splines)'='MARS',
                                                                                             'pls (partial least squares using VIP)'='pls',
                                                                                             'spls (sparse partial least squares)'='spls',
                                                                                             'o1pls (orthogonal partial least squares)'='o1pls'),
                                                                                           
                                                                                           selected=c('lmreg'))
                                                               ) 
                                                               
                                                        )
                                                      ),
                                                      conditionalPanel(
                                                        condition = "input.analysismode == 'regression' && input.pairedanalysis == 'TRUE'",
                                                        column(width=6,style='margin-top:25px;padding-left:0;',
                                                               
                                                               column(width=6,selectInput("featselmethodiv", width="500px","Select feature selection methods for regression", multiple = TRUE,
                                                                                           c('lmregrepeat (linear regression with mixed effects model)'='lmregrepeat',
                                                                                             
                                                                                             'plsrepeat (multilevel partial least squares using VIP)'='plsrepeat',
                                                                                             'splsrepeat (multilevel sparse partial least squares)'='splsrepeat'),
                                                                                           
                                                                                           
                                                                                           selected=c('lmregrepeat'))
                                                               )
                                                               
                                                        )
                                                      )
                                               ),
                                               column(width=12,style="margin-top:10px",
                                                      #column(width=6,style='padding-left:0;',numericInput(width="350px","pvalue_thresh", "P-value threshold:", 0.05, min = 0, max = 1)),
                                                      #column(width=6,numericInput(width="350px","foldchangethresh", "Fold change threshold (0-100 limit):", 0, min = 0, max = 100))
                                                      column(width=12,style='margin-top:25px;padding-left:0;',actionButton("argumentbutton2", "More arguments"))
                                               ),
                                               
                                               
                                               bsModal("argument_modal2", "More arguments for feature selection", "argumentbutton2", size = "large",
                                                       tags$div(
                                                         width=12,
                                                         style="height:575px",
                                                         
                                                         column(width=6,style='padding-left:0;',numericInput(width="350px","pvalue_thresh", "P-value threshold:", 0.05, min = 0, max = 1)),
                                                         column(width=6,numericInput(width="350px","foldchangethresh", "Fold change threshold (0-100 limit):", 0, min = 0, max = 100)),
                                                         column(width=6,style='padding-left:0;',numericInput(width="350px","fdrthresh", "False discovery threshold:", 0.05, min = 0, max = 1)),
                                                         column(width=6,selectInput(width="350px","fdr_method","Choose FDR correction method:",
                                                                                    c("BH (Benjamini-Hochberg 1995)"="BH","ST (Storey & Tibshirani 2001)"="ST","Strimmer (Strimmer 2008)"="Strimmer",
                                                                                      "BY (Benjamini-Yekutieli 2001)"="BY", "Bonferroni"="bonferroni","none"="none"))),
                                                         column(width=6,numericInput(width="350px","pls_vip_thresh", "VIP threshold (1-100 limit):", 2, min = 1, max = 100)),
                                                         #column(width=6,style='margin-top:25px;padding-left:0;',actionButton("argumentbutton", "More arguments")),
                                                         
                                                        # column(width=6,numericInput(width="350px","kfold", "k for k-fold Cross Validation (1-10000 limit):", 10, min = 1, max = 10000)),
                                                         
                                                         #  bsTooltip("pls_vip_thresh", "This argument only works for PLS or O1PLS models","bottom", options = list(container = "body")),
                                                         #   column(width=6,numericInput(width="350px","pls_ncomp", "Max number of components to consider:", 5, min = 1, max = 100000)),
                                                         
                                                         bsTooltip("pls_ncomp", "This argument only works for PLS, sPLS, or O1PLS models","bottom", options = list(container = "body")),
                                                         column(width=6,numericInput(width="350px","max_comp_sel", "Max number of components to use for feature selection:", 2, 
                                                                                     min = 1, max = 100000)),
                                                         
                                                         column(width=6,style='padding-left:0;',selectInput(width="350px","lme.modeltype","Model type for mixed effects model:",
                                                                                                            c("Random Intercept"="RI","Random Intercept & Random Slope"="RIRS"),selected="Random Intercept")),
                                                         column(width=12,style="margin-top:10px",
                                                                column(width=6,radioButtons("limmadecideTests", "Should the LIMMA decide tests [-1 (down),0 (no change),1 (up)] be performed?", inline=TRUE,c(True = "TRUE",False = "FALSE"),
                                                                                            selected = "FALSE")),
                                                                
                                                                column(width=6,selectInput(width="350px","limma.contrasts.type","Select LIMMA contrasts type:",c("Sum contrats as in ANOVA (contr.sum)"="contr.sum",
                                                                                                                                                                 "Use the first group as the reference (contr.treatment)"="contr.treatment"),
                                                                                           selected="contr.sum"))
                                                                
                                                         )
                                                         
                                                       )
                                                       
                                               )
                                        
                                        ),
                                        style = "primary"),
                        bsCollapsePanel("Network analysis",
                                        
                                        column(12,radioButtons("globalcor", "Perform correlation-based network analysis of selected features:", inline=TRUE,c(True = "TRUE",False = "FALSE"),selected = "FALSE")),
                                        conditionalPanel(
                                          condition = "input.globalcor == 'TRUE'",
                                          column(width=12, 
                                                 id="inputarea",
                                                 column(width=6,numericInput(width="350px","abs_cor_thresh", "Absolute correlation threshold:", 0.4, min = 0, max = 1)),
                                                 column(width=6,numericInput(width="350px","cor_fdrthresh", "FDR threshold for correlation analysis:", 0.05, min = 0, max = 1)),
                                                 column(width=6,selectInput(width="350px","cor_method","Correlation method:",c("spearman","pearson"))),
                                                 column(width=6,selectInput(width="350px","net_legend","Network legend:",c("TRUE","FALSE")))
                                          )
                                        ),
                                        style="primary"),
                        bsCollapsePanel("Other options (HCA, PCA, classification accuracy, and graphics)",
                        column(width=12,
                              # column(width=6,selectInput(width="350px","heatmap_color_scheme","Color theme for up or down status in heatmap/Manhattan plot/Volcano plot:",c("redblue","orangeblue","bluered","blueorange"),selected="redblue")),
                               #column(width=6,selectInput(width="350px","output_format","Output format (PDF or PNG for high-res):",c("png","pdf"))),
                               #column(width=6,numericInput(width="350px","pca_cex_val", "Size of points on PCA plots (1-20 limit):", 4, min = 1, max = 20))
                              column(width=6,
                                     id="inputarea",
                                     textInput(width="350px","ylabel.text", "Label for y-axis in boxplots, barplots, and lineplots:","Abundance",placeholder="Default: Abundance")),
                              
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
                               #   column(width=6,numericInput(width="350px","ellipse_conf_level", "Confidence interval for PCA ellipses (0-1 limit):", 0.95, min = 0, max = 1)),
                             #  column(width=6,selectInput(width="350px","pca_ellipse","Should ellipse be plotted on PCA plots?",c("TRUE","FALSE"),selected="TRUE")),
                               column(width=6,selectInput(width="350px","boxplot.type", "Boxplot type:", c("ggplot","simple"),selected="simple")),
                               column(width=6,selectInput(width="350px","boxplot.bool", "Generate boxplots:", c("TRUE","FALSE")))
                               # column(width=6,numericInput(width="350px","pca_cex_val", "Size of points on PCA plots (1-20 limit):", 4, min = 1, max = 20))
                        ),
                        column(width=12,
                              
                               column(width=6,numericInput(width="350px","hca.cex.legend", "Size of HCA legend text (set to -1 to turn off):", 0.7, min = 0.001, max = 20)),
                               column(width=6,numericInput(width="350px","cex.plots", "Size of x and y-axis text in the plots (used for all plots):", 0.6, min = 0.001))
                        ),
                        
                        column(width=12,
                               
                               column(width=6,
                                      id="inputarea",
                                      selectInput(width="350px","timeseries.lineplots", "Plot time series lineplots (for time-series data):",c("FALSE","TRUE"))
                                      
                               ),
                               column(width=6,
                                      id="inputarea",
                                      selectInput(width="350px","alphabetical.order", "Plot classes on the x-axis in alphabetical order:",
                                                  c("Automatically sort by alphabetical order"="TRUE",
                                               "Use the order in the class labels file"="FALSE"),selected="FALSE")
                                      
                               )),
                       
                        
                        column(width=12,
                               
                            
                               column(width=6,
                                      id="inputarea",
                                      selectInput(width="350px","pca.global.eval", "Perform PCA:",c("FALSE","TRUE"))
                                      
                               ),
                               
                               column(width=6,
                                      id="inputarea",
                                      selectInput(width="350px","evaluate.classification.accuracy", "Evaluate k-fold classification accuracy and ROC:",c("FALSE","TRUE"))
                                      
                               ),
                               
                               column(width=6,
                                      id="inputarea",
                                      selectInput(width="350px","hca_type", "Two-way or one-way hierarchicial clustering:",c("two-way","one-way"))
                                      
                               ),
                               column(width=6,
                                      id="inputarea",
                                      selectInput(width="350px","hca.labRow.value", "Show row and column names in HCA heatmaps:",c("TRUE","FALSE"))
                                      
                               )
                               
                        ),
                       
                        style="primary")
                        
             )))),
  

  column(12,style='padding-left:0px;margin-bottom:20px;',
    mainPanel(
      style='padding-left:0px;',
    verbatimTextOutput("nText2"),
    verbatimTextOutput("nText"),
    bsAlert("alert"),
    
    actionButton("go","Start processing",icon=icon("play-circle")),
    downloadButton("downloadData", label = "Download results"),
actionButton("resetAll", "Refresh app",icon=icon("refresh"))
  )),

  #column(12,style='padding-top:5px;padding-left:0;',tags$div(h4("Output"))),
  uiOutput("output_results")
  
)
  



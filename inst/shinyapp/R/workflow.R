source("R/load_packages.R")

workflow<-fluidRow(
  tags$div(
    id="maindiv2",
column(12,textOutput("txt")),
   
      column(width=12, 
             id="inputarea",
             column(width=6,selectInput(width="350px","analysismode","Select analysis mode:",c("classification","regression"))),
             column(width=6,radioButtons("pairedanalysis", "Is this a repeated-measurement design?", inline=TRUE,c(True = "TRUE",False = "FALSE"),selected = "FALSE")),
             column(width=12,
               conditionalPanel(
                 condition = "input.analysismode == 'classification' && input.pairedanalysis == 'FALSE'",
                 column(width=12,
# div(style="display: inline-block;vertical-align:top",tags$p(style="margin:0;padding-top:4px;padding-bottom:4px;font-weight:bold","Choose feature selection methods:  ")),
#div(style="display: inline-block;vertical-align:top",actionButton("methodbuttoni", "View")),
#bsModal("method_modali", "Method List", "methodbuttoni", size = "large",
                                column(width=12,selectInput(width="800px","featselmethodi", "Select feature selection methods for classification (click in the box to edit or select more options)", multiple = TRUE,c('limma (one-way ANOVA using LIMMA)'='limma',
                                                                     'pls (partial least squares using VIP)'='pls',
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
                                                                     'pamr (microarrays algorithm based on the nearest shrunked centroid method)'='pamr'),
                                                   selected=c('limma','pls')
                                                   )
                                )
#)
                 )
               ),
               conditionalPanel(
                 condition = "input.analysismode == 'classification' && input.pairedanalysis == 'TRUE'",
                 column(width=12,
#  div(style="display: inline-block;vertical-align:top",tags$p(style="margin:0;padding-top:4px;padding-bottom:4px;font-weight:bold","Choose feature selection methods:  ")),
# div(style="display: inline-block;vertical-align:top",actionButton("methodbuttonii", "View")),
# bsModal("method_modalii", "Method List", "methodbuttonii", size = "large",
                                column(width=12,selectInput(width="800px","featselmethodii", "Select feature selection methods for classification with repeat measures (click in the box to edit or select more options)", multiple = TRUE,
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
                                                   selected=c('lm1wayanovarepeat'))
                                                   
                                )
# )
                 )
               ),
               conditionalPanel(
                 condition = "input.analysismode == 'regression' && input.pairedanalysis == 'FALSE'",
                 column(width=12,

                                column(width=12,selectInput(width="800px","featselmethodiii", "Select feature selection methods for regression", multiple = TRUE,
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
                column(width=12,

                                               column(width=12,selectInput(width="800px","featselmethodiv", "Select feature selection methods for regression", multiple = TRUE,
                                                                  c('lmregrepeat (linear regression with mixed effects model)'='lmregrepeat',
                                                                                   
                                                                                    'plsrepeat (multilevel partial least squares using VIP)'='plsrepeat',
                                                                                    'splsrepeat (multilevel sparse partial least squares)'='splsrepeat'),
                                                                                    
                                                                 
                                                                  selected=c('lmregrepeat'))
                                               )

                                )
               )
             ),
             
             column(width=12,style="margin-top:10px",
               column(width=6,style='padding-left:0;',numericInput(width="350px","pvalue_thresh", "P-value threshold:", 0.05, min = 0, max = 1)),
                column(width=6,numericInput(width="350px","foldchangethresh", "Fold change threshold (0-100 limit):", 0, min = 0, max = 100))
               
),
column(width=12,style="margin-top:10px",
column(width=6,style='padding-left:0;',numericInput(width="350px","fdrthresh", "False discovery threshold:", 0.05, min = 0, max = 1)),
column(width=6,selectInput(width="350px","fdr_method","Choose FDR correction method:",
                           c("BH (Benjamini-Hochberg 1995)"="BH","ST (Storey & Tibshirani 2001)"="ST","Strimmer (Strimmer 2008)"="Strimmer",
                             "BY (Benjamini-Yekutieli 2001)"="BY", "Bonferroni"="bonferroni","none"="none")))
),
column(width=12,style="margin-top:10px",
column(width=6,numericInput(width="350px","pls_vip_thresh", "VIP threshold (1-100 limit):", 2, min = 1, max = 100)),
column(width=6,style='margin-top:25px;padding-left:0;',actionButton("argumentbutton", "More arguments"))
),
               bsModal("argument_modal", "More arguments for feature selection", "argumentbutton", size = "large",
                       tags$div(
                         width=12,
                         style="height:575px",
                         
                         column(width=6,numericInput(width="350px","kfold", "k for k-fold Cross Validation (1-10000 limit):", 10, min = 1, max = 10000)),
                         
                         bsTooltip("pls_vip_thresh", "This argument only works for PLS or O1PLS models","bottom", options = list(container = "body")),
                         column(width=6,numericInput(width="350px","pls_ncomp", "Max number of components to consider:", 5, min = 1, max = 100000)),
                         bsTooltip("pls_ncomp", "This argument only works for PLS, sPLS, or O1PLS models","bottom", options = list(container = "body")),
                         column(width=6,numericInput(width="350px","max_comp_sel", "Number of components to use for VIP selection:", 1, min = 1, max = 100000)),
                         bsTooltip("max_comp_sel", "This argument only works for PLS, sPLS, or O1PLS models","bottom", options = list(container = "body")),
                         column(width=6,selectInput(width="350px","optselect","Find optimal number of components:",c("TRUE","FALSE"))),
                         bsTooltip("optselect", "This argument only works for PLS, sPLS, or O1PLS models","bottom", options = list(container = "body")),
                         column(width=6,selectInput(width="350px","pls.vip.selection","VIP summarization across multiple components:",c("max","mean"))),
                         bsTooltip("optselect", "This argument only works for PLS, sPLS, or O1PLS models","bottom", options = list(container = "body")),
                         column(width=6,
                                tags$label("Number of permutations for calculating p-values:", `for` = "permu_switch"),
                                div(style="display: inline-block;vertical-align:top; width: 100px;", switchInput(inputId = "permu_switch",value = FALSE)),
                                div(style="display: inline-block;vertical-align:top; width: 250px;", numericInput(inputId = "pls_permut_count", label = NULL, value = 1000, min = 1, max = 100000))
                         ),
                         bsTooltip("pls_permut_count", "This argument only works for PLS, sPLS, or O1PLS models","bottom", options = list(container = "body")),
                         column(width=6,conditionalPanel(
                          condition = "output.checkmax_varsel",
                          numericInput(width="350px","max_varsel", "Max number of variables to be used:", 100, min = 1, max = 100000),
                          bsTooltip("max_varsel", "This argument only works for sPLS, spls1wayrepeat, spls2wayrepeat, rfesvm, and Random Forest","bottom", options = list(container = "body"))
                         )),

column(width=12,style="margin-top:10px",
column(width=6,radioButtons("limmadecideTests", "Should the LIMMA decide tests [-1 (down),0 (no change),1 (up)] be performed?", inline=TRUE,c(True = "TRUE",False = "FALSE"),
                            selected = "FALSE")),

column(width=6,selectInput(width="350px","limma.contrasts.type","Select contrast type:",c("Sum contrats as in ANOVA (contr.sum)"="contr.sum",
                                                                                          "Use the first group as the reference (contr.treatment)"="contr.treatment"),
                           selected="contr.sum"))
       
),
column(width=12,style="margin-top:10px",
     #  column(width=6,radioButtons("limmadecideTests", "Should the LIMMA decide tests [-1 (down),0 (no change),1 (up)] be performed?", inline=TRUE,c(True = "TRUE",False = "FALSE"),selected = "TRUE")),
       
     column(width=6,selectInput(width="350px","similarity.matrix","Select similarity matrix for hierarchical clustering:",c("correlation (similarity.matrix)"="correlation",
                                                                                                                            "topological overlap (TOM)"="TOM"),
                                selected="correlation (similarity.matrix)")),
     column(width=6,#style='padding-left:10;',
              tags$label("Network based feature ranking (differential network analysis):", `for` = "netbasedfeatranking_switch"),
              div(style="width: 100px;", switchInput(inputId = "netbasedfeatranking_switch",value = FALSE))
              
       )
),

                         
column(width=12,style="margin-top:10px",
column(width=6,
#conditionalPanel(
#                          condition = "output.checkaggregationmethod",
                           selectInput(width="350px","aggregation_method","Method for aggregating results:",
                                       c("None"="none","RankAggreg (cross entropy)"="RankAggreg","RankAggregGA (genetic algorithm)"="RankAggregGA","Consensus"="consensus"),selected="None")
                            ),
#column(width=6,radioButtons("limmadecideTests", "Should the LIMMA decide tests be performed?", inline=TRUE,c(True = "TRUE",False = "FALSE"),selected = "TRUE")),
column(width=6,radioButtons("balance.classes", "Should the ROSE method be used for balancing classes?", inline=TRUE,c(True = "TRUE",False = "FALSE"),selected = "FALSE"))
),

column(width=12,style="margin-top:10px",
          column(width=6,
                 #conditionalPanel(
                 #                          condition = "output.checkaggregationmethod",
                 selectInput(width="350px","outlier.method","Outlier detection method:",
                             c("Tukey's 1.5*interquartile range rule on summed values"="sumtukey","Tukey's 1.5*interquartile range rule on PCA scores"="pcatukey",
                               "Using the multivariate outliers detection method (pcout)"="pcout","Using robust Mahalanobis distance with a Chi-square test (pcachisq)"="pcachisq"),
                             selected="pcout")),
          column(width=6,
                 #conditionalPanel(
                 #                          condition = "output.checkaggregationmethod",
                 selectInput(width="350px","vcovHC.type","Estimation type for the robust covariance matrix estimation:",
                             c("HC0"="HC0","HC3"="HC3"),selected="HC3")
          )
),
column(width=12,style="margin-top:10px",
          column(width=6,style='padding-left:0;',selectInput(width="350px","lme.modeltype","Model type for mixed effects model:",
                                                             c("Random Intercept"="RI","Random Intercept & Random Slope"="RIRS"),selected="Random Intercept")),
          column(width=6,numericInput(width="350px","log2.transform.constant", "Small constant for log2 transformation:", 1, min = 0))

)




                       )
               
             )

)
      
    

))
  

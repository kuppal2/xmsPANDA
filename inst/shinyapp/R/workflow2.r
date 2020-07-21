library(shiny)
library(shinyWidgets)


workflow2<-fluidRow(
  tags$div(
    id="maindiv",
    #column(12,radioButtons("workflow", "Choose a statistical analysis workflow:", inline=TRUE,c('Workflow I (traditional approach)' = "workflowI",'Workflow II (advanced approach(under development))' = "workflowII"),selected = "workflowI")),
    #column(12,textOutput("txt")),
    #conditionalPanel(
      #condition = "input.workflow == 'workflowI'",
      column(width=12, 
             id="inputarea",
           
                 column(width=12,
                        div(style="display: inline-block;vertical-align:top",tags$p(style="margin:0;padding-top:4px;padding-bottom:4px;font-weight:bold","Choose feature selection methods:  ")),
                        div(style="display: inline-block;vertical-align:top",actionButton("methodbuttoni2", "View")),
                        bsModal("method_modali2", "Method List", "methodbuttoni2", size = "large",
                                checkboxGroupInput("featselmethodi2", "", inline = FALSE,
                                                   choiceNames =list('limma (one-way ANOVA using LIMMA)',
                                                                     'pls (partial least squares)',
                                                                     'lmreg (linear regression)',
                                                                     'lmregrobust (linear regression using robust sandwich estimator)',
                                                                     'logitreg (logistic regression)',
                                                                     'logitregrobust (logistic regression using robust sandwich estimator)',
                                                                      
                                                                     'rferadial (recursive feature elimination SVM with RBF kernel)',
                                                                        'rfe (recursive feature elimination SVM with linear kernel)',
                                                                     't.test (simple t-test)',
                                                                     'wilcox.test (wilcoxon test)',
                                                                     'welch.test (welch test)',
                                                                      'f.test (f test)',
                                                                        'kruskal.test (kruskal-wallis test)',
                                                                     'rf (random forest using the randomForest package)',
                                                                     'rfboruta (random forest using Boruta algorithm)',
'lasso',
'elasticnet',
                                                                     'spls (sparse partial least squares)',
                                                                     'o1pls (orthogonal partial least squares)'),
                                                                    
                                                   choiceValues =list('limma','pls','lmreg', 'lmregrobust',
                                                                      'logitreg','logitregrobust','rferadial','rfe','t.test','wilcox.test','welch.test','f.test','kruskal.test','rf','rfboruta','lasso','elasticnet','spls','o1pls'),
                                                   selected=c('limma','rfe','rf','pls','lasso'),
                                                   width='800px'
                                )
                        )
                 )
               ),
             
             column(width=12,style="margin-top:10px",
                column(width=6,style='padding-left:0;',numericInput(width="350px","num.var.sel", "Max number of variables to select from each method:", 10, min = 3, max = 1000)),
               
                column(width=6,style='padding-left:0;',selectInput(width="350px","split.train.test", "Randmly split the input data into training and test sets:", c("TRUE","FALSE"))),


                column(width=6,style='padding-left:0;',numericInput(width="350px","train.pct", "Proportion of samples to use for training:",  0.7, min = 0, max = 1)),

                column(width=6,style='padding-left:0;',selectInput(width="350px","learningsetmethod", "Method to generate learning sets:", choices=c("Cross-validation"="CV","Monte-carlo cross-validation"="MCCV","Bootstrap"="bootstrap"))),
                column(width=6,style='padding-left:0;',selectInput(width="350px","balance.classes", "Balance classes by generating synthetic samples:", choices=c("FALSE","TRUE"))),
                column(width=6,style='padding-left:0;',selectInput(width="350px","aggregation_method","Method for aggregating results:",choices=c("Probability-based"="consensus","RankAggreg (cross entropy)"="RankAggreg","RankAggregGA (genetic algorithm)"="RankAggregGA"))),

#conditionalPanel(
#condition = "input.num.var.sel == '10'",
#condition = "input.aggregation_method == 'consensus'",
#sliderInput("prop.select.thresh", "Min. proportion of subsets in which a variable is selected:", min=0, max=1, value=0.7)
#)
column(width=6,numericInput(width="350px","kfold", "k for k-fold Cross Validation (1-10000 limit):", 10, min = 1, max = 10000)),
column(width=6,style='padding-left:0;',numericInput(width="350px","prop.select.thresh", "Min. proportion of subsets in which a variable is selected (for probability-based aggregation):", 0.7, min = 0, max = 1))


                       )
               )
             )



help_page <- fluidPage(
  
  column(12,tags$h4("User Manual:")),
  column(12,tags$p(style="font-size: 15px;","Click ",tags$a(href='xmsPANDA-manual.pdf',target="_blank","here")," to see the user manual.")),
  column(12,tags$p(style="font-size: 15px;","Click ",tags$a(href='https://joneslabemory.github.io/March182019-metabolomics-analysis-workshop/',target="_blank","here")," to see the R code instruction.")),
  column(12,tags$h4("Input File Format:")),
  column(12,tags$p(style="margin-top:10px;font-weight:bold; font-size: 17px;","feature table file")),
  column(12,tags$p(style="text-align:justify;font-size: 15px;","For feature table file, two formats are accepted. Format 1: the first 2 columns are m/z and time; Format 2: the first column is labeled as Name and includes named metabolites and mz_time for unnamed features. The remaining columns should correspond to the samples in the class labels file and each column is the intensity profile of a sample. Here is an example for format 1 with m/z and time:",icon("hand-point-down", lib = "font-awesome", "fa-2x"))),
  column(12, div(style="display:block;margin-top:10px; margin-left: auto;margin-right: auto;width: 50%;",tableOutput('example_feat'))),

 column(12,tags$p(style="text-align:justify;font-size: 15px;","Here is an example for feature table format 2 with named and unnamed features:",icon("hand-point-down", lib = "font-awesome", "fa-2x"))),
 column(12, div(style="display:block;margin-top:10px; margin-left: auto;margin-right: auto;width: 50%;",tableOutput('example_feat2'))
         #div(style="display: inline-block;vertical-align:top",tableOutput('example_feat'))
         #div(style=";margin-left:40px; width:40%;display: inline-block;vertical-align:top",tags$p(style="text-align:justify;font-size: 15px;","The feature table should include m/z, retention time, and measured intensity in each sample for each analyte. The first 2 columns should be the m/z and time. The remaining columns should correspond to the samples in the class labels file with each column including the intensity profile of a sample"))
         ),
  column(12,div(style="margin-top:23px;margin-right:10px;display: inline-block;vertical-align:top",tags$p(style="font-weight:bold; font-size: 17px;","class label file"))
         #div(style="display: inline-block;vertical-align:top",selectInput(width="310px","classlabel_option","",c("multiclass comparison","multiclass comparison with covariates","regression","two-way anova","one factor repeatedmeasures","two factor repeatedmeasures")))
         ),
  column(12,
         column(4,align='center',style="display:block; margin-left: auto;margin-right:auto;",
                tableOutput('example_classlabel_table')
         ),
         column(8,
                div(style='margin-bottom:30px;',uiOutput("example_classlabel_text"),
                    tags$p(icon("hand-point-left", lib = "font-awesome", "fa-2x"),'the example format is on the left side.')
                ),
                div(
                  div(style='margin-bottom:10px;',tags$label("Class Label Options", `for` = "classlabel_option" )),
                  awesomeRadio(inputId = "classlabel_option",label = NULL, choices = c("multiclass comparison", "multiclass comparison with covariates", "regression", "two-way anova", "one factor repeatedmeasures", "two factor repeatedmeasures"), selected = "multiclass comparison", status = "primary")
                ) 
        )
            # conditionalPanel(
            #   condition = "input.classlabel_option == 'multiclass comparison'",
            #   column(12,
            #          div(style="display: inline-block;vertical-align:top",tableOutput('example_multiclass_comparison')),
            #          div(style=";margin-left:40px; width:75%;display: inline-block;vertical-align:top",tags$p(style="text-align:justify;font-size: 15px;","To do multiclass comparison, the class labels file should include two columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the factor1 (or class1) which contains the group information."))
            #   )
            # ),
            # conditionalPanel(
            #   condition = "input.classlabel_option == 'multiclass comparison with covariates'",
            #   column(12,
            #          div(style="display: inline-block;vertical-align:top",tableOutput('example_multiclass_comparison_covariates')),
            #          div(style=";margin-left:40px; width:75%;display: inline-block;vertical-align:top",tags$p(style="text-align:justify;font-size: 15px;","To do multiclass comparison adjusted for some covariates, the class labels file should include sampleID, group information and covariates. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the class which contains the group information. The remaing columns should be covariates information."))
            #   )
            # ),
            # conditionalPanel(
            #   condition = "input.classlabel_option == 'regression'",
            #   column(12,
            #          div(style="display: inline-block;vertical-align:top",tableOutput('example_regression')),
            #          div(style=";margin-left:40px; width:75%;display: inline-block;vertical-align:top",tags$p(style="text-align:justify;font-size: 15px;","To do the regression analysis (adjusted for covariates), the class labels file should include at least two columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the response variable which should be numeric variable. If you want to adjust for some covariates, put the covariates information after the first 2 columns."))
            #   )
            # ),
            # conditionalPanel(
            #   condition = "input.classlabel_option == 'two-way anova'",
            #   column(12,
            #          div(style="display: inline-block;vertical-align:top",tableOutput('example_two_way_anova')),
            #          div(style=";margin-left:40px; width:75%;display: inline-block;vertical-align:top",tags$p(style="text-align:justify;font-size: 15px;","To do the two-way ANOVA analysis, the class labels file should include three columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the factor1 information and the third column is the factor2 information."))
            #   )
            # ),
            # conditionalPanel(
            #   condition = "input.classlabel_option == 'one factor repeatedmeasures'",
            #   column(12,
            #          div(style="display: inline-block;vertical-align:top",tableOutput('example_one_factor_repeatedmeasures')),
            #          div(style=";margin-left:40px; width:75%;display: inline-block;vertical-align:top",tags$p(style="text-align:justify;font-size: 15px;","To do the one-way ANOVA analysis with repeated-measurement data, the class labels file should include three columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the subject ID and the third column is the factor1 which contains the group information."))
            #   )
            # ),
            # conditionalPanel(
            #   condition = "input.classlabel_option == 'two factor repeatedmeasures'",
            #   column(12,
            #          div(style="display: inline-block;vertical-align:top",tableOutput('example_two_factor_repeatedmeasures')),
            #          div(style=";margin-left:40px; width:75%;display: inline-block;vertical-align:top",tags$p(style="text-align:justify;font-size: 15px;","To do the two-way ANOVA analysis with repeated-measurement data, the class labels file should include four columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the subject ID. The third and the fourth column are the factor1 and factor2 information."))
            #   )
            # )
  ),
  column(12,tags$h4("Download Resources:")),
  column(12,tags$p(style="font-size: 15px;","Download ",tags$a(href='https://github.com/kuppal2/xmsPANDA',target="_blank","xmsPANDA")," from Github")),
  column(12,tags$p(style="font-size: 15px;","Download ",tags$a(href='https://github.com/kuppal2/xmsPANDA/tree/master/examples_and_manual/Example_feature_table_and_classlabels',target="_blank","example data")," from Github."))
  
)

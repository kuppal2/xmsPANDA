options(shiny.maxRequestSize=100*1024^2)
options(shiny.sanitize.errors=FALSE)
library('xmsPANDA')
#library('lsmeans')
#library('car')
#library('KEGGREST')
#.libPaths("R/source_codes/library/")
library(shiny)
library(shinyjs)
library(DT)
source("R/source_codes/xmsPANDA_v1.0.8.47.R")
#source("R/source_codes/xMSquant_v0.0.3.R")



# Server logic
server <- function(input, output, session) {
  
  ##################################  Introduction Page #################################################  
  
  ##################################  Main Analysis Page #################################################
  done <- reactiveValues(count = 0)
  go <- reactiveValues(count = 0)
  id1 <- NULL
  check <- reactiveValues(count = 0)
  output$checkpvalue <- reactive({
    if(input$analysismode == 'classification' && input$pairedanalysis == 'FALSE'){
      sum(input$featselmethodi%in%c('limma','limma2way','lm1wayanova','lm2wayanova','lmreg','logitreg','rfesvm','ttest','wilcox','RF','MARS','pamr'))>0
    }else{
      if(input$analysismode == 'classification' && input$pairedanalysis == 'TRUE'){
        sum(input$featselmethodii%in%c('limma1wayrepeat','limma2wayrepeat','lm1wayanovarepeat','lm2wayanovarepeat','ttestrepeat','wilcoxrepeat'))>0
      }else{
        if(input$analysismode == 'regression' && input$pairedanalysis == 'FALSE'){
          sum(input$featselmethodiii%in%c('lmreg','RF','MARS'))>0
        }else{
            if(input$analysismode == 'regression' && input$pairedanalysis == 'TRUE'){
                     sum(input$featselmethodiv%in%c('lmregrepeat','plsrepeat','splsrepeat'))>0
                   }
            
        }
        
        
      }}})
  outputOptions(output, "checkpvalue", suspendWhenHidden = FALSE)
  output$checkvip <- reactive({
    if(input$analysismode == 'classification' && input$pairedanalysis == 'FALSE'){
      sum(input$featselmethodi%in%c('pls','o1pls','spls'))>0
    }else{
      if(input$analysismode == 'classification' && input$pairedanalysis == 'TRUE'){
        sum(input$featselmethodi%in%c('spls1wayrepeat','spls2wayrepeat'))>0
      }else{
        if(input$analysismode == 'regression' && input$pairedanalysis == 'FALSE'){
          sum(input$featselmethodiii%in%c('pls','o1pls','spls'))>0
        }else{
            
            if(input$analysismode == 'regression' && input$pairedanalysis == 'TRUE'){
              sum(input$featselmethodiv%in%c('lmregrepeat','plsrepeat','splsrepeat'))>0
            }
        }
        
        
      }}})
  outputOptions(output, "checkvip", suspendWhenHidden = FALSE)
  output$checkmax_varsel <- reactive({
    if(input$analysismode == 'classification' && input$pairedanalysis == 'FALSE'){
      sum(input$featselmethodi%in%c('rfesvm','RF','spls'))>0
    }else{
      if(input$analysismode == 'classification' && input$pairedanalysis == 'TRUE'){
        sum(input$featselmethodi%in%c('spls1wayrepeat','spls2wayrepeat'))>0
      }else{
        if(input$analysismode == 'regression' && input$pairedanalysis == 'FALSE'){
          sum(input$featselmethodiii%in%c('RF','spls'))>0
        }}}})
  outputOptions(output, "checkmax_varsel", suspendWhenHidden = FALSE)
  output$checkaggregationmethod <- reactive({
    if(input$analysismode == 'classification' && input$pairedanalysis == 'FALSE'){
      length(input$featselmethodi)>1
    }else{
      if(input$analysismode == 'classification' && input$pairedanalysis == 'TRUE'){
        length(input$featselmethodii)>1
      }else{
        if(input$analysismode == 'regression' && input$pairedanalysis == 'FALSE'){
          length(input$featselmethodiii)>1
        }}}})
  outputOptions(output, "checkaggregationmethod", suspendWhenHidden = FALSE)
  #output$checkplspermutation <- reactive({
  #  if(input$analysismode == 'classification' && input$pairedanalysis == 'FALSE'){
  #    sum(input$featselmethodi%in%c('pls','spls','o1pls'))>0
  #  }else{
  #    if(input$analysismode == 'classification' && input$pairedanalysis == 'TRUE'){
  #      FALSE
  #    }else{
  #      if(input$analysismode == 'regression' && input$pairedanalysis == 'FALSE'){
  #        sum(input$featselmethodiii%in%c('pls','spls','o1pls'))>0
  #      }}}})
  #outputOptions(output, "checkplspermutation", suspendWhenHidden = FALSE)

  observe({
    if (input$permu_switch==TRUE) {
      shinyjs::enable("pls_permut_count")
    } else {
      shinyjs::disable("pls_permut_count")
    }
  })
  
  all_alert <- reactive({

    if(!is.integer(input$numreplicate)) {
      closeAlert(session, "numreplicateAlert")
      createAlert(session, "alert", "numreplicateAlert", title = "Argument Input Error", content = "'Number of technical replicates' argument should be a integer.", append = TRUE)
      checknumreplicateAlert <- FALSE
    } else if (input$numreplicate<1 | input$numreplicate>10) {
      closeAlert(session, "numreplicateAlert")
      createAlert(session, "alert", "numreplicateAlert", title = "Argument Input Error", content = "'Number of technical replicates' argument should be not smaller than 1 or larger than 10.", append = TRUE)
      checknumreplicateAlert <- FALSE
    } else {
      closeAlert(session, "numreplicateAlert")
      checknumreplicateAlert <- TRUE
    }
    
    if(is.na(input$summarization_ratio)){
      closeAlert(session, "summarization_ratioAlert")
      createAlert(session, "alert", "summarization_ratioAlert", title = "Argument Input Error", content = "'Maximum missing value ratio' argument can't be empty.", append = TRUE)
      checksummarization_ratioAlert <- FALSE
    } else if (input$summarization_ratio<0 | input$summarization_ratio>1) {
      closeAlert(session, "summarization_ratioAlert")
      createAlert(session, "alert", "summarization_ratioAlert", title = "Argument Input Error", content = "'Maximum missing value ratio' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      checksummarization_ratioAlert <- FALSE
    } else {
      closeAlert(session, "summarization_ratioAlert")
      checksummarization_ratioAlert <- TRUE
    }
    
    if(is.na(input$all_missing_thresh)){
      closeAlert(session, "all_missing_threshAlert")
      createAlert(session, "alert", "all_missing_threshAlert", title = "Argument Input Error", content = "'Minimum non-missing sample ratio' argument can't be empty.", append = TRUE)
      checkall_missing_threshAlert <- FALSE
    } else if (input$all_missing_thresh<0 | input$all_missing_thresh>1) {
      closeAlert(session, "all_missing_threshAlert")
      createAlert(session, "alert", "all_missing_threshAlert", title = "Argument Input Error", content = "'Minimum non-missing sample ratio' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      checkall_missing_threshAlert <- FALSE
    } else {
      closeAlert(session, "all_missing_threshAlert")
      checkall_missing_threshAlert <- TRUE
    }
    
    if(is.na(input$rsd_filt_list)){
      closeAlert(session, "rsd_filt_listAlert")
      createAlert(session, "alert", "rsd_filt_listAlert", title = "Argument Input Error", content = "'Minimum overall variance' argument can't be empty.", append = TRUE)
      checkrsd_filt_listAlert <- FALSE
    } else {
      closeAlert(session, "rsd_filt_listAlert")
      checkrsd_filt_listAlert <- TRUE
    }
    

    if(is.na(input$group_missing_thresh)){
      closeAlert(session, "group_missing_threshAlert")
      createAlert(session, "alert", "group_missing_threshAlert", title = "Argument Input Error", content = "'Minimum non-missing sample ratio for group' argument can't be empty.", append = TRUE)
      checkgroup_missing_threshAlert <- FALSE
    } else if (input$group_missing_thresh<0 | input$group_missing_thresh>1) {
      closeAlert(session, "group_missing_threshAlert")
      createAlert(session, "alert", "group_missing_threshAlert", title = "Argument Input Error", content = "'Minimum non-missing sample ratio for group' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      checkgroup_missing_threshAlert <- FALSE
    } else {
      closeAlert(session, "group_missing_threshAlert")
      checkgroup_missing_threshAlert <- TRUE
    }
    
    if(is.na(input$pvalue_thresh)){
      closeAlert(session, "pvalue_threshAlert")
      createAlert(session, "alert", "pvalue_threshAlert", title = "Argument Input Error", content = "'P-value threshold' argument can't be empty.", append = TRUE)
      checkpvalue_threshAlert <- FALSE
    } else if (input$pvalue_thresh<0 | input$pvalue_thresh>1) {
      closeAlert(session, "pvalue_threshAlert")
      createAlert(session, "alert", "pvalue_threshAlert", title = "Argument Input Error", content = "'P-value threshold' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      checkpvalue_threshAlert <- FALSE
    } else {
      closeAlert(session, "pvalue_threshAlert")
      checkpvalue_threshAlert <- TRUE
    }
    
    if(is.na(input$fdrthresh)){
      closeAlert(session, "fdrthreshAlert")
      createAlert(session, "alert", "fdrthreshAlert", title = "Argument Input Error", content = "'False discovery threshold' argument can't be empty.", append = TRUE)
      checkfdrthreshAlert <- FALSE
    } else if (input$fdrthresh<0 | input$fdrthresh>1) {
      closeAlert(session, "fdrthreshAlert")
      createAlert(session, "alert", "fdrthreshAlert", title = "Argument Input Error", content = "'False discovery threshold' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      checkfdrthreshAlert <- FALSE
    } else {
      closeAlert(session, "fdrthreshAlert")
      checkfdrthreshAlert <- TRUE
    }
    
    if(is.na(input$foldchangethresh)){
      closeAlert(session, "foldchangethreshAlert")
      createAlert(session, "alert", "foldchangethreshAlert", title = "Argument Input Error", content = "'Fold change threshold' argument can't be empty.", append = TRUE)
      checkfoldchangethreshAlert <- FALSE
    } else if (input$foldchangethresh<0 | input$foldchangethresh>100) {
      closeAlert(session, "foldchangethreshAlert")
      createAlert(session, "alert", "foldchangethreshAlert", title = "Argument Input Error", content = "'Fold change threshold' argument should be not smaller than 0 or larger than 100.", append = TRUE)
      checkfoldchangethreshAlert <- FALSE
    } else {
      closeAlert(session, "foldchangethreshAlert")
      checkfoldchangethreshAlert <- TRUE
    }
    
    if(is.na(input$kfold)){
      closeAlert(session, "kfoldAlert")
      createAlert(session, "alert", "kfoldAlert", title = "Argument Input Error", content = "'k for k-fold Cross Validation' argument can't be empty.", append = TRUE)
      checkkfoldAlert <- FALSE
    } else if (input$kfold<1 | input$kfold>10000) {
      closeAlert(session, "kfoldAlert")
      createAlert(session, "alert", "kfoldAlert", title = "Argument Input Error", content = "'k for k-fold Cross Validation' argument should be not smaller than 1 or larger than 10000.", append = TRUE)
      checkkfoldAlert <- FALSE
    } else {
      closeAlert(session, "kfoldAlert")
      checkkfoldAlert <- TRUE
    }
    
    if(is.na(input$pls_vip_thresh)){
      closeAlert(session, "pls_vip_threshAlert")
      createAlert(session, "alert", "pls_vip_threshAlert", title = "Argument Input Error", content = "'VIP threshold' argument can't be empty.", append = TRUE)
      checkpls_vip_threshAlert <- FALSE
    } else if (input$pls_vip_thresh<1 | input$pls_vip_thresh>100) {
      closeAlert(session, "pls_vip_threshAlert")
      createAlert(session, "alert", "pls_vip_threshAlert", title = "Argument Input Error", content = "'VIP threshold' argument should be not smaller than 1 or larger than 100.", append = TRUE)
      checkpls_vip_threshAlert <- FALSE
    } else {
      closeAlert(session, "pls_vip_threshAlert")
      checkpls_vip_threshAlert <- TRUE
    }
    
    if(is.na(input$pls_ncomp)){
      closeAlert(session, "pls_ncompAlert")
      createAlert(session, "alert", "pls_ncompAlert", title = "Argument Input Error", content = "'Max number of components to consider' argument can't be empty.", append = TRUE)
      checkpls_ncompAlert <- FALSE
    } else {
      closeAlert(session, "pls_ncompAlert")
      checkpls_ncompAlert <- TRUE
    }
    
    if(is.na(input$max_comp_sel)){
      closeAlert(session, "max_comp_selAlert")
      createAlert(session, "alert", "max_comp_selAlert", title = "Argument Input Error", content = "'Number of components to use for VIP selection' argument can't be empty.", append = TRUE)
      checkmax_comp_selAlert <- FALSE
    } else {
      closeAlert(session, "max_comp_selAlert")
      checkmax_comp_selAlert <- TRUE
    }
    
    if(is.na(input$pls_permut_count)){
      closeAlert(session, "pls_permut_countAlert")
      createAlert(session, "alert", "pls_permut_countAlert", title = "Argument Input Error", content = "'Number of permutations for calculating p-values' argument can't be empty.", append = TRUE)
      checkpls_permut_countAlert <- FALSE
    } else {
      closeAlert(session, "pls_permut_countAlert")
      checkpls_permut_countAlert <- TRUE
    }
    
    if(is.na(input$max_varsel)){
      closeAlert(session, "max_varselAlert")
      createAlert(session, "alert", "max_varselAlert", title = "Argument Input Error", content = "'Max number of variables to be used' argument can't be empty.", append = TRUE)
      checkmax_varselAlert <- FALSE
    } else {
      closeAlert(session, "max_varselAlert")
      checkmax_varselAlert <- TRUE
    }
    
    if(is.na(input$abs_cor_thresh)){
      closeAlert(session, "abs_cor_threshAlert")
      createAlert(session, "alert", "abs_cor_threshAlert", title = "Argument Input Error", content = "'Absolute correlation threshold' argument can't be empty.", append = TRUE)
      checkabs_cor_threshAlert <- FALSE
    } else if (input$abs_cor_thresh<0 | input$abs_cor_thresh>1) {
      closeAlert(session, "abs_cor_threshAlert")
      createAlert(session, "alert", "abs_cor_threshAlert", title = "Argument Input Error", content = "'Absolute correlation threshold' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      checkabs_cor_threshAlert <- FALSE
    } else {
      closeAlert(session, "abs_cor_threshAlert")
      checkabs_cor_threshAlert <- TRUE
    }
    
    if(is.na(input$cor_fdrthresh)){
      closeAlert(session, "abs_cor_threshAlert")
      createAlert(session, "alert", "cor_fdrthreshAlert", title = "Argument Input Error", content = "'FDR threshold for correlation analysis' argument can't be empty.", append = TRUE)
      checkcor_fdrthreshAlert <- FALSE
    } else if (input$cor_fdrthresh<0 | input$cor_fdrthresh>1) {
      closeAlert(session, "abs_cor_threshAlert")
      createAlert(session, "alert", "cor_fdrthreshAlert", title = "Argument Input Error", content = "'FDR threshold for correlation analysis' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      checkcor_fdrthreshAlert <- FALSE
    } else {
      closeAlert(session, "abs_cor_threshAlert")
      checkcor_fdrthreshAlert <- TRUE
    }
    
    if(is.na(input$pca_cex_val)){
      closeAlert(session, "pca_cex_valAlert")
      createAlert(session, "alert", "pca_cex_valAlert", title = "Argument Input Error", content = "'Size of points on PCA plots' argument can't be empty.", append = TRUE)
      checkpca_cex_valAlert <- FALSE
    } else if (input$pca_cex_val<1 | input$pca_cex_val>20) {
      closeAlert(session, "pca_cex_valAlert")
      createAlert(session, "alert", "pca_cex_valAlert", title = "Argument Input Error", content = "'Size of points on PCA plots' argument should be not smaller than 1 or larger than 20.", append = TRUE)
      checkpca_cex_valAlert <- FALSE
    } else {
      closeAlert(session, "pca_cex_valAlert")
      checkpca_cex_valAlert <- TRUE
    }
    
    if(is.na(input$ellipse_conf_level)){
      closeAlert(session, "ellipse_conf_levelAlert")
      createAlert(session, "alert", "ellipse_conf_levelAlert", title = "Argument Input Error", content = "'Confidence interval for PCA ellipses' argument can't be empty.", append = TRUE)
      checkellipse_conf_levelAlert <- FALSE
    } else if (input$ellipse_conf_level<0 | input$ellipse_conf_level>1) {
      closeAlert(session, "ellipse_conf_levelAlert")
      createAlert(session, "alert", "ellipse_conf_levelAlert", title = "Argument Input Error", content = "'Confidence interval for PCA ellipses' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      checkellipse_conf_levelAlert <- FALSE
    } else {
      closeAlert(session, "ellipse_conf_levelAlert")
      checkellipse_conf_levelAlert <- TRUE
    }
    
    all(checknumreplicateAlert, checksummarization_ratioAlert, checkall_missing_threshAlert, checkrsd_filt_listAlert, checkgroup_missing_threshAlert,
        checkpvalue_threshAlert, checkfdrthreshAlert, checkfoldchangethreshAlert, checkkfoldAlert, checkpls_vip_threshAlert,
        checkpls_ncompAlert, checkmax_comp_selAlert, checkpls_permut_countAlert, checkmax_varselAlert, checkabs_cor_threshAlert,checkcor_fdrthreshAlert,
        checkpca_cex_valAlert, checkellipse_conf_levelAlert)
    
  })
  
  #output$siderbar<-renderUI({sidebarPanel(style="margin-left:0;",sliderInput("obs", "Slide to go to next figure:", min = 1, max = 8,  value = 1),width=3)})
 

  ##############################################################
  #
  # reactive variables
  #
  rVal <- reactiveValues()
  rVal$process <- NULL
  rVal$msg <- NULL
  rVal$obs <- NULL
  counter <- 0
  results <- list()
  dfEmpty <- data.frame(results = numeric(0))

  observe({
    if (check$count == 0)
      return()
    isolate({
      # Your logic here
    })
  })


  observeEvent(input$resetAll, {
      #removeNotification(id1)
      #      id1 <<- showNotification("Refreshing the app now. ", duration=1)
      #  js$reset()
      
            
       runjs("history.go(0);")
       #return()
      
  })
  
  if(FALSE)
  {
  observeEvent(input$resetAll, {
       
       reset("maindiv")
        reset("maindiv1")
        reset("maindiv2")
        reset("maindiv3")
         reset("maindiv4")
        reset("featselmethodi")
        reset("featselmethodii")
        reset("analysismode")
        reset("pairedanalysis")
        reset("go")
         reset("check")
               reset("inputarea")
               reset("featuretable")
               reset("classlabel")
               
               check <- reactiveValues(count = 0)
               done <- reactiveValues(count = 1)
       
        id1 <<- showNotification(paste("input: ",input$go," : ",check$count,sep=""), duration=5)
        
   })
  }
  
  
  if(FALSE)
   {
  #
  # Stop the process
  #
  observeEvent(input$stop, {
      
      id1 <<- showNotification("Stopping processing now. Please refresh the webpage to reload the app.", duration=NULL)
      
      #observeEvent(input$resetAll, {
      #    reset("form")
        #})
 
 # stopApp(returnValue = reset("form"))
 	
    #  js$navigate()
    
  stop("Program terminated by the user.")
    rVal$result <- dfEmpty
    if (!is.null(rVal$process)) {
      tools::pskill(rVal$process$pid)
      rVal$msg <- sprintf("%1$s killed", rVal$process$pid)
      rVal$process <- NULL

      if (!is.null(rVal$obs)) {
        rVal$obs$destroy()
      }
    }
  })
   }

  #
  # Handle process event
  #
  observeEvent(rVal$process, {
    rVal$obs <- observe({
      invalidateLater(500, session)
      isolate({
      result <- mccollect(rVal$process, wait = FALSE)
      if (!is.null(result)) {
        rVal$result <- result
        rVal$obs$destroy()
        rVal$process <- NULL
      }
    })
    })
  })
  
  observeEvent(input$go, 
               {
                   reset("nText2")
                   reset("nText")
                   reset("id1")
                   output$nText2 <- renderText({shiny::validate(
                     need(input$featuretable, "No datasetA provided. Please upload dataset A in 'Choose Files'."),
                     need(input$featuretable$type=="text/csv" || input$featuretable$type=="text/plain", "The format of datasetA is not correct. Please upload the file with correct format."),
                     need(input$classlabel, "No class label file provided. Please upload class label file in 'Choose Files'."),
                     need(input$classlabel$type=="text/csv" || input$classlabel$type=="text/plain", "The format of class label file is not correct. Please upload the file with correct format."),
                     need(featselmethod_check(),"No feature selection method was selected. Please select at least one method.") 
                  )})
                   shiny::validate(
                     need(input$featuretable, "No datasetA provided. Please upload dataset A in 'Choose Files'."),
                     need(input$featuretable$type=="text/csv" || input$featuretable$type=="text/plain", "The format of datasetA is not correct. Please upload the file with correct format."),
                     need(input$classlabel, "No class label file provided. Please upload class label file in 'Choose Files'."),
                     need(input$classlabel$type=="text/csv" || input$classlabel$type=="text/plain", "The format of class label file is not correct. Please upload the file with correct format."),
                     need(featselmethod_check(),"No feature selection method was selected. Please select at least one method.")
                  )
                   check$count=1
                   id1 <<- showNotification("Starting processing now. Your results will be available for download shortly. The processing time depends on the number of methods you used.", duration=NULL)
      
               })
  
  
  
  #########################################
  
  featselmethod_check <-reactive({
    if(input$analysismode=='classification' && input$pairedanalysis=='FALSE'){
      featselmethod<-input$featselmethodi
    }else{
      if(input$analysismode=='classification' && input$pairedanalysis=='TRUE'){
        featselmethod<-input$featselmethodii
      }else{
        if(input$analysismode=='regression' && input$pairedanalysis=='FALSE'){
          featselmethod<-input$featselmethodiii
        }else{
            
            if(input$analysismode=='regression' && input$pairedanalysis=='TRUE'){
              featselmethod<-input$featselmethodiv
            }else{
                featselmethod<-NULL
            }
        }
      }
    }
    featselmethod
  })
  
  
  featuretable <- reactive({
    if(input$go!=0 & check$count==1 & !is.null(input$featuretable$name) ){
      
      if((input$featuretable$type=="text/csv" || input$featuretable$type=="text/plain")){
        req(input$featuretable)
        if(input$featuretable$type=="text/plain"){
          featuretable <- read.delim(input$featuretable$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          if(input$featuretable$type=="text/csv"){
            featuretable <- read.csv(input$featuretable$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
          }
        }
        featuretable 
      }
    }else{
      
      NA
    }
  })
  
  classlabel <- reactive({
    if(input$go!=0 & check$count==1 & !is.null(input$classlabel$name) ){
      
      if((input$classlabel$type=="text/csv" || input$classlabel$type=="text/plain")){
        req(input$classlabel)
        if(input$classlabel$type=="text/plain"){
          classlabel <- read.delim(input$classlabel$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          if(input$classlabel$type=="text/csv"){
            classlabel <- read.csv(input$classlabel$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
          }
        }
        classlabel 
      }
    }else{
      
      NA
    }
  })
  
  session_outloc <- reactive({
    if(input$go!=0 & check$count==1){
      cur_date<-Sys.time()
      cur_date<-gsub(x=cur_date,pattern="-",replacement="")
      cur_date<-gsub(x=cur_date,pattern=":",replacement="")
      cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
      if(input$outloc==""){
        outloc<-paste('~/xmsPANDAout',cur_date,sep="")
      }else{
        outloc<-paste('~/',input$outloc,cur_date,sep="")
      }
      outloc
    }else{
       outloc<-paste('~/xmsPANDAout',sep="")
       outloc
    }
  })
  
 
  
  ##########################################
  
  output$nText <- renderText({
    if(input$go!=0  & check$count==1 & !is.null(featselmethod_check()) & is.data.frame(featuretable()) & is.data.frame(classlabel()) & all_alert()==TRUE){
      
      #if(input$transformation=="Log2 transformation"){
        
        #}else{
        #if(input$transformation=="Median centering"){
        # log2transform=FALSE
        # medcenter=TRUE
        # znormtransform=FALSE
        #}else{
        # log2transform=FALSE
        # medcenter=FALSE
        #  znormtransform=TRUE
        #}
        #}
      
      #  if(input$normalization=="Quantile normalization"){
   
        # }else{
        #    if(input$normalization=="Lowess normalization"){
        #  quantile_norm=FALSE
        #  lowess_norm=TRUE
        # madscaling=FALSE
        # }else{
        #quantile_norm=FALSE
          #  lowess_norm=FALSE
          #madscaling=TRUE
          # }
          #}
      
      
      
      if(input$globalcor == 'TRUE'){
        globalcor=TRUE
      }else{
        globalcor=FALSE
      }
      if(input$WGCNAmodules=='TRUE'){
        WGCNAmodules=TRUE
      }else{
        WGCNAmodules=FALSE
      }
      if(input$globalclustering=='TRUE'){
        globalclustering=TRUE
      }else{
        globalclustering=FALSE
      }
        
      if(globalcor==TRUE){
        abs_cor_thresh=input$abs_cor_thresh
        cor_fdrthresh=input$cor_fdrthresh
        cor_method=input$cor_method
        networktype=input$networktype
      }else{
        abs_cor_thresh=0.4
        cor_fdrthresh=0.2
        cor_method="spearman"
        networktype="complete"
      }
      
      if(input$use_summarizion=="TRUE"){
        summarize.replicates=TRUE
      }else{
        summarize.replicates=FALSE
      }
      
      if(input$pairedanalysis=='TRUE'){
        pairedanalysis=TRUE
      }else{
        pairedanalysis=FALSE
      }
      
      if(input$missing_val=='0'){
        missing_val=0
      }else{
        missing_val=NA
      }
      
      if(input$pca_ellipse=='TRUE'){
        pca_ellipse=TRUE
      }else{
        pca_ellipse=FALSE
      }
         
      check_pvalue_thresh<-reactive({need(input$pvalue_thresh,"error")})
      if(is.null(check_pvalue_thresh())){
        pvalue_thresh=input$pvalue_thresh
      }else{
        pvalue_thresh=0.05
      }
      
      check_foldchangethresh<-reactive({need(input$foldchangethresh,"error")})
      if(is.null(check_foldchangethresh())){
        foldchangethresh=input$foldchangethresh
      }else{
        foldchangethresh=0
      }
      
      check_fdr_method<-reactive({need(input$fdr_method,"error")})
      if(is.null(check_fdr_method())){
        fdr_method=gsub(" \\(.*\\)","",input$fdr_method)
      }else{
        fdr_method="BH"
      }
      
      check_fdrthresh<-reactive({need(input$fdrthresh,"error")})
      if(is.null(check_fdrthresh())){
        fdrthresh=input$fdrthresh
      }else{
        fdrthresh=0.2
      }
      
      check_kfold<-reactive({need(input$kfold,"error")})
      if(is.null(check_kfold())){
        kfold=input$kfold
      }else{
        kfold=10
      }
      
      check_pls_vip_thresh<-reactive({need(input$pls_vip_thresh,"error")})
      if(is.null(check_pls_vip_thresh())){
        pls_vip_thresh=input$pls_vip_thresh
      }else{
        pls_vip_thresh=2
      }
      
      check_pls_permut_count<-reactive({need(input$pls_permut_count,"error")})
      if(is.null(check_pls_vip_thresh())){
        if(input$permu_switch==TRUE){
          pls_permut_count=input$pls_permut_count
        }else{
          pls_permut_count=NA
        }
      }else{
        pls_permut_count=NA
      }
      
      check_max_varsel<-reactive({need(input$max_varsel,"error")})
      if(is.null(check_max_varsel())){
        max_varsel=input$max_varsel
      }else{
        max_varsel=100
      }
      
      featselmethod <- featselmethod_check()
      
      check_aggregation_method<-reactive({need(input$aggregation_method,"error")})
      if(is.null(check_aggregation_method())){
        aggregation_method=gsub(" \\(.*\\)","",input$aggregation_method)
        if(aggregation_method=='None'){
          aggregation_method=NA
        }
      }else{
        aggregation_method=NA
      }
      
      if(input$netbasedfeatranking_switch==FALSE){
        degree_rank_method=NA
        differential.network.analysis=FALSE
      }else{
        degree_rank_method="DiffRank"
        differential.network.analysis=TRUE
      }
      
      # if(input$timeseries.lineplots==FALSE){
          
          #   timeseries.lineplots=FALSE
          #}
      
      ###################
    
    # if(input$workflow=='workflowI')
    # if(check$count==1)
    {
      
      #start: see manual for additional arguments and description
     #start: see manual for additional arguments and description
        demetabs_res<-try(diffexp(
          #1) arguments for input files
          Xmat=featuretable(),
          parentoutput_dir=session_outloc(), #paste('~/',input$outloc,sep="")
          Ymat=classlabel(),
          feature_table_file=NA,
          class_labels_file=NA,
          input.intensity.scale=input$input_intensity_scale,
          
          ##2) data preprocessing order: 1) summarization, 2) filtering by missing values, 3) imputation; 4) transformation and normalization: halffeaturemin
          num_replicates = input$numreplicate,
          
          summarize.replicates =summarize.replicates, summary.method=input$summarization_method,summary.na.replacement=input$summary_na_replacement,
          rep.max.missing.thresh=input$summarization_ratio,
          all.missing.thresh=input$all_missing_thresh, group.missing.thresh=input$group_missing_thresh, missing.val=missing_val,
          log2transform = FALSE, medcenter=FALSE, znormtransform = FALSE,
          quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE,
          normalization.method=input$normmethod,
          TIC_norm=FALSE,
          
          rsd.filt.list = input$rsd_filt_list,
          
          ##3) arguments for feature seletion: c("limma","pls","pamr","spls","pls","MARS","RF","rfesvm","logitreg","ttest","wilcox","o1pls","lmreg")
          #"rfesvm","pamr","MARS","RF","logitreg","ttest","wilcox","o1pls","lmreg","lm1wayanova"
          #c("limma","pls","spls","pls","MARS","RF","rfesvm","logitreg","ttest","wilcox","o1pls","lmreg","lm1wayanova")
          pairedanalysis = pairedanalysis, featselmethod=featselmethod,
          pvalue.thresh=pvalue_thresh,
          fdrthresh = fdrthresh, fdrmethod=fdr_method,
          kfold=kfold,networktype=networktype,
          samplermindex=NA,numtrees=5000,analysismode=input$analysismode, pls_vip_thresh = pls_vip_thresh, num_nodes = detectCores(),
          max_varsel = max_varsel, pls_ncomp = input$pls_ncomp, pred.eval.method="BER", rocfeatlist=seq(2,10,1),
          rocfeatincrement=TRUE,
          rocclassifier="svm",foldchangethresh=foldchangethresh,
          optselect=input$optselect,max_comp_sel=input$max_comp_sel,saveRda=FALSE,pls.permut.count=pls_permut_count,
          pca.ellipse=pca_ellipse,ellipse.conf.level=input$ellipse_conf_level,svm.acc.tolerance=5,pamr.threshold.select.max=FALSE,
          aggregation.method=aggregation_method,mars.gcv.thresh=1,pls.vip.selection=input$pls.vip.selection,limmadecideTests=TRUE,
          
          #4) arguments for WGCNA and global clustering analysis (HCA and EM clustering)
          wgcnarsdthresh=30,WGCNAmodules=WGCNAmodules,globalclustering=globalclustering,
          
          #5) arguments for correlation and network analysis using the selected features
          cor.method=cor_method, abs.cor.thresh = abs_cor_thresh, cor.fdrthresh=cor_fdrthresh,
          globalcor=globalcor,target.metab.file=NA,
          target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,
          degree_rank_method=degree_rank_method, #set to NA or "DiffRank"

          #6) arguments for graphical options: see manual for additional arguments
          output.device.type="png",pca.cex.val=input$pca_cex_val,legendlocation="bottomleft",
          net_node_colors=c("green","red"),net_legend=FALSE,
          manhattanplot.col.opt=c("darkblue","red3"),
          heatmap.col.opt=input$heatmap_color_scheme,sample.col.opt=input$sample.color.theme,
          boxplot.col.opt=input$boxplot.color.theme,barplot.col.opt=input$barplot.color.theme,
          aggregation.max.iter=100,
          plots.width=8,
          plots.height=8,
          plots.type="cairo",
          boxplot.type=input$boxplot.type,
          add.jitter=input$boxplot.jitter,
          timeseries.lineplots=input$timeseries.lineplots,
          alphabetical.order=input$alphabetical.order,
          ylab_text=input$ylabel.text,
          kegg_species_code="hsa",database="pathway",reference_set=NA,match_class_dist=TRUE,differential.network.analysis=differential.network.analysis
          
        ),silent=TRUE)
        
        if(is(demetabs_res,"try-error")){
             done$count=0
              go <- reactiveValues(count = 0)
              # done <- reactiveValues(count = 0)
            #file.copy(paste(getwd(),'matrix_centrality.txt',sep='/'),session_outloc())
            print("Error processing the data. Error message:")
            print(demetabs_res)
            
        }else{
            done$count=1
             go <- reactiveValues(count = 0)
            setwd(session_outloc())
            zip(zipfile=paste(basename(session_outloc()),'zip',sep='.'), files='.')
            print("Processing complete. Please click on download button to save the results.")
            
        }
    
    }
      
      #   input$go=0
      
      #  go <- reactiveValues(count = 0)
       # go <- reactiveValues(count = 0)
       # reset("go")
       
       #   done <- reactiveValues(count = 0)
       #go <- reactiveValues(count = 0)
      
      
      
    }else{
      
      NULL
    }
  })
  
  ##########################################
  
  observeEvent({if(done$count==1) TRUE else return()},{
    if (!is.null(id1)){
      removeNotification(id1)
      id1 <<- showNotification(paste("Processing complete. Please click on download button to save the results. Output location: ",session_outloc(),sep=""), duration=10)
    }
    
    if(length(featselmethod_check())>1 & !input$aggregation_method=="none"){
      featselmethodout <- c('AggregatedResults',featselmethod_check())
    }else{
      featselmethodout <-featselmethod_check()
    }
    
    #print(c(featselmethodout,input$method))
    
    output$output_results <- renderUI({
  
      column(12,
             column(12,style='padding-top:10px;padding-left:0;',tags$div(h4("Output"))),
             column(8,align='center',style="display:block; margin-left: auto;margin-right:auto;",
                    imageOutput("myImage",width="400px",height="400px",inline=TRUE)
             ),
             column(4,
                    div(style='margin-bottom:40px;', selectInput(width="250px","methodout","Choose method to display figures:",featselmethodout)),
                    uiOutput("figureradio")
             )
      )
    
    })
    
  })
  
  observeEvent(input$methodout,{
    
    if(input$methodout=="AggregatedResults"){
      l1 <- list.files(paste(session_outloc(),'AggregatedResults',sep="/"),".png",recursive=TRUE,full.names=FALSE)
      figurenum <- paste('Figure',seq(1:length(l1)))
    }else{
      folder <- grep(input$methodout,list.dirs(paste(session_outloc(),'Stage2',sep='/'),recursive=FALSE,full.names=FALSE),value=TRUE)
      l1 <- list.files(paste(session_outloc(),'Stage2',folder,sep="/"),".png",recursive=TRUE,full.names=FALSE)
      figurenum <- paste('Figure',seq(1:length(l1)))
    }
    
    if(length(l1)>=1){
      output$figureradio <- renderUI({
        div(
          div(style='margin-bottom:10px;',tags$label("Figure Choices", `for` = "figure_choices" )),
          awesomeRadio(inputId = "figure_choices",label = NULL, choices = figurenum, selected = "Figure 1", status = "primary")
        )
        })
    }
    
    if(!is.null(input$methodout) & length(l1)>=1){
      
      output$myImage <- renderImage({
        
        if(input$methodout=="AggregatedResults"){
          req(input$figure_choices)
          filename <- normalizePath(file.path(paste(session_outloc(),'AggregatedResults',sep='/'),l1[as.numeric(gsub('Figure ','',input$figure_choices))]))
        }else{
          req(input$figure_choices)
          filename <- normalizePath(file.path(paste(paste(session_outloc(),'Stage2',sep='/'),folder,sep="/"),l1[as.numeric(gsub('Figure ','',input$figure_choices))]))
        }
        list(src = filename,width=600,height=600,
             alt = "This is an image")
        
      }, deleteFile = FALSE)
      
    }
    
  })
  

  output$downloadData <- downloadHandler(

    #if(input$go!=0 && input$featselmethod!="-" && input$feature_table_file!="" && input$class_labels_file!=""){

    filename <- function() {
      paste(basename(session_outloc()), "zip", sep=".")
    },
    content <- function(file) {
      fname1<-paste(session_outloc(),"/",basename(session_outloc()), ".zip", sep="")
      file.copy(fname1, file)
    },
    contentType = "application/zip"
    #}
  )

 

  output$downloadDatatest <- downloadHandler(
  
   #id1 <<- showNotification("Downlad",duraiton=NULL) #paste("Results location: ",session_outloc(),sep=""), duration=NULL)
  
    #if(input$go!=0 && input$featselmethod!="-" && input$feature_table_file!="" && input$class_labels_file!=""){
    
    filename <- function() {
      #paste(basename(session_outloc()), "zip", sep=".")
    	paste("dataset-", Sys.Date(), ".csv", sep="")
    },
    content <- function(file) {
  	 #id1 <<- showNotification("Downloading results now.", duration=NULL)
	 #fname1<-paste(session_outloc(),"/",basename(session_outloc()), ".zip", sep="")
      	 #file.copy(fname1, file)
#	 setwd(session_outloc())

#	 zip(zipfile=paste(basename(session_outloc()),'zip',sep='.'), files='.')  
	#write.csv("InputParameters.csv", file) 
 	 write.csv(mtcars, file)  
  }
    #contentType = "application/zip"
		
#	}
  )
  
  ##################################  Additional Analysis Page ##########################
  #######start interactive analysis
  
  
  #####end interactive analysis
  

  #######compare normalization methods
  
  normcheck2 <- reactiveValues(count = 0)
  normdone2 <- reactiveValues(count = 0)
  #normdone2 <- reactiveValues(res = list())
  #update2 <- reactiveValues(count = 0)
  normid3 <-NULL
  
  observeEvent(input$normstart,{normcheck2$count=0})
  
  
  observeEvent(input$normstart,
               {
                 output$normText4 <- renderText({shiny::validate(
                   need(input$norminput1, "No data file was provided. Please upload your data file."),
                   need(input$norminput1$type=="text/csv" || input$norminput1$type=="text/plain", "The format of data file is not correct. Please upload the file with correct format.")
                 )})
                 shiny::validate(
                   need(input$norminput1, "No data file was provided. Please upload your data file."),
                   need(input$norminput1$type=="text/csv" || input$norminput1$type=="text/plain", "The format of data file is not correct. Please upload the file with correct format.")
                 )
                 normcheck2$count=1
                 normid3 <<- showNotification("Data is processing now.", duration=NULL)
               })
  
  
  norm_metab_data <- reactive({
    
    if(input$normstart!=0  & normcheck2$count==1){
      
      if(input$norminput1$type=="text/plain"){
        metab_data <- read.delim(input$norminput1$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        
      }else{
        if(input$norminput1$type=="text/csv"){
          metab_data <- read.csv(input$norminput1$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
          
        }else{
          
          metab_data <- NULL
          
        }
       }
       metab_data
    }
    
  })
  
  norm_class_data <- reactive({
     
     if(input$normstart!=0  & normcheck2$count==1){
       
       if(input$norminput2$type=="text/plain"){
         
         class_data <- read.delim(input$norminput2$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
       }else{
         if(input$norminput2$type=="text/csv"){
           
           class_data <- read.csv(input$norminput2$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
         }else{
           
           
           class_data<-NULL
         }
        }
        class_data
     }
     
   })
  
  session_outloc2 <- reactive({
     if(input$normstart!=0  & normcheck2$count==1){
       cur_date<-Sys.time()
       cur_date<-gsub(x=cur_date,pattern="-",replacement="")
       cur_date<-gsub(x=cur_date,pattern=":",replacement="")
       cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
       
         outloc<-paste('~/xmsPANDAnormresults',cur_date,sep="")
      
       outloc
     }else{
       cur_date<-Sys.time()
       cur_date<-gsub(x=cur_date,pattern="-",replacement="")
       cur_date<-gsub(x=cur_date,pattern=":",replacement="")
       cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
       outloc<-paste('~/xmsPANDAnormresults',cur_date,sep="")
       outloc
     }
     
     # normid3 <<- showNotification("Data is starting now.", duration=NULL)
   })
  
  output$normText5 <- renderText({
    
    if(input$normstart!=0  & normcheck2$count==1){
        
        if(input$abs.cor.thresh==(-1)){
            
            abs.cor.thresh=NA
        }else{
            
            abs.cor.thresh=input$abs.cor.thresh
        }

   
   res<-compare.normalization(Xmat=norm_metab_data(),Ymat=norm_class_data(),Zmat=NA,feature_table_file=NA,parentoutput_dir=session_outloc2(),class_labels_file=NA,num_replicates=1,feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method="mean",
       all.missing.thresh=0.5,group.missing.thresh=0.7,
   missing.val=0,samplermindex=NA, rep.max.missing.thresh=0.5,summary.na.replacement=input$summary_na_replacement,featselmethod,pairedanalysis=FALSE,
   normalization.method=input$normalization.method,input.intensity.scale="raw",
   abs.cor.thresh=abs.cor.thresh,pvalue.thresh=input$pvalue.thresh,cor.fdrthresh=input$cor.fdrthresh,cex.plots=0.7,plots.width=8,plots.height=8,plots.res=600,
   plots.type="cairo",heatmap.col.opt="RdBu",cor.method="spearman",pca.ellipse=FALSE,ground_truth_file=NA,cutree.method="default",rsd.filt.thresh=1,alphabetical.order=TRUE,
   analysistype=input$analysistype,lme.modeltype="RI",study.design=input$study.design,log2.transform.constant=input$log2.transform.constant)
   
   
      normdone2$count=1
      
      setwd(session_outloc2())
      zip(zipfile=paste(basename(session_outloc2()),'zip',sep='.'), files='.')
      print("Processing complete. Please click on download button to save the results.")

    reset("normstart")

    }else{

      NULL
    }

  })
  
  observeEvent({if(normdone2$count==1) TRUE else return()},{
  if (!is.null(normid3)){
    removeNotification(normid3)
    normid3 <<- showNotification("Processing complete. Please click on download button to save the results.", duration=NULL)
  }
  })
  
  output$downloadnormData <- downloadHandler(
    
    #if(input$go!=0 && input$featselmethod!="-" && input$feature_table_file!="" && input$class_labels_file!=""){
    
    filename <- function() {
      paste(basename(session_outloc2()), "zip", sep=".")
    },
    content <- function(file) {
      fname1<-paste(session_outloc2(),"/",basename(session_outloc2()), ".zip", sep="")
      file.copy(fname1, file)
    },
    contentType = "application/zip"
    #}
  )
  #####end compare normalizaito methods
  
  pcacheck2 <- reactiveValues(count = 0)
   pcadone2 <- reactiveValues(count = 0)
   #normdone2 <- reactiveValues(res = list())
   #update2 <- reactiveValues(count = 0)
   pcaid3 <-NULL
   
   observeEvent(input$pcastart,{pcacheck2$count=0})
   
   
   observeEvent(input$pcastart,
                {
                  output$pcaText4 <- renderText({shiny::validate(
                    need(input$pcainput1, "No data file was provided. Please upload your data file."),
                    need(input$pcainput1$type=="text/csv" || input$pcainput1$type=="text/plain", "The format of data file is not correct. Please upload the file with correct format.")
                  )})
                  shiny::validate(
                    need(input$pcainput1, "No data file was provided. Please upload your data file."),
                    need(input$pcainput1$type=="text/csv" || input$pcainput1$type=="text/plain", "The format of data file is not correct. Please upload the file with correct format.")
                  )
                 pcacheck2$count=1
                  pcaid3 <<- showNotification("Data is processing now.", duration=NULL)
                })
   
   
  pca_metab_data <- reactive({
     
     if(input$pcastart!=0  & pcacheck2$count==1){
       
       if(input$pcainput1$type=="text/plain"){
         metab_data <- read.delim(input$pcainput1$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
         
       }else{
         if(input$pcainput1$type=="text/csv"){
           metab_data <- read.csv(input$pcainput1$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
           
         }else{
           
           metab_data <- NULL
           
         }
        }
        metab_data
     }
     
   })
   
   pca_class_data <- reactive({
      
      if(input$pcastart!=0  & pcacheck2$count==1){
        
        if(input$pcainput2$type=="text/plain"){
          
          class_data <- read.delim(input$pcainput2$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          if(input$pcainput2$type=="text/csv"){
            
            class_data <- read.csv(input$pcainput2$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
          }else{
            
            
            class_data<-NULL
          }
         }
         class_data
      }
      
    })
   
   session_outloc3 <- reactive({
      if(input$pcastart!=0  & pcacheck2$count==1){
        cur_date<-Sys.time()
        cur_date<-gsub(x=cur_date,pattern="-",replacement="")
        cur_date<-gsub(x=cur_date,pattern=":",replacement="")
        cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
        
          outloc<-paste('~/xmsPANDApcaresults',cur_date,sep="")
       
        outloc
      }else{
        cur_date<-Sys.time()
        cur_date<-gsub(x=cur_date,pattern="-",replacement="")
        cur_date<-gsub(x=cur_date,pattern=":",replacement="")
        cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
        outloc<-paste('~/xmsPANDApcaresults',cur_date,sep="")
        outloc
      }
      
      # normid3 <<- showNotification("Data is starting now.", duration=NULL)
    })
   

  output$pcaText5 <- renderText({
    
    if(input$pcastart!=0  & pcacheck2$count==1){

   
  
   
   pcares<-get_pcascoredistplots(X=pca_metab_data(),Y=pca_class_data(),feature_table_file=NA,parentoutput_dir=session_outloc3(),class_labels_file=NA,sample.col.opt=input$pca.sample.color.theme,
   plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=NA,pairedanalysis=input$pca.pairedanalysis,pca.cex.val=2,legendlocation="topright",
   pca.ellipse=TRUE,ellipse.conf.level=0.95,filename="PCA_results",paireddesign=NA,newdevice=TRUE,
   lineplot.col.opt=input$pca.lineplot.color.theme,lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),timeseries.lineplots=input$pca.timeseries.lineplots,pcacenter=TRUE,pcascale=TRUE,alphabetical.order=input$pca.alphabetical.order,analysistype=input$pca.analysistype,lme.modeltype=lme.modeltype) #,silent=TRUE)
   
 
      
      pcadone2$count=1
      
      setwd(session_outloc3())
     
  
         zip(zipfile=paste(basename(session_outloc3()),'zip',sep='.'), files='.')
         
      print("Processing complete. Please click on download button to save the results.")

    reset("pcastart")
      

    }else{

      NULL
    }

  })
  
  
  observeEvent({if(pcadone2$count==1) TRUE else return()},{
    if (!is.null(pcaid3)){
      removeNotification(pcaid3)
      pcaid3 <<- showNotification("Processing complete. Click on Download Results to save the results.", duration=NULL)
      
      
    }
    
   
   # column(12,style='padding-top:10px;padding-left:0;',tags$div(h4("Output")))
    output$pcaoutput_results <- renderText({
  pcal1 <- list.files(paste(session_outloc3()),".pdf",recursive=FALSE,full.names=FALSE)
    
    print(pcal1)
    
    filename <- normalizePath(file.path(session_outloc3(),pcal1))
    print(filename)
    
                 PDFfile=filename
                 print(paste("file exists:",file.exists(PDFfile)))
                 list(src=PDFfile)
                 
             
             #     return(paste('<iframe style="height:600px; width:100%" src="', filename, '"></iframe>', sep = ""))
             
             
      
    
    }) #,deleteFile=FALSE)
    
  })
  
  
  
  output$downloadpcaData <- downloadHandler(
    
    #if(input$go!=0 && input$featselmethod!="-" && input$feature_table_file!="" && input$class_labels_file!=""){
    
    filename <- function() {
      paste(basename(session_outloc3()), "zip", sep=".")
    },
    content <- function(file) {
      fname1<-paste(session_outloc3(),"/",basename(session_outloc3()), ".zip", sep="")
      file.copy(fname1, file)
    },
    contentType = "application/zip"
    #}
  )
  
  
  ########################
  
  check2 <- reactiveValues(count = 0)
  #update2 <- reactiveValues(count = 0)
  done2 <- reactiveValues(cluster_table = matrix())
  id3 <- NULL
  #off <- reactiveVal(0)
  #observeEvent(input$clusterinput,{update2$count=1;start2$count=0})
  #observeEvent(input$kegg_species_code,{update2$count=1;start2$count=0})
  #observeEvent(input$database,{update2$count=1;start2$count=0})
  #observeEvent(input$type.statistic,{update2$count=1;start2$count=0})
  observeEvent(input$start2,{check2$count=0})
  
  
  observeEvent(input$start2, 
               {
                 output$nText4A <- renderText({shiny::validate(
                   need(input$clusterinput, "No data file was provided. Please upload your data file."),
                   need(input$clusterinput$type=="text/csv" || input$clusterinput$type=="text/plain", "The format of data file is not correct. Please upload the file with correct format.")
                 )})
                 shiny::validate(
                   need(input$clusterinput, "No data file was provided. Please upload your data file."),
                   need(input$clusterinput$type=="text/csv" || input$clusterinput$type=="text/plain", "The format of data file is not correct. Please upload the file with correct format.")
                 )
                 check2$count=1
                 id3 <<- showNotification(paste("Data is processing now using database ",input$fcs.database,sep=""), duration=NULL)
               })
  
  
  cluster_metab_data <- reactive({
    
    if(input$start2!=0  & check2$count==1){
      
      if(input$clusterinput$type=="text/plain"){
        metab_data <- read.delim(input$clusterinput$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
      }else{
        if(input$clusterinput$type=="text/csv"){
          metab_data <- read.csv(input$clusterinput$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          
          metab_data <- NULL
        }
       }
       metab_data
    }
    
  })
  
  cluster_metab_data2 <- reactive({
    
    if(input$start2!=0  & check2$count==1 & input$fcs.database=='custom'){
      
      if(input$clusterinput2$type=="text/plain"){
        metab_data <- read.delim(input$clusterinput2$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
      }else{
        if(input$clusterinput2$type=="text/csv"){
          metab_data <- read.csv(input$clusterinput2$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          
          metab_data <- NULL
        }
       }
       metab_data
    }
    
  })
  
  cluster_metab_data3 <- reactive({
    
    if(input$start2!=0  & check2$count==1 & input$fcs.upload.annotation.file=='TRUE'){
      
      if(input$clusterinput3$type=="text/plain"){
        metab_data <- read.delim(input$clusterinput3$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
      }else{
        if(input$clusterinput3$type=="text/csv"){
          metab_data <- read.csv(input$clusterinput3$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          
          metab_data <- NULL
        }
       }
       metab_data
    }
    
  })
  
  output$nText5A <- renderText({
    
    if(input$start2!=0  & check2$count==1){

      kegg_species_code <- switch(isolate(input$kegg_species_code), "Homo sapiens(default)" = "hsa",
                                  "Mus musculus" = "mmu",
                                  "Pan troglodytes" = "ptr",
                                  "Macaca mulatta" = "mcc",
                                  "Bos taurus" = "bta",
                                  "Rattus norvegicus" = "rno",
                                  "Danio rerio"= "dre",
                                  "C. elegans"= "cel",
                                  "Drosophila melanogaster"= "dme"
      )

      # if(isolate(input$database)=='pathway(default)'){
      #  database="pathway"
      #}else{
      # database="module"
      #}

      database=input$fcs.database
        
      if(isolate(input$type.statistic)=='TRUE'){
        type.statistic="p-value"
      }else{
        type.statistic="other"
      }
      
      
    
       
      
      if(isolate(input$fcs.upload.annotation.file) == 'TRUE'){
          target.data.annot=cluster_metab_data3()
          
      }else{
          target.data.annot="none"
          
      }
      #if(isolate(input$fcs.database)=='custom'){
      # checkfcsdatabase=1
      # }else{
        
        #     checkfcsdatabase=0
        #}
      #print("yes1")
      #print(head(cluster_metab_data()))
      #print(c(kegg_species_code,database,type.statistic))
      #done2$cluster_table <- ctable_module
      
      if(database=="custom"){
          done2$cluster_table <-get_fcs(target.data=cluster_metab_data(),target.data.annot=target.data.annot,kegg_species_code=kegg_species_code,database="custom",type.statistic=type.statistic,fcs.min.hits=input$fcs.min.hits,reference_set=cluster_metab_data2(),itrs=input$fcs.itrs)
          
      }else{
          done2$cluster_table <-get_fcs(target.data=cluster_metab_data(),target.data.annot=target.data.annot,kegg_species_code=kegg_species_code,database=input$fcs.database,type.statistic=type.statistic,fcs.min.hits=input$fcs.min.hits,reference_set=NA,itrs=input$fcs.itrs)
      }
      
      NULL
      #checktable1=1
      #print(dim(done2$cluster_table))
      #print("yes2")
      reset("start2")
    }else{

      NULL
    }

  })
  
  output$checkfcsdatabase <- reactive({
      if(input$fcs.database == 'custom'){
       1
     }else{
       0
     }
  })
  
  
  
  outputOptions(output, "checkfcsdatabase", suspendWhenHidden = FALSE)
  
  output$checkannotfile <- reactive({
         if(input$fcs.upload.annotation.file == 'TRUE'){
           1
         }else{
           0
         }
      })
      
       
       outputOptions(output, "checkannotfile", suspendWhenHidden = FALSE)
       
 
  observeEvent({if(dim(done2$cluster_table)[2]>1) TRUE else return()},{
    if (!is.null(id3)){
      removeNotification(id3)
      id3 <<-  showNotification("Processing complete. Click on View Results to view or download the results.", duration=NULL)
    }
    
    table1 = done2$cluster_table
    table1[,1] <- as.character(table1[,1])
    
    table1$XID<-gsub(as.character(table1$XID),pattern=";",replacement=";\t")
    if(input$fcs.database=="pathway"){
    map_cpd<-paste(table1[,1],"+",as.character(table1$XID),sep="")
    
    temp_pattern=paste(" ",input$path.bg.color,"+",sep="")
    map_cpd<-gsub(as.character(map_cpd),pattern=";",replacement=temp_pattern)
    map_cpd<-paste(map_cpd," ",input$path.bg.color,sep="")
    
    
     table1[,1] <- paste("<a target='_blank' href='https://www.genome.jp/kegg-bin/show_pathway?",map_cpd,"'>",table1[,1],"</a>",sep="")
    }else{
        
        if(input$fcs.database=="brite"){
            
            map_cpd<-paste(table1[,1],".keg+",as.character(table1$XID),sep="")
              
              temp_pattern=paste(" ",input$path.bg.color,"+",sep="")
              map_cpd<-gsub(as.character(map_cpd),pattern=";",replacement=temp_pattern)
              map_cpd<-paste(map_cpd," ",input$path.bg.color,sep="")
              
              
            table1[,1] <- paste("<a target='_blank' href='https://www.genome.jp/kegg-bin/get_htext?",map_cpd,"'>",table1[,1],"</a>",sep="")
        }else{
            
            
            
            if(input$fcs.database=="module"){
                       
                       map_cpd<-paste(table1[,1],"+",as.character(table1$XID),sep="")
                         
                         temp_pattern=paste("+",sep="")
                         map_cpd<-gsub(as.character(map_cpd),pattern=";",replacement=temp_pattern)
                         # map_cpd<-paste(map_cpd," ",input$path.bg.color,sep="")
                         
                         
                       table1[,1] <- paste("<a target='_blank' href='https://www.genome.jp/kegg-bin/show_module?",map_cpd,"'>",table1[,1],"</a>",sep="")
                   }
            
        }
        
    }
    all=dim(table1)[1]
               if(all>1000){
                 lines=c(5,10,50,100,500,1000,all)
               }else{
                 lines=c(5,10,50,100,500,1000)
               }
               
    output$pathwaytb = renderDT(
    #done2$cluster_table, options = list(pageLength = 5,lengthChange = FALSE)
    table1,options = list(dom='lrtip', lengthMenu = lines), rownames=FALSE, escape = FALSE)
     
    
  })
  
  
  
  output$pathwaytb1 <- renderDataTable(
    
    #if(dim(done2$cluster_table)[2]>1)
    {
        #print(done2$cluster_table)
      table = done2$cluster_table
      table[,1] <- as.character(table[,1])
      
      for(i in 1:dim(table)[1]){

          if(isolate(input$fcs.database)=='pathway'){
          table[i,1]=paste("<a target='_blank' href='https://www.genome.jp/kegg-bin/show_pathway?",table[i,1],"'>",table[i,1],"</a>",sep="")
        }else{

            #  moduleid = strsplit(table[i,1],split = ":")[[1]][length(strsplit(table[i,1],split = ":")[[1]])]
            #table[i,1]=paste("<a target='_blank' #href='https://www.kegg.jp/kegg-bin/show_module?",moduleid,"'>",table[i,1],"</a>",sep="")
        }

      }
      
      all=dim(table)[1]
      if(all>1000){
        lines=c(5,10,50,100,500,1000,all)
      }else{
        lines=c(5,10,50,100,500,1000)
      }
      
      datatable(table,options = list(dom='lrtip', lengthMenu = lines), rownames=FALSE, escape = FALSE)
    }
    
  )
  
  output$downloadtableData <- downloadHandler(
    filename = function() {
        if(input$fcs.database=='pathway'){
          paste("pathway_table", ".csv", sep = "")
        }else{
          paste("custom_table", ".csv", sep = "")
        }
      },
      content = function(file) {
        write.csv(done2$cluster_table, file, row.names = FALSE)
      }
  )
  # output$downloadPlot3 <- downloadHandler(
  #   
  #   filename <- "barplot",
  #   content = function(file) {
  #     
  #     if(is.na(input$figurewidth)){
  #       width=10
  #     }else{
  #       width=input$figurewidth
  #     }
  #     
  #     if(is.na(input$figureheight)){
  #       height=6
  #     }else{
  #       height=input$figureheight
  #     }
  #     ggsave(file, plot = plot5(), device = "png", width = width, height = height, units='in')
  #   }
  # )
  
  
  #output$checkplot3 <- reactive({!is.null(pplot3())})
  #outputOptions(output, "checkplot3", suspendWhenHidden = FALSE)
  output$checktable1 <- reactive(
    {
        if(dim(done2$cluster_table)[2]>1){
            
            1
        }else{
            
            0
        }
        
    }
  )
  outputOptions(output, "checktable1", suspendWhenHidden = FALSE)
  
  
  output$downloadbutton <- renderUI({
        column(12,style="margin-top:20px;text-align:right;",
               downloadButton(style = "background-color:#417ee0;color:#ffffff;",outputId = "downloadtableData", label = "Download Table")
               )
    })
  
  #output$pathwaytb <- renderTable({pathway_table()[,-1]}, striped = TRUE) 
  
  #output$hover <- renderPrint(list(input$figurewidth,input$figureheight))
  
  ########################

  check3 <- reactiveValues(count = 0)
  done3 <- reactiveValues(count = 0)
  id4 <- NULL
  #off <- reactiveVal(0)
  observeEvent(input$start3,{check3$count=0})
  
  
  observeEvent(input$start3, 
               {
                 output$nText6 <- renderText({shiny::validate(
                   need(input$featuretable_file, "No feature table provided. Please upload your feature table."),
                   need(input$featuretable_file$type=="text/csv" || input$featuretable_file$type=="text/plain", "The format of feature table is not correct. Please upload the file with correct format."),
                   need(input$classlabel_file, "No class label file was provided. Please upload your class label file."),
                   need(input$classlabel_file$type=="text/csv" || input$classlabel_file$type=="text/plain", "The format of class label file is not correct. Please upload the file with correct format."),
                   need(input$ref_meta_file, "No standard metabolite library file was provided. Please upload your data file."),
                   need(input$ref_meta_file$type=="text/csv" || input$ref_meta_file$type=="text/plain", "The format of standard metabolite library file is not correct. Please upload the file with correct format."),
                   need(any(c(input$step1,input$step2,input$step3)),"Please select at least one step."),
                   if(input$step3){
                     need(input$foldchange_file, "No fold change file was provided. Please upload your data file.")
                   },
                   if(input$step3){
                     need(input$foldchange_file$type=="text/csv" || input$foldchange_file$type=="text/plain", "The format of fold change file is not correct. Please upload the file with correct format.")
                   }
                 )})
                 shiny::validate(
                   need(input$featuretable_file, "No feature table provided. Please upload your feature table."),
                   need(input$featuretable_file$type=="text/csv" || input$featuretable_file$type=="text/plain", "The format of feature table is not correct. Please upload the file with correct format."),
                   need(input$classlabel_file, "No class label file was provided. Please upload your class label file."),
                   need(input$classlabel_file$type=="text/csv" || input$classlabel_file$type=="text/plain", "The format of class label file is not correct. Please upload the file with correct format."),
                   need(input$ref_meta_file, "No standard metabolite library file was provided. Please upload your data file."),
                   need(input$ref_meta_file$type=="text/csv" || input$ref_meta_file$type=="text/plain", "The format of standard metabolite library file is not correct. Please upload the file with correct format."),
                   need(any(c(input$step1,input$step2,input$step3)),"Please select at least one step."),
                   if(input$step3){
                     need(input$foldchange_file, "No fold change file was provided. Please upload your data file.")
                   },
                   if(input$step3){
                     need(input$foldchange_file$type=="text/csv" || input$foldchange_file$type=="text/plain", "The format of fold change file is not correct. Please upload the file with correct format.")
                   }
                 )
                 check3$count=1
                 id4 <<- showNotification("Data is processing now.", duration=NULL)
               })
  
  
  featuretable_file <- reactive({
    if(input$start3!=0  & check3$count==1 & !is.null(input$featuretable_file$name) ){
      
      if((input$featuretable_file$type=="text/csv" || input$featuretable_file$type=="text/plain")){
        req(input$featuretable_file)
        if(input$featuretable_file$type=="text/plain"){
          featuretable_file <- read.delim(input$featuretable_file$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          if(input$featuretable_file$type=="text/csv"){
            featuretable_file <- read.csv(input$featuretable_file$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
          }
        }
        featuretable_file 
        
      }
    }
  })
  
  classlabel_file <- reactive({
    if(input$start3!=0 & check3$count==1 & !is.null(input$classlabel_file$name) ){
      
      if((input$classlabel_file$type=="text/csv" || input$classlabel_file$type=="text/plain")){
        req(input$classlabel_file)
        if(input$classlabel_file$type=="text/plain"){
          classlabel_file <- read.delim(input$classlabel_file$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          if(input$classlabel_file$type=="text/csv"){
            classlabel_file <- read.csv(input$classlabel_file$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
          }
        }
        classlabel_file 
        
      }
    }
  })
  
  ref_meta_file <- reactive({
    if(input$start3!=0 & check3$count==1 & !is.null(input$ref_meta_file$name) ){
      
      if((input$ref_meta_file$type=="text/csv" || input$ref_meta_file$type=="text/plain")){
        req(input$ref_meta_file)
        if(input$ref_meta_file$type=="text/plain"){
          ref_meta_file <- read.delim(input$ref_meta_file$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          if(input$ref_meta_file$type=="text/csv"){
            ref_meta_file <- read.csv(input$ref_meta_file$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
          }
        }
        ref_meta_file 
        
      }
    }
  })
  
  foldchange_file <- reactive({
    if(input$start3!=0 & check3$count==1 & input$step3 & !is.null(input$foldchange_file$name) ){
      
      if((input$foldchange_file$type=="text/csv" || input$foldchange_file$type=="text/plain")){
        req(input$foldchange_file)
        if(input$foldchange_file$type=="text/plain"){
          foldchange_file <- read.delim(input$foldchange_file$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          if(input$foldchange_file$type=="text/csv"){
            foldchange_file <- read.csv(input$foldchange_file$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
          }
        }
        foldchange_file 
        
      }
    }else{
      
      NA
    }
  })
  
  session_outloc_quant <- reactive({
    if(input$start3!=0 & check3$count==1){
      cur_date<-Sys.time()
      cur_date<-gsub(x=cur_date,pattern="-",replacement="")
      cur_date<-gsub(x=cur_date,pattern=":",replacement="")
      cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
      outloc<-paste('~/metabolite_quantification_analysis_results',cur_date,sep="")
      outloc
    }else{
      NULL
    }
  })
  
  
  output$nText7 <- renderText({
    
    if(input$start3!=0  & check3$count==1){
      
    isolate({
      steps=""
      if(input$step1){
        steps=paste(steps,"1",sep="")
        if(input$groupcheck=="TRUE"){
          groupcheck=TRUE
        }else{
          groupcheck=FALSE
        }
        if(input$targetID==""){
          targetID=NA
        }else{
          targetID=unlist(strsplit(input$targetID,split=","))
        }
      }else{
        groupcheck=FALSE
        targetID=NA
      }
      if(input$step2){
        steps=paste(steps,"2",sep="")
      }
      if(input$step3){
        steps=paste(steps,"3",sep="")
        highcolor=input$highcolor
        lowcolor=input$lowcolor
        if(!is.na(input$minhit)){
          minhit <- input$minhit
        }else{
          stop("Please enter the correct value for 'Minimum #metablites hits in KEGG map'.")
        }
      }else{
        minhit=3
        highcolor="red"
        lowcolor="blue"
      }
      
      if(input$summarize_replicates=="TRUE"){
        summarize_replicates=TRUE
      }else{
        summarize_replicates=FALSE
      }
      
      if(!is.na(input$num_replicate2)){
        num_replicate2 <- input$num_replicate2
      }else{
        stop("Please enter the correct value for 'Number of technical replicates'.")
      }
      
      if(!is.na(input$rep_max_missing_thresh)){
        rep_max_missing_thresh <- input$rep_max_missing_thresh
      }else{
        stop("Please enter the correct value for 'Maximum missing value ratio'.")
      }
      
      if(!is.na(input$mass_error)){
        mass_error <- input$mass_error
      }else{
        stop("Please enter the correct value for 'Mass-to-charge tolerance'.")
      }
      
      if(!is.na(input$time_error)){
        time_error <- input$time_error
      }else{
        stop("Please enter the correct value for 'Retention time tolerance'.")
      }
      
      quant(Xmat=featuretable_file(),Ymat=classlabel_file(),Wmat=ref_meta_file(),Zmat=foldchange_file(),
            feature_table=NA,class_file=NA,ref_list=NA,foldchange_list=NA,
            outloc=session_outloc_quant(),
            num_replicates=num_replicate2,
            summarize_replicates=summarize_replicates,
            rep.max.missing.thresh=rep_max_missing_thresh,
            summary.method=input$summary_method,
            mass_error= mass_error,
            time_error= time_error,
            percent_node=1,
            steps=steps,
            min_num_nonmissing=3,
            targetID=targetID,
            minhit=minhit,
            groupcheck=groupcheck,
            highcolor=highcolor,
            lowcolor=lowcolor) 
      
      
      done3$count=1
      #file.copy(paste(getwd(),'matrix_centrality.txt',sep='/'),session_outloc())
      setwd(session_outloc_quant())
      zip(zipfile=paste(basename(session_outloc_quant()),'zip',sep='.'), files='.')
      print("Processing complete. Please click on download button to save the results.")
      reset("start3")
    })  
      
    }else{
      
      NULL
    }
      
  })
  
  
  observeEvent({if(done3$count==1) TRUE else return()},{
    if (!is.null(id4)){
      removeNotification(id4)
      id4 <<- showNotification("Processing complete. Click on download button to save the results.", duration=NULL)
    }
  })
  
  output$downloadQdata <- downloadHandler(
    
    #if(input$go!=0 && input$featselmethod!="-" && input$feature_table_file!="" && input$class_labels_file!=""){
    
    filename <- function() {
      paste(basename(session_outloc_quant()), "zip", sep=".")
    },
    content <- function(file) {
      fname1<-paste(session_outloc_quant(),"/",basename(session_outloc_quant()), ".zip", sep="")
      file.copy(fname1, file)
    },
    contentType = "application/zip"
    #}
  )
  
  output$done <- reactive({done3$count==1})
  outputOptions(output, "done", suspendWhenHidden = FALSE)
  
  ##################################  Help Page #################################################  
 
  example_feat <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/feature_table_one_or_two_factor_analysis.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_feat <- example_feat[1:5,1:7]
  colnames(example_feat) <- c(colnames(example_feat)[1:6],"...")
  output$example_feat <- renderTable({ example_feat }, striped = TRUE)
  
  example_feat2 <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/feature_table_one_or_two_factor_analysis_withNamesorIDs.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_feat2 <- example_feat2[1:5,1:7]
  colnames(example_feat2) <- c(colnames(example_feat2)[1:6],"...")
  output$example_feat2 <- renderTable({ example_feat2 }, striped = TRUE)
  
  example_multiclass_comparison <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/classlabels_multiclass_comparison.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_multiclass_comparison <- rbind(head(example_multiclass_comparison[example_multiclass_comparison$Factor1=="Group1",],8),head(example_multiclass_comparison[example_multiclass_comparison$Factor1=="Group2",],8))
  example_multiclass_comparison_covariates <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/classlabels_with_covariates.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_multiclass_comparison_covariates <- rbind(head(example_multiclass_comparison_covariates[example_multiclass_comparison_covariates$Class=="NonSmoker",],8),head(example_multiclass_comparison_covariates[example_multiclass_comparison_covariates$Class=="Smoker",],8))
  example_regression <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/classlabels_regression.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_regression <- example_regression[1:16,]
  example_two_way_anova <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/classlabels_two_way_anova.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_two_way_anova <- rbind(head(example_two_way_anova[example_two_way_anova$Factor1=="Group1",],8),head(example_two_way_anova[example_two_way_anova$Factor1=="Group2",],8))
  example_one_factor_repeatedmeasures <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/classlabels_one_factor_repeatedmeasures.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_one_factor_repeatedmeasures <- example_one_factor_repeatedmeasures[1:16,]
  example_two_factor_repeatedmeasures <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/classlabels_two_factor_repeatedmeasures.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_two_factor_repeatedmeasures <- example_two_factor_repeatedmeasures[1:16,]
  
  
  
  output$example_classlabel_text <- renderUI({
    txt <- switch(input$classlabel_option, 
                'multiclass comparison' = tags$p(style="text-align:justify;font-size: 15px;","To do multiclass comparison, the class labels file should include two columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the factor1 (or class1) which contains the group information."),
                'multiclass comparison with covariates'= tags$p(style="text-align:justify;font-size: 15px;","To do multiclass comparison adjusted for some covariates, the class labels file should include sampleID, group information and covariates. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the class which contains the group information. The remaing columns should be covariates information."),
                'regression' = tags$p(style="text-align:justify;font-size: 15px;","To do the regression analysis (adjusted for covariates), the class labels file should include at least two columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the response variable which should be numeric variable. If you want to adjust for some covariates, put the covariates information after the first 2 columns."),
                'two-way anova' = tags$p(style="text-align:justify;font-size: 15px;","To do the two-way ANOVA analysis, the class labels file should include three columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the factor1 information and the third column is the factor2 information."),
                'one factor repeatedmeasures' = tags$p(style="text-align:justify;font-size: 15px;","To do the one-way ANOVA analysis with repeated-measurement data, the class labels file should include three columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the subject ID and the third column is the factor1 which contains the group information."),
                'two factor repeatedmeasures' = tags$p(style="text-align:justify;font-size: 15px;","To do the two-way ANOVA analysis with repeated-measurement data, the class labels file should include four columns. The first column is the the sample ID (or filename) which should correspond to the column name of each sample in feature table. The second column is the subject ID. The third and the fourth column are the factor1 and factor2 information.")
    )
    txt
  })
  
  output$example_classlabel_table <- renderTable({
    tb <- switch(input$classlabel_option, 
           'multiclass comparison' = example_multiclass_comparison,
           'multiclass comparison with covariates'= example_multiclass_comparison_covariates,
           'regression' = example_regression,
           'two-way anova' = example_two_way_anova,
           'one factor repeatedmeasures' = example_one_factor_repeatedmeasures,
           'two factor repeatedmeasures' = example_two_factor_repeatedmeasures
    )
    tb
  }, striped = TRUE) 
  
  #output$example_multiclass_comparison <- renderTable({ example_multiclass_comparison }, striped = TRUE) 
  #output$example_multiclass_comparison_covariates <- renderTable({ example_multiclass_comparison_covariates }, striped = TRUE) 
  #output$example_regression <- renderTable({ example_regression }, striped = TRUE) 
  #output$example_two_way_anova <- renderTable({ example_two_way_anova }, striped = TRUE) 
  #output$example_one_factor_repeatedmeasures <- renderTable({ example_one_factor_repeatedmeasures }, striped = TRUE)
  #output$example_two_factor_repeatedmeasures <- renderTable({ example_two_factor_repeatedmeasures }, striped = TRUE)
  
  #example_one_factor_repeatedmeasures_2timepoints <- read.delim("example_data/classlabels_one_factor_repeatedmeasures_2timepoints.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  #example_one_factor_repeatedmeasures_2timepoints <- example_one_factor_repeatedmeasures_2timepoints[1:16,]
  #output$example_one_factor_repeatedmeasures_2timepoints <- renderTable({ example_one_factor_repeatedmeasures_2timepoints }, striped = TRUE)
  
  #output$table_dataset <- renderPrint({
  
  #library(htmlTable)
  #library('xMWAS')        
  #data(exnci60)        
  
  #output<-exnci60$mrna[1:5,1:6]
  #output<-htmlTable(output,align='c',caption="<strong>Dataset File Format:</strong>",css.table="width:100%;")
  #output
  
  #})
  
  #################### interactive plot start
  #################### interactive plot start
  
  check_interactive <- reactiveValues(count = 1)
  search_interactive <- reactiveValues(count = 0)
  id_interactive <- NULL
  #update1 <- reactiveValues(count = 0)
  off_interactive <- reactiveVal(0)
  observeEvent(input$input_file_interactive,{check_interactive$count=1})
  observeEvent(input$feat_inte_interactive,{check_interactive$count=1})
  observeEvent(input$classlabel_inte_interactive,{check_interactive$count=1})
  
  observeEvent(input$start1_interactive,
               {
                 if(isolate(input$choose_graph)=='Manhattan Plot only' || isolate(input$choose_graph)=='Volcano Plot only'){
                   output$nText_interactive <- renderText({shiny::validate(
                     need(isolate(input$input_file_interactive), "No input file provided. Please upload your input file."),
                     need(isolate(input$input_file_interactive$type)=="text/csv" || isolate(input$input_file_interactive$type)=="text/plain", "The format of input file is not correct. Please upload the file with correct format.")
                   )})
                   shiny::validate(
                     need(isolate(input$input_file_interactive), "No input file provided. Please upload your input file."),
                     need(isolate(input$input_file_interactive$type)=="text/csv" || isolate(input$input_file_interactive$type)=="text/plain", "The format of input file is not correct. Please upload the file with correct format.")
                   )
                 }
                 if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot' || isolate(input$choose_graph)=='Volcano Plot with Box Plot'){
                   output$nText_interactive <- renderText({shiny::validate(
                     need(isolate(input$feat_inte_interactive), "No feature table provided. Please upload your feature table."),
                     need(isolate(input$feat_inte_interactive$type)=="text/csv" || isolate(input$feat_inte_interactive$type)=="text/plain", "The format of feature table is not correct. Please upload the file with correct format."),
                     need(isolate(input$classlabel_inte_interactive), "No class label file provided. Please upload class label file."),
                     need(isolate(input$classlabel_inte_interactive$type)=="text/csv" || isolate(input$classlabel_inte_interactive$type)=="text/plain", "The format of class label file is not correct. Please upload the file with correct format.")
                   )})
                   shiny::validate(
                     need(isolate(input$feat_inte_interactive), "No feature table provided. Please upload your feature table."),
                     need(isolate(input$feat_inte_interactive$type)=="text/csv" || isolate(input$feat_inte_interactive$type)=="text/plain", "The format of feature table is not correct. Please upload the file with correct format."),
                     need(isolate(input$classlabel_inte_interactive), "No class label file provided. Please upload class label file."),
                     need(isolate(input$classlabel_inte_interactive$type)=="text/csv" || isolate(input$classlabel_inte_interactive$type)=="text/plain", "The format of class label file is not correct. Please upload the file with correct format.")
                   )
                 }
                 check_interactive$count=0
                 id_interactive <<- showNotification("Figure is being generated now.", duration=NULL)
               })
  
  input_file_interactive <- reactive({
    if(input$start1_interactive!=0  & check_interactive$count==0 & !is.null(isolate(input$input_file_interactive$name)) ){
      
      if((isolate(input$input_file_interactive$type)=="text/csv" || isolate(input$input_file_interactive$type)=="text/plain")){
        req(isolate(input$input_file_interactive))
        if(isolate(input$input_file_interactive$type)=="text/plain"){
          input_file <- read.delim(isolate(input$input_file_interactive$datapath),header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          input_file <- read.csv(isolate(input$input_file_interactive$datapath),header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
        }
        input_file
        
      }
    }
  })
  
  feat_inte_interactive <- reactive({
    if(isolate(input$start1_interactive)!=0  & check_interactive$count==0 & !is.null(isolate(input$feat_inte_interactive$name)) ){
      
      if((isolate(input$feat_inte_interactive$type)=="text/csv" || isolate(input$feat_inte_interactive$type)=="text/plain")){
        req(isolate(input$feat_inte_interactive))
        if(isolate(input$feat_inte_interactive$type)=="text/plain"){
          feature_table_file <- read.delim(isolate(input$feat_inte_interactive$datapath),header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          feature_table_file <- read.csv(isolate(input$feat_inte_interactive$datapath),header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
        }
        feature_table_file
        
      }
    }
  })
  
  classlabel_inte_interactive <- reactive({
    if(isolate(input$start1_interactive)!=0  & check_interactive$count==0 & !is.null(isolate(input$classlabel_inte_interactive$name))){
      
      if((isolate(input$classlabel_inte_interactive$type)=="text/csv" || isolate(input$classlabel_inte_interactive$type=="text/plain"))){
        req(isolate(input$classlabel_inte_interactive))
        if(isolate(input$classlabel_inte_interactive$type)=="text/plain"){
          class_labels_file <- read.delim(isolate(input$classlabel_inte_interactive$datapath),header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          class_labels_file <- read.csv(isolate(input$classlabel_inte_interactive$datapath),header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
        }
        class_labels_file
      }
    }
  })
  
  
  all_alert <- reactive({
    
    if(!is.null(input_file_interactive())){
      
      if(isolate(input$choose_graph)=='Manhattan Plot only'){
        
        if(isolate(input$yaxislabel_manhattan_only)=='pvalue' && isolate(input$adjdashline_manhattan_only)=='no'){
          
          if(ncol(input_file_interactive())==3 || ncol(input_file_interactive())==4){
            closeAlert(session, "check_input_interactive_single_alert")
            check_input_interactive_single <- TRUE
          }else{
            closeAlert(session, "check_input_interactive_single_alert")
            createAlert(session, "alert_interactive", "check_input_interactive_single_alert", title = "Input File Error", content = "The input file is only allowed to have 3 or 4 columns separately for 'Name', 'log2foldchange', 'p-value', 'adjusted p-value'(optional).", append = TRUE)
            check_input_interactive_single <- FALSE
          }
          
        } else if(isolate(input$yaxislabel_manhattan_only)=='pvalue' && isolate(input$adjdashline_manhattan_only)=='yes') {
          
          if(ncol(input_file_interactive())==4){
            closeAlert(session, "check_input_interactive_single_alert")
            check_input_interactive_single <- TRUE
          }else{
            closeAlert(session, "check_input_interactive_single_alert")
            createAlert(session, "alert_interactive", "check_input_interactive_single_alert", title = "Input File Error", content = "The input file is only allowed to have 4 columns separately for 'Name', 'log2foldchange', 'p-value', 'adjusted p-value'.", append = TRUE)
            check_input_interactive_single <- FALSE
          }
          
        } else{
          
          if(ncol(input_file_interactive())==3){
            closeAlert(session, "check_input_interactive_single_alert")
            check_input_interactive_single <- TRUE
          }else{
            closeAlert(session, "check_input_interactive_single_alert")
            createAlert(session, "alert_interactive", "check_input_interactive_single_alert", title = "Input File Error", content = "The input file is only allowed to have 3 columns separately for 'Name', 'log2foldchange', 'VIP'.", append = TRUE)
            check_input_interactive_single <- FALSE
          }
        }
        
      }else{
        
        if(isolate(input$adjdashline_volcano_only)=='no'){
          
          if(ncol(input_file_interactive())==3 || ncol(input_file_interactive())==4){
            closeAlert(session, "check_input_interactive_single_alert")
            check_input_interactive_single <- TRUE
          }else{
            closeAlert(session, "check_input_interactive_single_alert")
            createAlert(session, "alert_interactive", "check_input_interactive_single_alert", title = "Input File Error", content = "The input file is only allowed to have 3 or 4 columns separately for 'Name', 'log2foldchange', 'p-value', 'adjusted p-value'(optional).", append = TRUE)
            check_input_interactive_single <- FALSE
          }
          
        }else{
          
          if(ncol(input_file_interactive())==4){
            closeAlert(session, "check_input_interactive_single_alert")
            check_input_interactive_single <- TRUE
          }else{
            closeAlert(session, "check_input_interactive_single_alert")
            createAlert(session, "alert_interactive", "check_input_interactive_single_alert", title = "Input File Error", content = "The input file is only allowed to have 4 columns separately for 'Name', 'log2foldchange', 'p-value', 'adjusted p-value'.", append = TRUE)
            check_input_interactive_single <- FALSE
          }
          
        }
        
      }
      
    }else{
      
      check_input_interactive_single <- TRUE
    }
    
    if(is.na(isolate(input$x_axis_spacing_type1_manhattan_only))){
      closeAlert(session, "x_axis_spacing_type1_manhattan_only_alert")
      createAlert(session, "alert_interactive", "x_axis_spacing_type1_manhattan_only_alert", title = "Argument Input Error", content = "'Type1 X axis spacing' argument can't be empty.", append = TRUE)
      check_x_axis_spacing_type1_manhattan_only <- FALSE
    } else if (isolate(input$x_axis_spacing_type1_manhattan_only)<0 | isolate(input$x_axis_spacing_type1_manhattan_only)>1000) {
      closeAlert(session, "x_axis_spacing_type1_manhattan_only_alert")
      createAlert(session, "alert_interactive", "x_axis_spacing_type1_manhattan_only_alert", title = "Argument Input Error", content = "'Type1 X axis spacing' argument should be not smaller than 0 or larger than 1000.", append = TRUE)
      check_x_axis_spacing_type1_manhattan_only <- FALSE
    } else {
      closeAlert(session, "x_axis_spacing_type1_manhattan_only_alert")
      check_x_axis_spacing_type1_manhattan_only <- TRUE
    }
    
    if(is.na(isolate(input$x_axis_spacing_type2_manhattan_only))){
      closeAlert(session, "x_axis_spacing_type2_manhattan_only_alert")
      createAlert(session, "alert_interactive", "x_axis_spacing_type2_manhattan_only_alert", title = "Argument Input Error", content = "'Type2 X axis spacing' argument can't be empty.", append = TRUE)
      check_x_axis_spacing_type2_manhattan_only <- FALSE
    } else if (isolate(input$x_axis_spacing_type2_manhattan_only)<0 | isolate(input$x_axis_spacing_type2_manhattan_only)>1000) {
      closeAlert(session, "x_axis_spacing_type2_manhattan_only_alert")
      createAlert(session, "alert_interactive", "x_axis_spacing_type2_manhattan_only_alert", title = "Argument Input Error", content = "'Type2 X axis spacing' argument should be not smaller than 0 or larger than 1000.", append = TRUE)
      check_x_axis_spacing_type2_manhattan_only <- FALSE
    } else {
      closeAlert(session, "x_axis_spacing_type2_manhattan_only_alert")
      check_x_axis_spacing_type2_manhattan_only <- TRUE
    }
    
    if(is.na(isolate(input$y_axis_spacing_manhattan_only))){
      closeAlert(session, "y_axis_spacing_manhattan_only_alert")
      createAlert(session, "alert_interactive", "y_axis_spacing_manhattan_only_alert", title = "Argument Input Error", content = "'Y axis spacing' argument can't be empty.", append = TRUE)
      check_y_axis_spacing_manhattan_only <- FALSE
    } else if (isolate(input$y_axis_spacing_manhattan_only)<0 | isolate(input$y_axis_spacing_manhattan_only)>5) {
      closeAlert(session, "y_axis_spacing_manhattan_only_alert")
      createAlert(session, "alert_interactive", "y_axis_spacing_manhattan_only_alert", title = "Argument Input Error", content = "'Y axis spacing' argument should be not smaller than 0 or larger than 5.", append = TRUE)
      check_y_axis_spacing_manhattan_only <- FALSE
    } else {
      closeAlert(session, "y_axis_spacing_manhattan_only_alert")
      check_y_axis_spacing_manhattan_only <- TRUE
    }
    
    if(is.na(isolate(input$pvaluecutoff_manhattan_only))){
      closeAlert(session, "pvaluecutoff_manhattan_only_alert")
      createAlert(session, "alert_interactive", "pvaluecutoff_manhattan_only_alert", title = "Argument Input Error", content = "'Threshold for p-value' argument can't be empty.", append = TRUE)
      check_pvaluecutoff_manhattan_only <- FALSE
    } else if (isolate(input$pvaluecutoff_manhattan_only)<0 | isolate(input$pvaluecutoff_manhattan_only)>1) {
      closeAlert(session, "pvaluecutoff_manhattan_only_alert")
      createAlert(session, "alert_interactive", "pvaluecutoff_manhattan_only_alert", title = "Argument Input Error", content = "'Threshold for p-value' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      check_pvaluecutoff_manhattan_only <- FALSE
    } else {
      closeAlert(session, "pvaluecutoff_manhattan_only_alert")
      check_pvaluecutoff_manhattan_only <- TRUE
    }
    
    if(is.na(isolate(input$adjpvaluecutoff_manhattan_only))){
      closeAlert(session, "adjpvaluecutoff_manhattan_only_alert")
      createAlert(session, "alert_interactive", "adjpvaluecutoff_manhattan_only_alert", title = "Argument Input Error", content = "'Threshold for adjusted p-value' argument can't be empty.", append = TRUE)
      check_adjpvaluecutoff_manhattan_only <- FALSE
    } else if (isolate(input$adjpvaluecutoff_manhattan_only)<0 | isolate(input$adjpvaluecutoff_manhattan_only)>1) {
      closeAlert(session, "adjpvaluecutoff_manhattan_only_alert")
      createAlert(session, "alert_interactive", "adjpvaluecutoff_manhattan_only_alert", title = "Argument Input Error", content = "'Threshold for adjusted p-value' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      check_adjpvaluecutoff_manhattan_only <- FALSE
    } else {
      closeAlert(session, "adjpvaluecutoff_manhattan_only_alert")
      check_adjpvaluecutoff_manhattan_only <- TRUE
    }
    
    if(is.na(isolate(input$vipcutoff_manhattan_only))){
      closeAlert(session, "vipcutoff_manhattan_only_alert")
      createAlert(session, "alert_interactive", "vipcutoff_manhattan_only_alert", title = "Argument Input Error", content = "'Threshold for VIP' argument can't be empty.", append = TRUE)
      check_vipcutoff_manhattan_only <- FALSE
    } else if (isolate(input$vipcutoff_manhattan_only)<0 | isolate(input$vipcutoff_manhattan_only)>10) {
      closeAlert(session, "vipcutoff_manhattan_only_alert")
      createAlert(session, "alert_interactive", "vipcutoff_manhattan_only_alert", title = "Argument Input Error", content = "'Threshold for VIP' argument should be not smaller than 0 or larger than 10.", append = TRUE)
      check_vipcutoff_manhattan_only <- FALSE
    } else {
      closeAlert(session, "vipcutoff_manhattan_only_alert")
      check_vipcutoff_manhattan_only <- TRUE
    }
    
    if(is.na(isolate(input$pvaluecutoff_volcano_only))){
      closeAlert(session, "pvaluecutoff_volcano_only_alert")
      createAlert(session, "alert_interactive", "pvaluecutoff_volcano_only_alert", title = "Argument Input Error", content = "'Threshold for p-value' argument can't be empty.", append = TRUE)
      check_pvaluecutoff_volcano_only <- FALSE
    } else if (isolate(input$pvaluecutoff_volcano_only)<0 | isolate(input$pvaluecutoff_volcano_only)>1) {
      closeAlert(session, "pvaluecutoff_volcano_only_alert")
      createAlert(session, "alert_interactive", "pvaluecutoff_volcano_only_alert", title = "Argument Input Error", content = "'Threshold for p-value' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      check_pvaluecutoff_volcano_only <- FALSE
    } else {
      closeAlert(session, "pvaluecutoff_volcano_only_alert")
      check_pvaluecutoff_volcano_only <- TRUE
    }
    
    if(is.na(isolate(input$lfc_volcano_only))){
      closeAlert(session, "lfc_volcano_only_alert")
      createAlert(session, "alert_interactive", "lfc_volcano_only_alert", title = "Argument Input Error", content = "'Left side threshold for fold change' argument can't be empty.", append = TRUE)
      check_lfc_volcano_only <- FALSE
    } else if (isolate(input$lfc_volcano_only)<0 | isolate(input$lfc_volcano_only)>1) {
      closeAlert(session, "lfc_volcano_only_alert")
      createAlert(session, "alert_interactive", "lfc_volcano_only_alert", title = "Argument Input Error", content = "'Left side threshold for fold change' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      check_lfc_volcano_only <- FALSE
    } else {
      closeAlert(session, "lfc_volcano_only_alert")
      check_lfc_volcano_only <- TRUE
    }
    
    if(is.na(isolate(input$rfc_volcano_only))){
      closeAlert(session, "rfc_volcano_only_alert")
      createAlert(session, "alert_interactive", "rfc_volcano_only_alert", title = "Argument Input Error", content = "'Right side threshold for fold change' argument can't be empty.", append = TRUE)
      check_rfc_volcano_only <- FALSE
    } else if (isolate(input$rfc_volcano_only)<1 | isolate(input$rfc_volcano_only)>10) {
      closeAlert(session, "rfc_volcano_only_alert")
      createAlert(session, "alert_interactive", "rfc_volcano_only_alert", title = "Argument Input Error", content = "'Right side threshold for fold change' argument should be not smaller than 1 or larger than 10.", append = TRUE)
      check_rfc_volcano_only <- FALSE
    } else {
      closeAlert(session, "rfc_volcano_only_alert")
      check_rfc_volcano_only <- TRUE
    }
    
    if(is.na(isolate(input$x_axis_boundary_volcano_only))){
      closeAlert(session, "x_axis_boundary_volcano_only_alert")
      createAlert(session, "alert_interactive", "x_axis_boundary_volcano_only_alert", title = "Argument Input Error", content = "'X axis boundary' argument can't be empty.", append = TRUE)
      check_x_axis_boundary_volcano_only <- FALSE
    } else if (isolate(input$x_axis_boundary_volcano_only)<0 | isolate(input$x_axis_boundary_volcano_only)>20) {
      closeAlert(session, "x_axis_boundary_volcano_only_alert")
      createAlert(session, "alert_interactive", "x_axis_boundary_volcano_only_alert", title = "Argument Input Error", content = "'X axis boundary' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      check_x_axis_boundary_volcano_only <- FALSE
    } else {
      closeAlert(session, "x_axis_boundary_volcano_only_alert")
      check_x_axis_boundary_volcano_only <- TRUE
    }
    
    if(is.na(isolate(input$y_axis_spacing_volcano_only))){
      closeAlert(session, "y_axis_spacing_volcano_only_alert")
      createAlert(session, "alert_interactive", "y_axis_spacing_volcano_only_alert", title = "Argument Input Error", content = "'Y axis spacing' argument can't be empty.", append = TRUE)
      check_y_axis_spacing_volcano_only <- FALSE
    } else if (isolate(input$y_axis_spacing_volcano_only)<0 | isolate(input$y_axis_spacing_volcano_only)>5) {
      closeAlert(session, "y_axis_spacing_volcano_only_alert")
      createAlert(session, "alert_interactive", "y_axis_spacing_volcano_only_alert", title = "Argument Input Error", content = "'Y axis spacing' argument should be not smaller than 0 or larger than 5.", append = TRUE)
      check_y_axis_spacing_volcano_only <- FALSE
    } else {
      closeAlert(session, "y_axis_spacing_volcano_only_alert")
      check_y_axis_spacing_volcano_only <- TRUE
    }
    
    if(is.na(isolate(input$adjpvaluecutoff_volcano_only))){
      closeAlert(session, "adjpvaluecutoff_volcano_only_alert")
      createAlert(session, "alert_interactive", "adjpvaluecutoff_volcano_only_alert", title = "Argument Input Error", content = "'Threshold for adjusted p-value' argument can't be empty.", append = TRUE)
      check_adjpvaluecutoff_volcano_only <- FALSE
    } else if (isolate(input$adjpvaluecutoff_volcano_only)<0 | isolate(input$adjpvaluecutoff_volcano_only)>1) {
      closeAlert(session, "adjpvaluecutoff_volcano_only_alert")
      createAlert(session, "alert_interactive", "adjpvaluecutoff_volcano_only_alert", title = "Argument Input Error", content = "'Threshold for adjusted p-value' argument should be not smaller than 0 or larger than 1.", append = TRUE)
      check_adjpvaluecutoff_volcano_only <- FALSE
    } else {
      closeAlert(session, "adjpvaluecutoff_volcano_only_alert")
      check_adjpvaluecutoff_volcano_only <- TRUE
    }
    
    if(is.na(isolate(input$x_axis_spacing_type1_manhattan_box))){
      closeAlert(session, "x_axis_spacing_type1_manhattan_box_alert")
      createAlert(session, "alert_interactive", "x_axis_spacing_type1_manhattan_box_alert", title = "Argument Input Error", content = "'Type1 X axis spacing' argument can't be empty.", append = TRUE)
      check_x_axis_spacing_type1_manhattan_box <- FALSE
    } else if (isolate(input$x_axis_spacing_type1_manhattan_box)<0 | isolate(input$x_axis_spacing_type1_manhattan_box)>1000) {
      closeAlert(session, "x_axis_spacing_type1_manhattan_box_alert")
      createAlert(session, "alert_interactive", "x_axis_spacing_type1_manhattan_box_alert", title = "Argument Input Error", content = "'Type1 X axis spacing' argument should be not smaller than 0 or larger than 1000.", append = TRUE)
      check_x_axis_spacing_type1_manhattan_box <- FALSE
    } else {
      closeAlert(session, "x_axis_spacing_type1_manhattan_box_alert")
      check_x_axis_spacing_type1_manhattan_box <- TRUE
    }
    
    if(is.na(isolate(input$x_axis_spacing_type2_manhattan_box))){
      closeAlert(session, "x_axis_spacing_type2_manhattan_box_alert")
      createAlert(session, "alert_interactive", "x_axis_spacing_type2_manhattan_box_alert", title = "Argument Input Error", content = "'Type2 X axis spacing' argument can't be empty.", append = TRUE)
      check_x_axis_spacing_type2_manhattan_box <- FALSE
    } else if (isolate(input$x_axis_spacing_type2_manhattan_box)<0 | isolate(input$x_axis_spacing_type2_manhattan_box)>1000) {
      closeAlert(session, "x_axis_spacing_type2_manhattan_box_alert")
      createAlert(session, "alert_interactive", "x_axis_spacing_type2_manhattan_box_alert", title = "Argument Input Error", content = "'Type2 X axis spacing' argument should be not smaller than 0 or larger than 1000.", append = TRUE)
      check_x_axis_spacing_type2_manhattan_box <- FALSE
    } else {
      closeAlert(session, "x_axis_spacing_type2_manhattan_box_alert")
      check_x_axis_spacing_type2_manhattan_box <- TRUE
    }
    
    if(is.na(isolate(input$y_axis_spacing_manhattan_box))){
      closeAlert(session, "y_axis_spacing_manhattan_box_alert")
      createAlert(session, "alert_interactive", "y_axis_spacing_manhattan_box_alert", title = "Argument Input Error", content = "'Y axis spacing' argument can't be empty.", append = TRUE)
      check_y_axis_spacing_manhattan_box <- FALSE
    } else if (isolate(input$y_axis_spacing_manhattan_box)<0 | isolate(input$y_axis_spacing_manhattan_box)>5) {
      closeAlert(session, "y_axis_spacing_manhattan_box_alert")
      createAlert(session, "alert_interactive", "y_axis_spacing_manhattan_box_alert", title = "Argument Input Error", content = "'Y axis spacing' argument should be not smaller than 0 or larger than 5", append = TRUE)
      check_y_axis_spacing_manhattan_box <- FALSE
    } else {
      closeAlert(session, "y_axis_spacing_manhattan_box_alert")
      check_y_axis_spacing_manhattan_box <- TRUE
    }
    
    if(is.na(isolate(input$pvaluecutoff_manhattan_box))){
      closeAlert(session, "pvaluecutoff_manhattan_box_alert")
      createAlert(session, "alert_interactive", "pvaluecutoff_manhattan_box_alert", title = "Argument Input Error", content = "'Y axis spacing' argument can't be empty.", append = TRUE)
      check_pvaluecutoff_manhattan_box <- FALSE
    } else if (isolate(input$pvaluecutoff_manhattan_box)<0 | isolate(input$pvaluecutoff_manhattan_box)>1) {
      closeAlert(session, "pvaluecutoff_manhattan_box_alert")
      createAlert(session, "alert_interactive", "pvaluecutoff_manhattan_box_alert", title = "Argument Input Error", content = "'Y axis spacing' argument should be not smaller than 0 or larger than 1", append = TRUE)
      check_pvaluecutoff_manhattan_box <- FALSE
    } else {
      closeAlert(session, "pvaluecutoff_manhattan_box_alert")
      check_pvaluecutoff_manhattan_box <- TRUE
    }
    
    if(is.na(isolate(input$adjpvaluecutoff_manhattan_box))){
      closeAlert(session, "adjpvaluecutoff_manhattan_box_alert")
      createAlert(session, "alert_interactive", "adjpvaluecutoff_manhattan_box_alert", title = "Argument Input Error", content = "'Threshold for adjusted p-value' argument can't be empty.", append = TRUE)
      check_adjpvaluecutoff_manhattan_box <- FALSE
    } else if (isolate(input$adjpvaluecutoff_manhattan_box)<0 | isolate(input$adjpvaluecutoff_manhattan_box)>1) {
      closeAlert(session, "adjpvaluecutoff_manhattan_box_alert")
      createAlert(session, "alert_interactive", "adjpvaluecutoff_manhattan_box_alert", title = "Argument Input Error", content = "'Threshold for adjusted p-value' argument should be not smaller than 0 or larger than 1", append = TRUE)
      check_adjpvaluecutoff_manhattan_box <- FALSE
    } else {
      closeAlert(session, "adjpvaluecutoff_manhattan_box_alert")
      check_adjpvaluecutoff_manhattan_box <- TRUE
    }
    
    if(is.na(isolate(input$vipcutoff_manhattan_box))){
      closeAlert(session, "vipcutoff_manhattan_box_alert")
      createAlert(session, "alert_interactive", "vipcutoff_manhattan_box_alert", title = "Argument Input Error", content = "'Threshold for VIP' argument can't be empty.", append = TRUE)
      check_vipcutoff_manhattan_box <- FALSE
    } else if (isolate(input$vipcutoff_manhattan_box)<0 | isolate(input$vipcutoff_manhattan_box)>10) {
      closeAlert(session, "vipcutoff_manhattan_box_alert")
      createAlert(session, "alert_interactive", "vipcutoff_manhattan_box_alert", title = "Argument Input Error", content = "'Threshold for VIP' argument should be not smaller than 0 or larger than 10", append = TRUE)
      check_vipcutoff_manhattan_box <- FALSE
    } else {
      closeAlert(session, "vipcutoff_manhattan_box_alert")
      check_vipcutoff_manhattan_box <- TRUE
    }
    
    if(is.na(isolate(input$pvaluecutoff_volcano_box))){
      closeAlert(session, "pvaluecutoff_volcano_box_alert")
      createAlert(session, "alert_interactive", "pvaluecutoff_volcano_box_alert", title = "Argument Input Error", content = "'Threshold for p-value' argument can't be empty.", append = TRUE)
      check_pvaluecutoff_volcano_box <- FALSE
    } else if (isolate(input$pvaluecutoff_volcano_box)<0 | isolate(input$pvaluecutoff_volcano_box)>1) {
      closeAlert(session, "pvaluecutoff_volcano_box_alert")
      createAlert(session, "alert_interactive", "pvaluecutoff_volcano_box_alert", title = "Argument Input Error", content = "'Threshold for p-value' argument should be not smaller than 0 or larger than 1", append = TRUE)
      check_pvaluecutoff_volcano_box <- FALSE
    } else {
      closeAlert(session, "pvaluecutoff_volcano_box_alert")
      check_pvaluecutoff_volcano_box <- TRUE
    }
    
    if(is.na(isolate(input$lfc_volcano_box))){
      closeAlert(session, "lfc_volcano_box_alert")
      createAlert(session, "alert_interactive", "lfc_volcano_box_alert", title = "Argument Input Error", content = "'Left side threshold for fold change' argument can't be empty.", append = TRUE)
      check_lfc_volcano_box <- FALSE
    } else if (isolate(input$lfc_volcano_box)<0 | isolate(input$lfc_volcano_box)>1) {
      closeAlert(session, "lfc_volcano_box_alert")
      createAlert(session, "alert_interactive", "lfc_volcano_box_alert", title = "Argument Input Error", content = "'Left side threshold for fold change' argument should be not smaller than 0 or larger than 1", append = TRUE)
      check_lfc_volcano_box <- FALSE
    } else {
      closeAlert(session, "lfc_volcano_box_alert")
      check_lfc_volcano_box <- TRUE
    }
    
    if(is.na(isolate(input$rfc_volcano_box))){
      closeAlert(session, "rfc_volcano_box_alert")
      createAlert(session, "alert_interactive", "rfc_volcano_box_alert", title = "Argument Input Error", content = "'Right side threshold for fold change' argument can't be empty.", append = TRUE)
      check_rfc_volcano_box <- FALSE
    } else if (isolate(input$rfc_volcano_box)<1 | isolate(input$rfc_volcano_box)>10) {
      closeAlert(session, "rfc_volcano_box_alert")
      createAlert(session, "alert_interactive", "rfc_volcano_box_alert", title = "Argument Input Error", content = "'Right side threshold for fold change' argument should be not smaller than 1 or larger than 10", append = TRUE)
      check_rfc_volcano_box <- FALSE
    } else {
      closeAlert(session, "rfc_volcano_box_alert")
      check_rfc_volcano_box <- TRUE
    }
    
    if(is.na(isolate(input$x_axis_boundary_volcano_box))){
      closeAlert(session, "x_axis_boundary_volcano_box_alert")
      createAlert(session, "alert_interactive", "x_axis_boundary_volcano_box_alert", title = "Argument Input Error", content = "'X axis boundary' argument can't be empty.", append = TRUE)
      check_x_axis_boundary_volcano_box <- FALSE
    } else if (isolate(input$x_axis_boundary_volcano_box)<0 | isolate(input$x_axis_boundary_volcano_box)>20) {
      closeAlert(session, "x_axis_boundary_volcano_box_alert")
      createAlert(session, "alert_interactive", "x_axis_boundary_volcano_box_alert", title = "Argument Input Error", content = "'X axis boundary' argument should be not smaller than 0 or larger than 20", append = TRUE)
      check_x_axis_boundary_volcano_box <- FALSE
    } else {
      closeAlert(session, "x_axis_boundary_volcano_box_alert")
      check_x_axis_boundary_volcano_box <- TRUE
    }
    
    if(is.na(isolate(input$y_axis_spacing_volcano_box))){
      closeAlert(session, "y_axis_spacing_volcano_box_alert")
      createAlert(session, "alert_interactive", "y_axis_spacing_volcano_box_alert", title = "Argument Input Error", content = "'Y axis spacing' argument can't be empty.", append = TRUE)
      check_y_axis_spacing_volcano_box <- FALSE
    } else if (isolate(input$y_axis_spacing_volcano_box)<0 | isolate(input$y_axis_spacing_volcano_box)>5) {
      closeAlert(session, "y_axis_spacing_volcano_box_alert")
      createAlert(session, "alert_interactive", "y_axis_spacing_volcano_box_alert", title = "Argument Input Error", content = "'Y axis spacing' argument should be not smaller than 0 or larger than 5", append = TRUE)
      check_y_axis_spacing_volcano_box <- FALSE
    } else {
      closeAlert(session, "y_axis_spacing_volcano_box_alert")
      check_y_axis_spacing_volcano_box <- TRUE
    }
    
    if(is.na(isolate(input$adjpvaluecutoff_volcano_box))){
      closeAlert(session, "adjpvaluecutoff_volcano_box_alert")
      createAlert(session, "alert_interactive", "adjpvaluecutoff_volcano_box_alert", title = "Argument Input Error", content = "'Threshold for adjusted p-value' argument can't be empty.", append = TRUE)
      check_adjpvaluecutoff_volcano_box <- FALSE
    } else if (isolate(input$adjpvaluecutoff_volcano_box)<0 | isolate(input$adjpvaluecutoff_volcano_box)>1) {
      closeAlert(session, "adjpvaluecutoff_volcano_box_alert")
      createAlert(session, "alert_interactive", "adjpvaluecutoff_volcano_box_alert", title = "Argument Input Error", content = "'Threshold for adjusted p-value' argument should be not smaller than 0 or larger than 1", append = TRUE)
      check_adjpvaluecutoff_volcano_box <- FALSE
    } else {
      closeAlert(session, "adjpvaluecutoff_volcano_box_alert")
      check_adjpvaluecutoff_volcano_box <- TRUE
    }
    
    all(check_input_interactive_single,check_x_axis_spacing_type1_manhattan_only,check_x_axis_spacing_type2_manhattan_only,check_y_axis_spacing_manhattan_only,check_pvaluecutoff_manhattan_only,check_adjpvaluecutoff_manhattan_only,check_vipcutoff_manhattan_only,
        check_pvaluecutoff_volcano_only,check_lfc_volcano_only,check_rfc_volcano_only,check_x_axis_boundary_volcano_only,check_y_axis_spacing_volcano_only,check_adjpvaluecutoff_volcano_only,
        check_x_axis_spacing_type1_manhattan_box,check_x_axis_spacing_type2_manhattan_box,check_y_axis_spacing_manhattan_box,check_pvaluecutoff_manhattan_box,check_adjpvaluecutoff_manhattan_box,check_vipcutoff_manhattan_box,
        check_pvaluecutoff_volcano_box,check_lfc_volcano_box,check_rfc_volcano_box,check_x_axis_boundary_volcano_box,check_y_axis_spacing_volcano_box,check_adjpvaluecutoff_volcano_box)
    
  })
  
  
  plot1_single_plot <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot only' || isolate(input$choose_graph)=='Volcano Plot only'){
      
      if(input$start1_interactive!=0  & check_interactive$count==0 & !is.null(input_file_interactive()) & all_alert()){
        
        if(isolate(input$choose_graph)=='Manhattan Plot only'){
          
          if(isolate(input$yaxislabel_manhattan_only)=='pvalue'){
            
            if(ncol(input_file_interactive())==3){
              manhattan_only = input_file_interactive()
              colnames(manhattan_only) = c('name','log2foldchange','p.value')
            }else{
              manhattan_only = input_file_interactive()
              colnames(manhattan_only) = c('name','log2foldchange','p.value','adjusted.p.value')
            }
            
            manhattan_only = manhattan_only[grep('_',manhattan_only$name),]
            manhattan_only = manhattan_only %>% separate('name',c("mz","time"),remove=FALSE,sep="_")
            manhattan_only$mz=as.numeric(manhattan_only$mz)
            manhattan_only$time=as.numeric(manhattan_only$time)
            
            if(isolate(input$adjdashline_manhattan_only)=='yes'){
              
              if(isolate(input$psignif_manhattan_only)=='pvalue'){
                
                cutoff_manhattan_only <- isolate(input$pvaluecutoff_manhattan_only)
              }else{
                
                cutoff_manhattan_only <- max(manhattan_only[manhattan_only$p.value<isolate(input$pvaluecutoff_manhattan_only) & manhattan_only$adjusted.p.value<isolate(input$adjpvaluecutoff_manhattan_only),'p.value'])
              }
              
              adjdashline_manhattan_only <- geom_hline(yintercept=-log10(max(manhattan_only[manhattan_only$p.value<isolate(input$pvaluecutoff_manhattan_only) & manhattan_only$adjusted.p.value<isolate(input$adjpvaluecutoff_manhattan_only),'p.value'])),linetype='dotted',color= input$dottedlinecol_manhattan_only)
            }else{
              
              cutoff_manhattan_only <- isolate(input$pvaluecutoff_manhattan_only)
              adjdashline_manhattan_only <- NULL
            }
            
            
            if(isolate(input$labelexpression_manhattan_only)=='yes'){
              
              manhattan_only[manhattan_only$p.value<cutoff_manhattan_only & manhattan_only$log2foldchange>0,'color']<- input$poscol_manhattan_only
              manhattan_only[manhattan_only$p.value<cutoff_manhattan_only & manhattan_only$log2foldchange<0,'color']<- input$negcol_manhattan_only
            }else{
              
              manhattan_only[manhattan_only$p.value<cutoff_manhattan_only,'color']<- input$sigcol_manhattan_only
            }
            
            manhattan_only[manhattan_only$p.value>=cutoff_manhattan_only,'color']<- input$insigcol_manhattan_only
            manhattan_only[manhattan_only$color== input$insigcol_manhattan_only,'size']=0.2
            manhattan_only[!(manhattan_only$color== input$insigcol_manhattan_only),'size']=1.5
            pal_manhattan_only <- levels(factor(manhattan_only$color))
            names(pal_manhattan_only)<-pal_manhattan_only
            
            xincrement_manhattan_only = isolate(input$x_axis_spacing_type1_manhattan_only)
            xmin_val_manhattan_only <- min(0,manhattan_only$mz)
            xmax_val_manhattan_only <- max(manhattan_only$mz)
            xlabel_manhattan_only <- xlab('mass-to-charge (m/z)')
            maintitle_manhattan_only <- ggtitle('Type 1 manhattan plot (-log10p vs mz)')
            gplotline_manhattan_only <- ggplot(manhattan_only, aes(x=mz, y=-log10(p.value), color=color, text=paste0('mz:',manhattan_only$mz," time:",manhattan_only$time)))
            
            yvec_manhattan_only <- -log10(manhattan_only$p.value)
            if(max(yvec_manhattan_only)>20){yvec_manhattan_only[yvec_manhattan_only>20]<-20}
            ymax_val_manhattan_only <- max(yvec_manhattan_only)
            ymax_manhattan_only =floor(max(yvec_manhattan_only))+1
            ylabel_manhattan_only <- ylab('-log10p')
            
            dashline_manhattan_only <- geom_hline(yintercept=-log10(isolate(input$pvaluecutoff_manhattan_only)),linetype='dashed',color= input$dashedlinecol_manhattan_only)
            
          }else{
            
            manhattan_only = input_file_interactive()
            colnames(manhattan_only) = c('name','log2foldchange','VIP')
            manhattan_only=manhattan_only[grep("_",manhattan_only$name),]
            rownames(manhattan_only)=1:nrow(manhattan_only)
            manhattan_only = manhattan_only %>% separate(name,c("mz","time"),remove=FALSE,sep="_")
            manhattan_only$mz=as.numeric(manhattan_only$mz)
            manhattan_only$time=as.numeric(manhattan_only$time)
            
            if(isolate(input$labelexpression_manhattan_only)=='yes'){
              
              manhattan_only[manhattan_only$VIP > isolate(input$vipcutoff_manhattan_only) & manhattan_only$log2foldchange>0,'color']<- input$poscol_manhattan_only
              manhattan_only[manhattan_only$VIP > isolate(input$vipcutoff_manhattan_only) & manhattan_only$log2foldchange<0,'color']<- input$negcol_manhattan_only
            }else{
              
              manhattan_only[manhattan_only$VIP > isolate(input$vipcutoff_manhattan_only),'color']<- input$sigcol_manhattan_only
            }
            
            manhattan_only[manhattan_only$VIP<= isolate(input$vipcutoff_manhattan_only),'color']<- input$insigcol_manhattan_only
            manhattan_only[manhattan_only$color== input$insigcol_manhattan_only,'size']=0.2
            manhattan_only[!(manhattan_only$color== input$insigcol_manhattan_only),'size']=1.5
            pal_manhattan_only <- levels(factor(manhattan_only$color))
            names(pal_manhattan_only)<-pal_manhattan_only
            
            xincrement_manhattan_only = isolate(input$x_axis_spacing_type1_manhattan_only)
            xmin_val_manhattan_only <- min(0,manhattan_only$mz)
            xmax_val_manhattan_only <- max(manhattan_only$mz)
            xlabel_manhattan_only <- xlab('mass-to-charge (m/z)')
            maintitle_manhattan_only <- ggtitle('Type 1 manhattan plot (VIP vs mz)')
            gplotline_manhattan_only <- ggplot(manhattan_only, aes(x=mz, y=VIP, color=color, text=paste0('mz:',manhattan_only$mz," time:",manhattan_only$time)))
            
            
            yvec_manhattan_only <- manhattan_only$VIP
            if(max(yvec_manhattan_only)>20){yvec_manhattan_only[yvec_manhattan_only>20]<-20}
            ymax_val_manhattan_only <- max(yvec_manhattan_only)
            ymax_manhattan_only =floor(max(yvec_manhattan_only))+1
            ylabel_manhattan_only <- ylab('VIP')
            
            dashline_manhattan_only <- geom_hline(yintercept= isolate(input$vipcutoff_manhattan_only),linetype='dashed',color= input$dashedlinecol_manhattan_only)
            
            adjdashline_manhattan_only <- NULL
            
          }
          
          
          p1 <- gplotline_manhattan_only +
            geom_point(show.legend=F,size=manhattan_only$size,alpha=0.7) +
            xlabel_manhattan_only + ylabel_manhattan_only +
            scale_color_manual(values=pal_manhattan_only) +
            scale_x_continuous(breaks = seq(xmin_val_manhattan_only, xmax_val_manhattan_only, by = xincrement_manhattan_only)) +
            scale_y_continuous(breaks = seq(0, (ymax_val_manhattan_only + 2), by = isolate(input$y_axis_spacing_manhattan_only)), limits = c(0,ymax_manhattan_only)) +
            dashline_manhattan_only + adjdashline_manhattan_only +
            maintitle_manhattan_only +
            theme(axis.title=element_text(size=12),
                  axis.text =element_text(size=10),
                  axis.text.x = element_text(angle=0),
                  axis.line = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                  plot.title = element_text(size=14, hjust = 0.5),
                  legend.position = "none")
          
          closeAlert(session, "check_search_interactive_alert")
          p1
          
        }else{
          
          if(ncol(input_file_interactive())==3){
            volcano_only = input_file_interactive()
            colnames(volcano_only) = c('name','log2foldchange','p.value')
          }else{
            volcano_only = input_file_interactive()
            colnames(volcano_only) = c('name','log2foldchange','p.value','adjusted.p.value')
          }
          
          
          if(isolate(input$adjdashline_volcano_only)=='yes'){
            
            if(isolate(input$psignif_volcano_only)=='pvalue'){
              
              cutoff_volcano_only <- isolate(input$pvaluecutoff_volcano_only)
            }else{
              
              cutoff_volcano_only <- max(volcano_only[volcano_only$p.value<isolate(input$pvaluecutoff_volcano_only) & volcano_only$adjusted.p.value<isolate(input$adjpvaluecutoff_volcano_only),'p.value'])
            }
            
            adjdashline_volcano_only <- geom_hline(yintercept=-log10(max(volcano_only[volcano_only$p.value<isolate(input$pvaluecutoff_volcano_only) & volcano_only$adjusted.p.value<isolate(input$adjpvaluecutoff_volcano_only),'p.value'])),linetype='dotted',color= input$dottedlinecol_volcano_only)
          }else{
            
            cutoff_volcano_only <- isolate(input$pvaluecutoff_volcano_only)
            adjdashline_volcano_only <- NULL
          }
          
          
          if(isolate(input$labelexpression_volcano_only)=='yes'){
            
            volcano_only[volcano_only$log2foldchange > log2(isolate(input$rfc_volcano_only)) & volcano_only$p.value < cutoff_volcano_only, 'color'] <- input$poscol_volcano_only
            volcano_only[volcano_only$log2foldchange < log2(isolate(input$lfc_volcano_only)) & volcano_only$p.value < cutoff_volcano_only, 'color'] <- input$negcol_volcano_only
            
          }else{
            
            volcano_only[volcano_only$log2foldchange > log2(isolate(input$rfc_volcano_only)) & volcano_only$p.value < cutoff_volcano_only, 'color'] <- input$sigcol_volcano_only
            volcano_only[volcano_only$log2foldchange < log2(isolate(input$lfc_volcano_only)) & volcano_only$p.value < cutoff_volcano_only, 'color'] <- input$sigcol_volcano_only
            
          }
          
          volcano_only[is.na(volcano_only$color),'color'] <- input$insigcol_volcano_only
          volcano_only[volcano_only$color== input$insigcol_volcano_only,'size']=0.2
          volcano_only[!(volcano_only$color== input$insigcol_volcano_only),'size']=1.5
          pal_volcano_only <- levels(factor(volcano_only$color))
          names(pal_volcano_only)<-pal_volcano_only
          volcano_only$log10pvalue <- -log10(volcano_only$p.value)
          
          if(isolate(input$set_x_boundary_volcano_only)=='yes'){
            
            xmin_val_volcano_only <- - isolate(input$x_axis_boundary_volcano_only)
            xmax_val_volcano_only <- isolate(input$x_axis_boundary_volcano_only)
            x_scale_volcano_only <- scale_x_continuous(name="log2(Fold Change)", limits = c(xmin_val_volcano_only,xmax_val_volcano_only))
          }else{
            x_scale_volcano_only <- scale_x_continuous(name="log2(Fold Change)")
          }
          
          yvec_volcano_only <- -log10(volcano_only$p.value)
          if(max(yvec_volcano_only)>20){yvec_volcano_only[yvec_volcano_only>20]<-20}
          ymax_val_volcano_only <- max(yvec_volcano_only)
          ymax_volcano_only =floor(max(yvec_volcano_only))+1
          
          dashline_volcano_only <- geom_hline(yintercept=-log10(isolate(input$pvaluecutoff_volcano_only)),linetype='dashed',color= input$dashedlinecol_volcano_only)
          
          
          p1 <- ggplot(data=volcano_only,aes(x=log2foldchange, y=log10pvalue, color=color, text=name)) +
            geom_point(show.legend=F,size=volcano_only$size,alpha=0.7) +
            dashline_volcano_only + adjdashline_volcano_only +
            geom_vline(xintercept = log2(isolate(input$lfc_volcano_only)), linetype=2, colour= input$lfcdashed_volcano_only) +
            geom_vline(xintercept = log2(isolate(input$rfc_volcano_only)), linetype=2, colour= input$rfcdashed_volcano_only) +
            x_scale_volcano_only +
            scale_y_continuous(name="-log10(p-value)", breaks = seq(0, (ymax_val_volcano_only + 2), by = isolate(input$y_axis_spacing_volcano_only)), limits = c(0,ymax_volcano_only)) +
            scale_color_manual(values=pal_volcano_only) +
            theme(axis.title=element_text(size=12),
                  axis.text =element_text(size=10),
                  axis.text.x = element_text(angle=0),
                  axis.line = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                  plot.title = element_text(size=14, hjust = 0.5),
                  legend.position = "none")
          
          closeAlert(session, "check_search_interactive_alert")
          p1
          
        }
        
        
      }
      
    }else{
      
      if(input$start1_interactive!=0  & check_interactive$count==0 & !is.null(feat_inte_interactive()) & !is.null(classlabel_inte_interactive()) & all_alert()){
        
        if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot'){
          
          if(isolate(input$yaxislabel_manhattan_box)=='pvalue'){
            
            if(isolate(input$adjdashline_manhattan_box)=='no'){
              manhattan_box = feat_inte_interactive()
              colnames(manhattan_box)[1:3] = c('name','log2foldchange','p.value')
            }else{
              manhattan_box = feat_inte_interactive()
              colnames(manhattan_box)[1:4] = c('name','log2foldchange','p.value','adjusted.p.value')
            }
            
            manhattan_box = manhattan_box[grep('_',manhattan_box$name),]
            manhattan_box = manhattan_box %>% separate('name',c("mz","time"),remove=FALSE,sep="_")
            manhattan_box$mz=as.numeric(manhattan_box$mz)
            manhattan_box$time=as.numeric(manhattan_box$time)
            
            if(isolate(input$adjdashline_manhattan_box)=='yes'){
              
              if(isolate(input$psignif_manhattan_box)=='pvalue'){
                
                cutoff_manhattan_box <- isolate(input$pvaluecutoff_manhattan_box)
              }else{
                
                cutoff_manhattan_box <- max(manhattan_box[manhattan_box$p.value<isolate(input$pvaluecutoff_manhattan_box) & manhattan_box$adjusted.p.value<isolate(input$adjpvaluecutoff_manhattan_box),'p.value'])
              }
              
              adjdashline_manhattan_box <- geom_hline(yintercept=-log10(max(manhattan_box[manhattan_box$p.value<isolate(input$pvaluecutoff_manhattan_box) & manhattan_box$adjusted.p.value<isolate(input$adjpvaluecutoff_manhattan_box),'p.value'])),linetype='dotted',color= input$dottedlinecol_manhattan_box)
            }else{
              
              cutoff_manhattan_box <- isolate(input$pvaluecutoff_manhattan_box)
              adjdashline_manhattan_box <- NULL
            }
            
            
            if(isolate(input$labelexpression_manhattan_box)=='yes'){
              
              manhattan_box[manhattan_box$p.value<cutoff_manhattan_box & manhattan_box$log2foldchange>0,'color']<- input$poscol_manhattan_box
              manhattan_box[manhattan_box$p.value<cutoff_manhattan_box & manhattan_box$log2foldchange<0,'color']<- input$negcol_manhattan_box
            }else{
              
              manhattan_box[manhattan_box$p.value<cutoff_manhattan_box,'color']<- input$sigcol_manhattan_box
            }
            
            manhattan_box[manhattan_box$p.value>=cutoff_manhattan_box,'color']<- input$insigcol_manhattan_box
            manhattan_box[manhattan_box$color== input$insigcol_manhattan_box,'size']=0.2
            manhattan_box[!(manhattan_box$color== input$insigcol_manhattan_box),'size']=1.5
            pal_manhattan_box <- levels(factor(manhattan_box$color))
            names(pal_manhattan_box)<-pal_manhattan_box
            
            if(isolate(input$plottype_manhattan_box)=='type1'){
              
              xincrement_manhattan_box = isolate(input$x_axis_spacing_type1_manhattan_box)
              xmin_val_manhattan_box <- min(0,manhattan_box$mz)
              xmax_val_manhattan_box <- max(manhattan_box$mz)
              xlabel_manhattan_box <- xlab('mass-to-charge (m/z)')
              maintitle_manhattan_box <- ggtitle('Type 1 manhattan plot (-log10p vs mz)')
              gplotline_manhattan_box <- ggplot(manhattan_box, aes(x=mz, y=-log10(p.value), color=color, text=paste0('mz:',manhattan_box$mz," time:",manhattan_box$time)))
              
            }else{
              
              xincrement_manhattan_box = isolate(input$x_axis_spacing_type2_manhattan_box)
              xmin_val_manhattan_box <- min(0,manhattan_box$time)
              xmax_val_manhattan_box <- max(manhattan_box$time)
              xlabel_manhattan_box <- xlab('Retention time (s)')
              maintitle_manhattan_box <- ggtitle('Type 2 manhattan plot (-log10p vs time)')
              gplotline_manhattan_box <- ggplot(manhattan_box, aes(x=time, y=-log10(p.value), color=color, text=paste0('mz:',manhattan_box$mz," time:",manhattan_box$time)))
              
            }
            

            yvec_manhattan_box <- -log10(manhattan_box$p.value)
            if(max(yvec_manhattan_box)>20){yvec_manhattan_box[yvec_manhattan_box>20]<-20}
            ymax_val_manhattan_box <- max(yvec_manhattan_box)
            ymax_manhattan_box =floor(max(yvec_manhattan_box))+1
            ylabel_manhattan_box <- ylab('-log10p')
            
            dashline_manhattan_box <- geom_hline(yintercept=-log10(isolate(input$pvaluecutoff_manhattan_box)),linetype='dashed',color= input$dashedlinecol_manhattan_box)
            
          }else{
            
            manhattan_box = feat_inte_interactive()
            colnames(manhattan_box)[1:3] = c('name','log2foldchange','VIP')
            manhattan_box=manhattan_box[grep("_",manhattan_box$name),]
            rownames(manhattan_box)=1:nrow(manhattan_box)
            manhattan_box = manhattan_box %>% separate(name,c("mz","time"),remove=FALSE,sep="_")
            manhattan_box$mz=as.numeric(manhattan_box$mz)
            manhattan_box$time=as.numeric(manhattan_box$time)
            
            if(isolate(input$labelexpression_manhattan_box)=='yes'){
              
              manhattan_box[manhattan_box$VIP > isolate(input$vipcutoff_manhattan_box) & manhattan_box$log2foldchange>0,'color']<- input$poscol_manhattan_box
              manhattan_box[manhattan_box$VIP > isolate(input$vipcutoff_manhattan_box) & manhattan_box$log2foldchange<0,'color']<- input$negcol_manhattan_box
            }else{
              
              manhattan_box[manhattan_box$VIP > isolate(input$vipcutoff_manhattan_box),'color']<- input$sigcol_manhattan_box
            }
            
            manhattan_box[manhattan_box$VIP<= isolate(input$vipcutoff_manhattan_box),'color']<- input$insigcol_manhattan_box
            manhattan_box[manhattan_box$color== input$insigcol_manhattan_box,'size']=0.2
            manhattan_box[!(manhattan_box$color== input$insigcol_manhattan_box),'size']=1.5
            pal_manhattan_box <- levels(factor(manhattan_box$color))
            names(pal_manhattan_box)<-pal_manhattan_box
            
            if(isolate(input$plottype_manhattan_box)=='type1'){
              
              xincrement_manhattan_box = isolate(input$x_axis_spacing_type1_manhattan_box)
              xmin_val_manhattan_box <- min(0,manhattan_box$mz)
              xmax_val_manhattan_box <- max(manhattan_box$mz)
              xlabel_manhattan_box <- xlab('mass-to-charge (m/z)')
              maintitle_manhattan_box <- ggtitle('Type 1 manhattan plot (VIP vs mz)')
              gplotline_manhattan_box <- ggplot(manhattan_box, aes(x=mz, y=VIP, color=color, text=paste0('mz:',manhattan_box$mz," time:",manhattan_box$time)))
              
            }else{
              
              xincrement_manhattan_box = isolate(input$x_axis_spacing_type2_manhattan_box)
              xmin_val_manhattan_box <- min(0,manhattan_box$time)
              xmax_val_manhattan_box <- max(manhattan_box$time)
              xlabel_manhattan_box <- xlab('Retention time (s)')
              maintitle_manhattan_box <- ggtitle('Type 2 manhattan plot (VIP vs time)')
              gplotline_manhattan_box <- ggplot(manhattan_box, aes(x=time, y=VIP, color=color, text=paste0('mz:',manhattan_box$mz," time:",manhattan_box$time)))
              
            }
            
            yvec_manhattan_box <- manhattan_box$VIP
            if(max(yvec_manhattan_box)>20){yvec_manhattan_box[yvec_manhattan_box>20]<-20}
            ymax_val_manhattan_box <- max(yvec_manhattan_box)
            ymax_manhattan_box =floor(max(yvec_manhattan_box))+1
            ylabel_manhattan_box <- ylab('VIP')
            
            dashline_manhattan_box <- geom_hline(yintercept= isolate(input$vipcutoff_manhattan_box),linetype='dashed',color= input$dashedlinecol_manhattan_box)
            
            adjdashline_manhattan_box <- NULL
            
          }
          
          
          p1 <- gplotline_manhattan_box +
            geom_point(show.legend=F,size=manhattan_box$size,alpha=0.7) +
            xlabel_manhattan_box + ylabel_manhattan_box +
            scale_color_manual(values=pal_manhattan_box) +
            scale_x_continuous(breaks = seq(xmin_val_manhattan_box, xmax_val_manhattan_box, by = xincrement_manhattan_box)) +
            scale_y_continuous(breaks = seq(0, (ymax_val_manhattan_box + 2), by = isolate(input$y_axis_spacing_manhattan_box)), limits = c(0,ymax_manhattan_box)) +
            dashline_manhattan_box + adjdashline_manhattan_box +
            maintitle_manhattan_box +
            theme(axis.title=element_text(size=12),
                  axis.text =element_text(size=10),
                  axis.text.x = element_text(angle=0),
                  axis.line = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                  plot.title = element_text(size=14, hjust = 0.5),
                  legend.position = "none")
          
          closeAlert(session, "check_search_interactive_alert")
          p1
          
        }else{
          
          
          if(ncol(feat_inte_interactive())==3){
            volcano_box = feat_inte_interactive()
            colnames(volcano_box)[1:3] = c('name','log2foldchange','p.value')
          }else{
            volcano_box = feat_inte_interactive()
            colnames(volcano_box)[1:4] = c('name','log2foldchange','p.value','adjusted.p.value')
          }
          
          volcano_box = volcano_box[grep('_',volcano_box$name),]
          volcano_box = volcano_box %>% separate('name',c("mz","time"),remove=FALSE,sep="_")
          volcano_box$mz=as.numeric(volcano_box$mz)
          volcano_box$time=as.numeric(volcano_box$time)
          
          if(isolate(input$adjdashline_volcano_box)=='yes'){
            
            if(isolate(input$psignif_volcano_box)=='pvalue'){
              
              cutoff_volcano_box <- isolate(input$pvaluecutoff_volcano_box)
            }else{
              
              cutoff_volcano_box <- max(volcano_box[volcano_box$p.value<isolate(input$pvaluecutoff_volcano_box) & volcano_box$adjusted.p.value<isolate(input$adjpvaluecutoff_volcano_box),'p.value'])
            }
            
            adjdashline_volcano_box <- geom_hline(yintercept=-log10(max(volcano_box[volcano_box$p.value<isolate(input$pvaluecutoff_volcano_box) & volcano_box$adjusted.p.value<isolate(input$adjpvaluecutoff_volcano_box),'p.value'])),linetype='dotted',color= input$dottedlinecol_volcano_box)
          }else{
            
            cutoff_volcano_box <- isolate(input$pvaluecutoff_volcano_box)
            adjdashline_volcano_box <- NULL
          }
          
          
          if(isolate(input$labelexpression_volcano_box)=='yes'){
            
            volcano_box[volcano_box$log2foldchange > log2(isolate(input$rfc_volcano_box)) & volcano_box$p.value < cutoff_volcano_box, 'color'] <- input$poscol_volcano_box
            volcano_box[volcano_box$log2foldchange < log2(isolate(input$lfc_volcano_box)) & volcano_box$p.value < cutoff_volcano_box, 'color'] <- input$negcol_volcano_box
            
          }else{
            
            volcano_box[volcano_box$log2foldchange > log2(isolate(input$rfc_volcano_box)) & volcano_box$p.value < cutoff_volcano_box, 'color'] <- input$sigcol_volcano_box
            volcano_box[volcano_box$log2foldchange < log2(isolate(input$lfc_volcano_box)) & volcano_box$p.value < cutoff_volcano_box, 'color'] <- input$sigcol_volcano_box
            
          }
          
          volcano_box[is.na(volcano_box$color),'color'] <- input$insigcol_volcano_box
          volcano_box[volcano_box$color== input$insigcol_volcano_box,'size']=0.2
          volcano_box[!(volcano_box$color== input$insigcol_volcano_box),'size']=1.5
          pal_volcano_box <- levels(factor(volcano_box$color))
          names(pal_volcano_box)<-pal_volcano_box
          volcano_box$log10pvalue <- -log10(volcano_box$p.value)
          
          if(isolate(input$set_x_boundary_volcano_box)=='yes'){
            
            xmin_val_volcano_box <- - isolate(input$x_axis_boundary_volcano_box)
            xmax_val_volcano_box <- isolate(input$x_axis_boundary_volcano_box)
            x_scale_volcano_box <- scale_x_continuous(name="log2(Fold Change)", limits = c(xmin_val_volcano_box,xmax_val_volcano_box))
          }else{
            x_scale_volcano_box <- scale_x_continuous(name="log2(Fold Change)")
          }
          
          yvec_volcano_box <- -log10(volcano_box$p.value)
          if(max(yvec_volcano_box)>20){yvec_volcano_box[yvec_volcano_box>20]<-20}
          ymax_val_volcano_box <- max(yvec_volcano_box)
          ymax_volcano_box =floor(max(yvec_volcano_box))+1
          
          dashline_volcano_box <- geom_hline(yintercept=-log10(isolate(input$pvaluecutoff_volcano_box)),linetype='dashed',color= input$dashedlinecol_volcano_box)
          
          
          p1 <- ggplot(data=volcano_box,aes(x=log2foldchange, y=log10pvalue, color=color, text=paste0('mz:',volcano_box$mz," time:",volcano_box$time))) +
            geom_point(show.legend=F,size=volcano_box$size,alpha=0.7) +
            dashline_volcano_box + adjdashline_volcano_box +
            geom_vline(xintercept = log2(isolate(input$lfc_volcano_box)), linetype=2, colour= input$lfcdashed_volcano_box) +
            geom_vline(xintercept = log2(isolate(input$rfc_volcano_box)), linetype=2, colour= input$rfcdashed_volcano_box) +
            x_scale_volcano_box +
            scale_y_continuous(name="-log10(p-value)", breaks = seq(0, (ymax_val_volcano_box + 2), by = isolate(input$y_axis_spacing_volcano_box)), limits = c(0,ymax_volcano_box)) +
            scale_color_manual(values=pal_volcano_box) +
            theme(axis.title=element_text(size=12),
                  axis.text =element_text(size=10),
                  axis.text.x = element_text(angle=0),
                  axis.line = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                  plot.title = element_text(size=14, hjust = 0.5),
                  legend.position = "none")
          
          closeAlert(session, "check_search_interactive_alert")
          p1
          
        }
        
      }
      
    }
    
  })
  
  
  plot2_single_plot <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot only' || isolate(input$choose_graph)=='Volcano Plot only'){
      
      if(input$start1_interactive!=0  & check_interactive$count==0 & !is.null(input_file_interactive()) & all_alert()){
        
        if(isolate(input$choose_graph)=='Manhattan Plot only'){
          
          if(isolate(input$yaxislabel_manhattan_only)=='pvalue'){
            
            if(ncol(input_file_interactive())==3){
              manhattan_only = input_file_interactive()
              colnames(manhattan_only) = c('name','log2foldchange','p.value')
            }else{
              manhattan_only = input_file_interactive()
              colnames(manhattan_only) = c('name','log2foldchange','p.value','adjusted.p.value')
            }
            
            manhattan_only = manhattan_only[grep('_',manhattan_only$name),]
            manhattan_only = manhattan_only %>% separate('name',c("mz","time"),remove=FALSE,sep="_")
            manhattan_only$mz=as.numeric(manhattan_only$mz)
            manhattan_only$time=as.numeric(manhattan_only$time)
            
            if(isolate(input$adjdashline_manhattan_only)=='yes'){
              
              if(isolate(input$psignif_manhattan_only)=='pvalue'){
                
                cutoff_manhattan_only <- isolate(input$pvaluecutoff_manhattan_only)
              }else{
                
                cutoff_manhattan_only <- max(manhattan_only[manhattan_only$p.value<isolate(input$pvaluecutoff_manhattan_only) & manhattan_only$adjusted.p.value<isolate(input$adjpvaluecutoff_manhattan_only),'p.value'])
              }
              
              adjdashline_manhattan_only <- geom_hline(yintercept=-log10(max(manhattan_only[manhattan_only$p.value<isolate(input$pvaluecutoff_manhattan_only) & manhattan_only$adjusted.p.value<isolate(input$adjpvaluecutoff_manhattan_only),'p.value'])),linetype='dotted',color= input$dottedlinecol_manhattan_only)
            }else{
              
              cutoff_manhattan_only <- isolate(input$pvaluecutoff_manhattan_only)
              adjdashline_manhattan_only <- NULL
            }
            
            
            if(isolate(input$labelexpression_manhattan_only)=='yes'){
              
              manhattan_only[manhattan_only$p.value<cutoff_manhattan_only & manhattan_only$log2foldchange>0,'color']<- input$poscol_manhattan_only
              manhattan_only[manhattan_only$p.value<cutoff_manhattan_only & manhattan_only$log2foldchange<0,'color']<- input$negcol_manhattan_only
            }else{
              
              manhattan_only[manhattan_only$p.value<cutoff_manhattan_only,'color']<- input$sigcol_manhattan_only
            }
            
            manhattan_only[manhattan_only$p.value>=cutoff_manhattan_only,'color']<- input$insigcol_manhattan_only
            manhattan_only[manhattan_only$color== input$insigcol_manhattan_only,'size']=0.2
            manhattan_only[!(manhattan_only$color== input$insigcol_manhattan_only),'size']=1.5
            pal_manhattan_only <- levels(factor(manhattan_only$color))
            names(pal_manhattan_only)<-pal_manhattan_only
            
            xincrement_manhattan_only = isolate(input$x_axis_spacing_type2_manhattan_only) #round_any(max(manhattan$time)/10,10,f=floor)
            xmin_val_manhattan_only <- min(0,manhattan_only$time)
            xmax_val_manhattan_only <- max(manhattan_only$time)
            xlabel_manhattan_only <- xlab('Retention time (s)')
            maintitle_manhattan_only <- ggtitle('Type 2 manhattan plot (-log10p vs time)')
            gplotline_manhattan_only <- ggplot(manhattan_only, aes(x=time, y=-log10(p.value), color=color, text=paste0('mz:',manhattan_only$mz," time:",manhattan_only$time)))
            
            yvec_manhattan_only <- -log10(manhattan_only$p.value)
            if(max(yvec_manhattan_only)>20){yvec_manhattan_only[yvec_manhattan_only>20]<-20}
            ymax_val_manhattan_only <- max(yvec_manhattan_only)
            ymax_manhattan_only =floor(max(yvec_manhattan_only))+1
            ylabel_manhattan_only <- ylab('-log10p')
            
            dashline_manhattan_only <- geom_hline(yintercept=-log10(isolate(input$pvaluecutoff_manhattan_only)),linetype='dashed',color= input$dashedlinecol_manhattan_only)
            
          }else{
            
            manhattan_only = input_file_interactive()
            colnames(manhattan_only) = c('name','log2foldchange','VIP')
            manhattan_only=manhattan_only[grep("_",manhattan_only$name),]
            rownames(manhattan_only)=1:nrow(manhattan_only)
            manhattan_only = manhattan_only %>% separate(name,c("mz","time"),remove=FALSE,sep="_")
            manhattan_only$mz=as.numeric(manhattan_only$mz)
            manhattan_only$time=as.numeric(manhattan_only$time)
            
            if(isolate(input$labelexpression_manhattan_only)=='yes'){
              
              manhattan_only[manhattan_only$VIP > isolate(input$vipcutoff_manhattan_only) & manhattan_only$log2foldchange>0,'color']<- input$poscol_manhattan_only
              manhattan_only[manhattan_only$VIP > isolate(input$vipcutoff_manhattan_only) & manhattan_only$log2foldchange<0,'color']<- input$negcol_manhattan_only
            }else{
              
              manhattan_only[manhattan_only$VIP > isolate(input$vipcutoff_manhattan_only),'color']<- input$sigcol_manhattan_only
            }
            
            manhattan_only[manhattan_only$VIP<= isolate(input$vipcutoff_manhattan_only),'color']<- input$insigcol_manhattan_only
            manhattan_only[manhattan_only$color== input$insigcol_manhattan_only,'size']=0.2
            manhattan_only[!(manhattan_only$color== input$insigcol_manhattan_only),'size']=1.5
            pal_manhattan_only <- levels(factor(manhattan_only$color))
            names(pal_manhattan_only)<-pal_manhattan_only
            
            xincrement_manhattan_only = isolate(input$x_axis_spacing_type2_manhattan_only)
            xmin_val_manhattan_only <- min(0,manhattan_only$time)
            xmax_val_manhattan_only <- max(manhattan_only$time)
            xlabel_manhattan_only <- xlab('Retention time (s)')
            maintitle_manhattan_only <- ggtitle('Type 2 manhattan plot (VIP vs time)')
            gplotline_manhattan_only <- ggplot(manhattan_only, aes(x=time, y=VIP, color=color, text=paste0('mz:',manhattan_only$mz," time:",manhattan_only$time)))
            
            yvec_manhattan_only <- manhattan_only$VIP
            if(max(yvec_manhattan_only)>20){yvec_manhattan_only[yvec_manhattan_only>20]<-20}
            ymax_val_manhattan_only <- max(yvec_manhattan_only)
            ymax_manhattan_only =floor(max(yvec_manhattan_only))+1
            ylabel_manhattan_only <- ylab('VIP')
            
            dashline_manhattan_only <- geom_hline(yintercept= isolate(input$vipcutoff_manhattan_only),linetype='dashed',color= input$dashedlinecol_manhattan_only)
            
            adjdashline_manhattan_only <- NULL
            
          }
          
          
          p2 <- gplotline_manhattan_only +
            geom_point(show.legend=F,size=manhattan_only$size,alpha=0.7) +
            xlabel_manhattan_only + ylabel_manhattan_only +
            scale_color_manual(values=pal_manhattan_only) +
            scale_x_continuous(breaks = seq(xmin_val_manhattan_only, xmax_val_manhattan_only, by = xincrement_manhattan_only)) +
            scale_y_continuous(breaks = seq(0, (ymax_val_manhattan_only + 2), by = isolate(input$y_axis_spacing_manhattan_only)), limits = c(0,ymax_manhattan_only)) +
            dashline_manhattan_only + adjdashline_manhattan_only +
            maintitle_manhattan_only +
            theme(axis.title=element_text(size=12),
                  axis.text =element_text(size=10),
                  axis.text.x = element_text(angle=0),
                  axis.line = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                  plot.title = element_text(size=14, hjust = 0.5),
                  legend.position = "none")
          
          closeAlert(session, "check_search_interactive_alert")
          p2
          
        }else{
          
          closeAlert(session, "check_search_interactive_alert")
          NULL
        }
      
      }
      
    }else{
        
        if(input$start1_interactive!=0  & check_interactive$count==0 & !is.null(feat_inte_interactive()) & !is.null(classlabel_inte_interactive()) & all_alert()){
          
          if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot'){
            
            if(!is.null(pplot1())){
              
              d <- event_data('plotly_click',source = 'scatterplot')
              
              if(!is.null(d)){
                
                if(isolate(input$yaxislabel_manhattan_box)=='pvalue'){
                  
                  if(isolate(input$adjdashline_manhattan_box)=='no'){
                    manhattan_box = feat_inte_interactive()
                    colnames(manhattan_box)[1:3] = c('name','log2foldchange','p.value')
                  }else{
                    manhattan_box = feat_inte_interactive()
                    colnames(manhattan_box)[1:4] = c('name','log2foldchange','p.value','adjusted.p.value')
                  }
                  
                }else{
                  
                  manhattan_box = feat_inte_interactive()
                  colnames(manhattan_box)[1:3] = c('name','log2foldchange','VIP')
                  
                }
                
                manhattan_box=manhattan_box[grep("_",manhattan_box$name),]
                rownames(manhattan_box)=1:nrow(manhattan_box)
                
                x1_manhattan_box <-as.numeric(d['curveNumber']+1);x2_manhattan_box <-as.numeric(d['pointNumber']+1)
                text_manhattan_box <- pplot1()$x$data[[x1_manhattan_box]]$text[x2_manhattan_box]
                mz_manhattan_box <- unlist(strsplit(unlist(strsplit(text_manhattan_box,split = ' '))[1],split = ':'))[2]
                time_manhattan_box <- unlist(strsplit(unlist(strsplit(text_manhattan_box,split = ' '))[2],split = ':'))[2]
                
                if(!is.na(mz_manhattan_box) & !is.na(time_manhattan_box)){
                  
                  name_manhattan_box  = paste0(mz_manhattan_box,"_",time_manhattan_box)
                  group_manhattan_box <- names(table(classlabel_inte_interactive()[,2]))
                  boxplotdata_manhattan_box <- data.frame(matrix(NA,ncol=2,nrow=1))
                  colnames(boxplotdata_manhattan_box) <- c('group','intensity')
                  
                  for(i in 1:length(group_manhattan_box)){
                    
                    tmp2_manhattan_box <- as.numeric(manhattan_box[manhattan_box$name==name_manhattan_box,classlabel_inte_interactive()[classlabel_inte_interactive()[,2]==group_manhattan_box[i],1]])
                    tmp1_manhattan_box <- rep(group_manhattan_box[i],length(tmp2_manhattan_box))
                    tmpdata_manhattan_box <- cbind(tmp1_manhattan_box,tmp2_manhattan_box)
                    colnames(tmpdata_manhattan_box) <- c('group','intensity')
                    boxplotdata_manhattan_box <- rbind(boxplotdata_manhattan_box,tmpdata_manhattan_box)
                  }
                  boxplotdata_manhattan_box <- boxplotdata_manhattan_box[-1,]
                  boxplotdata_manhattan_box$intensity <- as.numeric(boxplotdata_manhattan_box$intensity)
                  
                  if(isolate(input$boxplotcolor_manhattan_box)=='yes'){
                    
                    if(isolate(input$xaxis_name_manhattan_box_1)==""){
                      xlab_manhattan_box <- xlab('group')
                    }else{
                      xlab_manhattan_box <- xlab(isolate(input$xaxis_name_manhattan_box_1))
                    }
                    if(isolate(input$yaxis_name_manhattan_box_1)==""){
                      ylab_manhattan_box <- ylab('intensity')
                    }else{
                      ylab_manhattan_box <-  ylab(isolate(input$yaxis_name_manhattan_box_1))
                    }
                    
                  }else if(isolate(input$boxplotcolor_manhattan_box)=='no'){
                    
                    if(isolate(input$xaxis_name_manhattan_box_2)==""){
                      xlab_manhattan_box <- xlab('group')
                    }else{
                      xlab_manhattan_box <- xlab(isolate(input$xaxis_name_manhattan_box_2))
                    }
                    if(isolate(input$yaxis_name_manhattan_box_2)==""){
                      ylab_manhattan_box <- ylab('intensity')
                    }else{
                      ylab_manhattan_box <-  ylab(isolate(input$yaxis_name_manhattan_box_2))
                    }
                    
                  }
                  
                  if(isolate(input$boxplotcolor_manhattan_box)=='yes'){
                    p2 <- ggboxplot(boxplotdata_manhattan_box, x = "group", y = "intensity", color = "group", bxp.errorbar = TRUE, bxp.errorbar.width=0.5, size = 0.5) +
                      ggtitle(text_manhattan_box) + xlab_manhattan_box + ylab_manhattan_box +
                      theme(axis.title=element_text(size=12),
                            axis.text =element_text(size=10),
                            axis.text.x = element_text(angle=0),
                            axis.line = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                            plot.title = element_text(size=14, hjust = 0.5))
                  }else{
                    p2 <- ggboxplot(boxplotdata_manhattan_box, x = "group", y = "intensity", bxp.errorbar = TRUE, bxp.errorbar.width=0.5, size = 0.5) +
                      ggtitle(text_manhattan_box) + xlab_manhattan_box + ylab_manhattan_box +
                      theme(axis.title=element_text(size=12),
                            axis.text =element_text(size=10),
                            axis.text.x = element_text(angle=0),
                            axis.line = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                            plot.title = element_text(size=14, hjust = 0.5))
                  }
                  
                  closeAlert(session, "check_search_interactive_alert")
                  p2
                  
                }else{
                  
                  closeAlert(session, "check_search_interactive_alert")
                  NULL
                }
                
              }
            }
            
            
          }else{
            
            if(!is.null(pplot1())){
              
              d <- event_data('plotly_click',source = 'scatterplot')
              
              if(!is.null(d)){
                  
                  if(isolate(input$adjdashline_volcano_box)=='no'){
                    volcano_box = feat_inte_interactive()
                    colnames(volcano_box)[1:3] = c('name','log2foldchange','p.value')
                  }else{
                    volcano_box = feat_inte_interactive()
                    colnames(volcano_box)[1:4] = c('name','log2foldchange','p.value','adjusted.p.value')
                  }
                
                volcano_box=volcano_box[grep("_",volcano_box$name),]
                rownames(volcano_box)=1:nrow(volcano_box)
                
                x1_volcano_box <-as.numeric(d['curveNumber']+1);x2_volcano_box <-as.numeric(d['pointNumber']+1)
                text_volcano_box <- pplot1()$x$data[[x1_volcano_box]]$text[x2_volcano_box]
                mz_volcano_box <- unlist(strsplit(unlist(strsplit(text_volcano_box,split = ' '))[1],split = ':'))[2]
                time_volcano_box <- unlist(strsplit(unlist(strsplit(text_volcano_box,split = ' '))[2],split = ':'))[2]
                
                if(!is.na(mz_volcano_box) & !is.na(time_volcano_box)){
                  
                  name_volcano_box  = paste0(mz_volcano_box,"_",time_volcano_box)
                  group_volcano_box <- names(table(classlabel_inte_interactive()[,2]))
                  boxplotdata_volcano_box <- data.frame(matrix(NA,ncol=2,nrow=1))
                  colnames(boxplotdata_volcano_box) <- c('group','intensity')
                  
                  for(i in 1:length(group_volcano_box)){
                    
                    tmp2_volcano_box <- as.numeric(volcano_box[volcano_box$name==name_volcano_box,classlabel_inte_interactive()[classlabel_inte_interactive()[,2]==group_volcano_box[i],1]])
                    tmp1_volcano_box <- rep(group_volcano_box[i],length(tmp2_volcano_box))
                    tmpdata_volcano_box <- cbind(tmp1_volcano_box,tmp2_volcano_box)
                    colnames(tmpdata_volcano_box) <- c('group','intensity')
                    boxplotdata_volcano_box <- rbind(boxplotdata_volcano_box,tmpdata_volcano_box)
                  }
                  boxplotdata_volcano_box <- boxplotdata_volcano_box[-1,]
                  boxplotdata_volcano_box$intensity <- as.numeric(boxplotdata_volcano_box$intensity)
                  
                  if(isolate(input$boxplotcolor_volcano_box)=='yes'){
                    
                    if(isolate(input$xaxis_name_volcano_box_1)==""){
                      xlab_volcano_box <- xlab('group')
                    }else{
                      xlab_volcano_box <- xlab(isolate(input$xaxis_name_volcano_box_1))
                    }
                    if(isolate(input$yaxis_name_volcano_box_1)==""){
                      ylab_volcano_box <- ylab('intensity')
                    }else{
                      ylab_volcano_box <-  ylab(isolate(input$yaxis_name_volcano_box_1))
                    }
                    
                  }else if(isolate(input$boxplotcolor_volcano_box)=='no'){
                    
                    if(isolate(input$xaxis_name_volcano_box_2)==""){
                      xlab_volcano_box <- xlab('group')
                    }else{
                      xlab_volcano_box <- xlab(isolate(input$xaxis_name_volcano_box_2))
                    }
                    if(isolate(input$yaxis_name_volcano_box_2)==""){
                      ylab_volcano_box <- ylab('intensity')
                    }else{
                      ylab_volcano_box <-  ylab(isolate(input$yaxis_name_volcano_box_2))
                    }
                    
                  }
                  
                  if(isolate(input$boxplotcolor_volcano_box)=='yes'){
                    p2 <- ggboxplot(boxplotdata_volcano_box, x = "group", y = "intensity", color = "group", bxp.errorbar = TRUE, bxp.errorbar.width=0.5, size = 0.5) +
                      ggtitle(text_volcano_box) + xlab_volcano_box + ylab_volcano_box +
                      theme(axis.title=element_text(size=12),
                            axis.text =element_text(size=10),
                            axis.text.x = element_text(angle=0),
                            axis.line = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                            plot.title = element_text(size=14, hjust = 0.5))
                  }else{
                    p2 <- ggboxplot(boxplotdata_volcano_box, x = "group", y = "intensity", bxp.errorbar = TRUE, bxp.errorbar.width=0.5, size = 0.5) +
                      ggtitle(text_volcano_box) + xlab_volcano_box + ylab_volcano_box +
                      theme(axis.title=element_text(size=12),
                            axis.text =element_text(size=10),
                            axis.text.x = element_text(angle=0),
                            axis.line = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                            plot.title = element_text(size=14, hjust = 0.5))
                  }
                  
                  closeAlert(session, "check_search_interactive_alert")
                  p2
                  
                }else{
                  
                  closeAlert(session, "check_search_interactive_alert")
                  NULL
                }
                
              }
            }
            
          }
          
        }
    }
    
  })
  
  
  observeEvent(input$search_manhattan_box,{search_interactive$count=1})
  observeEvent(input$search_volcano_box,{search_interactive$count=1})
  
  
  plot3_single_plot <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot' || isolate(input$choose_graph)=='Volcano Plot with Box Plot'){
    
    if(input$start1_interactive!=0  & check_interactive$count==0 & !is.null(feat_inte_interactive()) & !is.null(classlabel_inte_interactive()) & all_alert()){
      
      if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot'){
        
        if(!is.null(pplot1())){
          
          if(input$search_manhattan_box!=0 & search_interactive$count==1){
            
            if(!str_trim(tolower(isolate(input$name_manhattan_box)))==""){
            
              if(isolate(input$yaxislabel_manhattan_box)=='pvalue'){
                
                if(isolate(input$adjdashline_manhattan_box)=='no'){
                  manhattan_box = feat_inte_interactive()
                  colnames(manhattan_box)[1:3] = c('name','log2foldchange','p.value')
                }else{
                  manhattan_box = feat_inte_interactive()
                  colnames(manhattan_box)[1:4] = c('name','log2foldchange','p.value','adjusted.p.value')
                }
                
              }else{
                
                manhattan_box = feat_inte_interactive()
                colnames(manhattan_box)[1:3] = c('name','log2foldchange','VIP')
                
              }
              
              if(sum(tolower(manhattan_box$name)==str_trim(tolower(isolate(input$name_manhattan_box))))>0){
                
                group_manhattan_box <- names(table(classlabel_inte_interactive()[,2]))
                boxplotdata_manhattan_box <- data.frame(matrix(NA,ncol=2,nrow=1))
                colnames(boxplotdata_manhattan_box) <- c('group','intensity')
                
                for(i in 1:length(group_manhattan_box)){
                  
                  tmp2_manhattan_box <- as.numeric(manhattan_box[tolower(manhattan_box$name)==str_trim(tolower(isolate(input$name_manhattan_box))),classlabel_inte_interactive()[classlabel_inte_interactive()[,2]==group_manhattan_box[i],1]])
                  tmp1_manhattan_box <- rep(group_manhattan_box[i],length(tmp2_manhattan_box))
                  tmpdata_manhattan_box <- cbind(tmp1_manhattan_box,tmp2_manhattan_box)
                  colnames(tmpdata_manhattan_box) <- c('group','intensity')
                  boxplotdata_manhattan_box <- rbind(boxplotdata_manhattan_box,tmpdata_manhattan_box)
                }
                boxplotdata_manhattan_box <- boxplotdata_manhattan_box[-1,]
                boxplotdata_manhattan_box$intensity <- as.numeric(boxplotdata_manhattan_box$intensity)
                
                if(length(grep("_",str_trim(tolower(isolate(input$name_manhattan_box)))))>0){
                  
                  text_manhattan_box <- paste('mz:',strsplit(str_trim(tolower(isolate(input$name_manhattan_box))),"_")[[1]][1],' time:',strsplit(str_trim(tolower(isolate(input$name_manhattan_box))),"_")[[1]][2],sep='')
                }else{
                  
                  text_manhattan_box <- as.character(manhattan_box[tolower(manhattan_box$name)==str_trim(tolower(isolate(input$name_manhattan_box))),'name'])
                }
                
                if(isolate(input$boxplotcolor_manhattan_box)=='yes'){
                  
                  if(isolate(input$xaxis_name_manhattan_box_1)==""){
                    xlab_manhattan_box <- xlab('group')
                  }else{
                    xlab_manhattan_box <- xlab(isolate(input$xaxis_name_manhattan_box_1))
                  }
                  if(isolate(input$yaxis_name_manhattan_box_1)==""){
                    ylab_manhattan_box <- ylab('intensity')
                  }else{
                    ylab_manhattan_box <-  ylab(isolate(input$yaxis_name_manhattan_box_1))
                  }
                  
                }else if(isolate(input$boxplotcolor_manhattan_box)=='no'){
                  
                  if(isolate(input$xaxis_name_manhattan_box_2)==""){
                    xlab_manhattan_box <- xlab('group')
                  }else{
                    xlab_manhattan_box <- xlab(isolate(input$xaxis_name_manhattan_box_2))
                  }
                  if(isolate(input$yaxis_name_manhattan_box_2)==""){
                    ylab_manhattan_box <- ylab('intensity')
                  }else{
                    ylab_manhattan_box <-  ylab(isolate(input$yaxis_name_manhattan_box_2))
                  }
                  
                }
                
                
                if(isolate(input$boxplotcolor_manhattan_box)=='yes'){
                  p3 <- ggboxplot(boxplotdata_manhattan_box, x = "group", y = "intensity", color = "group", bxp.errorbar = TRUE, bxp.errorbar.width=0.5, size = 0.5) +
                    ggtitle(text_manhattan_box) + xlab_manhattan_box + ylab_manhattan_box +
                    theme(axis.title=element_text(size=12),
                          axis.text =element_text(size=10),
                          axis.text.x = element_text(angle=0),
                          axis.line = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                          plot.title = element_text(size=14, hjust = 0.5))
                }else{
                  p3 <- ggboxplot(boxplotdata_manhattan_box, x = "group", y = "intensity", bxp.errorbar = TRUE, bxp.errorbar.width=0.5, size = 0.5) +
                    ggtitle(text_manhattan_box) + xlab_manhattan_box + ylab_manhattan_box +
                    theme(axis.title=element_text(size=12),
                          axis.text =element_text(size=10),
                          axis.text.x = element_text(angle=0),
                          axis.line = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                          plot.title = element_text(size=14, hjust = 0.5))
                }
                
                closeAlert(session, "check_search_interactive_alert")
                p3
                
              }else{
                
                closeAlert(session, "check_search_interactive_alert")
                createAlert(session, "alert_interactive", "check_search_interactive_alert", title = "Search Error", content = paste0("No metabolite with ",str_trim(isolate(input$name_manhattan_box))," was found, please check your input file."), append = TRUE)
                NULL
              }
              
            }else{
              
              closeAlert(session, "check_search_interactive_alert")
              createAlert(session, "alert_interactive", "check_search_interactive_alert", title = "Search Error", content = paste0("No search text was entered."), append = TRUE)
              NULL
            }
          }
        }
        
      }else{
        
        if(!is.null(pplot1())){
          
          if(input$search_volcano_box!=0 & search_interactive$count==1){
            
            if(!str_trim(tolower(isolate(input$name_volcano_box)))==""){
                
                if(isolate(input$adjdashline_volcano_box)=='no'){
                  volcano_box = feat_inte_interactive()
                  colnames(volcano_box)[1:3] = c('name','log2foldchange','p.value')
                }else{
                  volcano_box = feat_inte_interactive()
                  colnames(volcano_box)[1:4] = c('name','log2foldchange','p.value','adjusted.p.value')
                }
              
              
                if(sum(tolower(volcano_box$name)==str_trim(tolower(isolate(input$name_volcano_box))))>0){
                  
                  group_volcano_box <- names(table(classlabel_inte_interactive()[,2]))
                  boxplotdata_volcano_box <- data.frame(matrix(NA,ncol=2,nrow=1))
                  colnames(boxplotdata_volcano_box) <- c('group','intensity')
                  
                  
                  for(i in 1:length(group_volcano_box)){
                    
                    tmp2_volcano_box <- as.numeric(volcano_box[tolower(volcano_box$name)==str_trim(tolower(isolate(input$name_volcano_box))),classlabel_inte_interactive()[classlabel_inte_interactive()[,2]==group_volcano_box[i],1]])
                    tmp1_volcano_box <- rep(group_volcano_box[i],length(tmp2_volcano_box))
                    tmpdata_volcano_box <- cbind(tmp1_volcano_box,tmp2_volcano_box)
                    colnames(tmpdata_volcano_box) <- c('group','intensity')
                    boxplotdata_volcano_box <- rbind(boxplotdata_volcano_box,tmpdata_volcano_box)
                  }
                  boxplotdata_volcano_box <- boxplotdata_volcano_box[-1,]
                  boxplotdata_volcano_box$intensity <- as.numeric(boxplotdata_volcano_box$intensity)
                  
                  if(length(grep("_",str_trim(tolower(isolate(input$name_volcano_box)))))>0){
                    
                    text_volcano_box <- paste('mz:',strsplit(str_trim(tolower(isolate(input$name_volcano_box))),"_")[[1]][1],' time:',strsplit(str_trim(tolower(isolate(input$name_volcano_box))),"_")[[1]][2],sep='')
                  }else{
                    
                    text_volcano_box <- as.character(volcano_box[tolower(volcano_box$name)==str_trim(tolower(isolate(input$name_volcano_box))),'name'])
                  }
                  
                  if(isolate(input$boxplotcolor_volcano_box)=='yes'){
                    
                    if(isolate(input$xaxis_name_volcano_box_1)==""){
                      xlab_volcano_box <- xlab('group')
                    }else{
                      xlab_volcano_box <- xlab(isolate(input$xaxis_name_volcano_box_1))
                    }
                    if(isolate(input$yaxis_name_volcano_box_1)==""){
                      ylab_volcano_box <- ylab('intensity')
                    }else{
                      ylab_volcano_box <-  ylab(isolate(input$yaxis_name_volcano_box_1))
                    }
                    
                  }else if(isolate(input$boxplotcolor_volcano_box)=='no'){
                    
                    if(isolate(input$xaxis_name_volcano_box_2)==""){
                      xlab_volcano_box <- xlab('group')
                    }else{
                      xlab_volcano_box <- xlab(isolate(input$xaxis_name_volcano_box_2))
                    }
                    if(isolate(input$yaxis_name_volcano_box_2)==""){
                      ylab_volcano_box <- ylab('intensity')
                    }else{
                      ylab_volcano_box <-  ylab(isolate(input$yaxis_name_volcano_box_2))
                    }
                    
                  }
                  
                  
                  if(isolate(input$boxplotcolor_volcano_box)=='yes'){
                    p3 <- ggboxplot(boxplotdata_volcano_box, x = "group", y = "intensity", color = "group", bxp.errorbar = TRUE, bxp.errorbar.width=0.5, size = 0.5) +
                      ggtitle(text_volcano_box) + xlab_volcano_box + ylab_volcano_box +
                      theme(axis.title=element_text(size=12),
                            axis.text =element_text(size=10),
                            axis.text.x = element_text(angle=0),
                            axis.line = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                            plot.title = element_text(size=14, hjust = 0.5))
                  }else{
                    p3 <- ggboxplot(boxplotdata_volcano_box, x = "group", y = "intensity", bxp.errorbar = TRUE, bxp.errorbar.width=0.5, size = 0.5) +
                      ggtitle(text_volcano_box) + xlab_volcano_box + ylab_volcano_box +
                      theme(axis.title=element_text(size=12),
                            axis.text =element_text(size=10),
                            axis.text.x = element_text(angle=0),
                            axis.line = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                            plot.title = element_text(size=14, hjust = 0.5))
                  }
                  
                  closeAlert(session, "check_search_interactive_alert")
                  p3
                  
                }else{
                  
                  closeAlert(session, "check_search_interactive_alert")
                  createAlert(session, "alert_interactive", "check_search_interactive_alert", title = "Search Error", content = paste0("No metabolite with ",str_trim(isolate(input$name_volcano_box))," was found, please check your input file."), append = TRUE)
                  NULL
                }
              
            }else{
              
              closeAlert(session, "check_search_interactive_alert")
              createAlert(session, "alert_interactive", "check_search_interactive_alert", title = "Search Error", content = paste0("No search text was entered."), append = TRUE)
              NULL
            }
          }
        }
        
      }
      
    }
    
    }
  
  })
  
  
  
  pplot1 <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot only' || isolate(input$choose_graph)=='Volcano Plot only'){
      
      if(!is.null(plot1_single_plot())){
        
        ggp1 <- ggplotly(plot1_single_plot(),tooltip = c("text"),source="scatterplot") %>% config(displayModeBar = F)
        ggp1
        
      }
      
    }else{
      
      if(!is.null(plot1_single_plot())){
        
        ggp1 <- ggplotly(plot1_single_plot(),tooltip = c("text"),source="scatterplot") %>% config(displayModeBar = F)
        ggp1
        
      }
      
    }
    
    
  })
  
  
  observeEvent({if(!is.null(event_data('plotly_click',source = 'scatterplot'))) TRUE else return()},{off_interactive(0)})
  observeEvent(input$search_manhattan_box,{off_interactive(1)})
  observeEvent(input$search_volcano_box,{off_interactive(1)})
  
  
  plot4_single_plot <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot'){
      
      if(!is.null(plot3_single_plot()) & off_interactive()==1){
        plot3_single_plot()
      }else if(!is.null(plot2_single_plot()) & off_interactive()==0){
        plot2_single_plot()
      }else{
        NULL
      }
      
    } else if(isolate(input$choose_graph)=='Volcano Plot with Box Plot'){
      
      if(!is.null(plot3_single_plot()) & off_interactive()==1){
        plot3_single_plot()
      }else if(!is.null(plot2_single_plot()) & off_interactive()==0){
        plot2_single_plot()
      }else{
        NULL
      }
      
      
    }
    
  })
  
  
  pplot2 <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot only' || isolate(input$choose_graph)=='Volcano Plot only'){
      
      if(!is.null(plot2_single_plot())){
        
        ggp2 <- ggplotly(plot2_single_plot(),tooltip = c("text"),source="scatterplot") %>% config(displayModeBar = F)
        ggp2
        
      }
      
    }else{
      
      if(!is.null(plot4_single_plot())){
        
        ggp2 <- ggplotly(plot4_single_plot()) %>% config(displayModeBar = F)
        ggp2
      }
      
    }
    
    
  })
  
  
  output$check_plot_manhattan_only <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot only' & !is.null(pplot1()) & !is.null(pplot2())){
      TRUE
    }else{
      FALSE
    }
    
  })
  outputOptions(output, "check_plot_manhattan_only", suspendWhenHidden = FALSE)
  
  
  output$check_plot_volcano_only <- reactive({
    
    if(isolate(input$choose_graph)=='Volcano Plot only' & !is.null(pplot1()) & is.null(pplot2())){
      TRUE
    }else{
      FALSE
    }
    
  })
  outputOptions(output, "check_plot_volcano_only", suspendWhenHidden = FALSE)
  
  
  output$check_plot_manhattan_box1 <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot' & !is.null(pplot1()) ){
      TRUE
    }else{
      FALSE
    }
    
  })
  outputOptions(output, "check_plot_manhattan_box1", suspendWhenHidden = FALSE)
  
  output$check_plot_manhattan_box2_1 <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot' & !is.null(pplot2()) & isolate(input$boxplotcolor_manhattan_box)=='yes'){
      TRUE
    }else{
      FALSE
    }
    
  })
  outputOptions(output, "check_plot_manhattan_box2_1", suspendWhenHidden = FALSE)
  
  output$check_plot_manhattan_box2_2 <- reactive({
    
    if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot' & !is.null(pplot2()) & isolate(input$boxplotcolor_manhattan_box)=='no'){
      TRUE
    }else{
      FALSE
    }
    
  })
  outputOptions(output, "check_plot_manhattan_box2_2", suspendWhenHidden = FALSE)
  
  output$check_plot_volcano_box1 <- reactive({
    
    if(isolate(input$choose_graph)=='Volcano Plot with Box Plot' & !is.null(pplot1()) ){
      TRUE
    }else{
      FALSE
    }
    
  })
  outputOptions(output, "check_plot_volcano_box1", suspendWhenHidden = FALSE)
  
  output$check_plot_volcano_box2_1 <- reactive({
    
    if(isolate(input$choose_graph)=='Volcano Plot with Box Plot' & !is.null(pplot2()) & isolate(input$boxplotcolor_volcano_box)=='yes'){
      TRUE
    }else{
      FALSE
    }
    
  })
  outputOptions(output, "check_plot_volcano_box2_1", suspendWhenHidden = FALSE)
  
  output$check_plot_volcano_box2_2 <- reactive({
    
    if(isolate(input$choose_graph)=='Volcano Plot with Box Plot' & !is.null(pplot2()) & isolate(input$boxplotcolor_volcano_box)=='no'){
      TRUE
    }else{
      FALSE
    }
    
  })
  outputOptions(output, "check_plot_volcano_box2_2", suspendWhenHidden = FALSE)
  
  output$interactive_plot <- renderUI({
    
    if(isolate(input$choose_graph)=='Manhattan Plot only' & !is.null(pplot1()) & !is.null(pplot2())){
      
        if (!is.null(id_interactive)){
          removeNotification(id_interactive)
          id_interactive <<- NULL
        }
      
      output$manhattan_only_plot1 <- renderPlotly({if(!is.null(pplot1())){pplot1()}})
      output$manhattan_only_plot2 <- renderPlotly({if(!is.null(pplot2())){pplot2()}})
      
      column(width=12,
             style='margin-bottom:10px',
             column(width=6,plotlyOutput("manhattan_only_plot1", width = "480px", height = "400px")),
             column(width=6,plotlyOutput("manhattan_only_plot2", width = "480px", height = "400px"))
             
      )
      
    } else if(isolate(input$choose_graph)=='Volcano Plot only' & !is.null(pplot1()) & is.null(pplot2())){
      
      if (!is.null(id_interactive)){
        removeNotification(id_interactive)
        id_interactive <<- NULL
      }
      
      output$volcano_only_plot1 <- renderPlotly({if(!is.null(pplot1())){pplot1()}})
      
      column(width=12,
             style='margin-bottom:10px',
             column(width=12,style='text-align:center',plotlyOutput("volcano_only_plot1",width = "550px", height = "400px"))
             
      )
      
    } else if(isolate(input$choose_graph)=='Manhattan Plot with Box Plot' & !is.null(pplot1()) ){
      
      if (!is.null(id_interactive)){
        removeNotification(id_interactive)
        id_interactive <<- NULL
      }
      
      output$manhattan_box_plot1 <- renderPlotly({if(!is.null(pplot1())){pplot1()}})
      output$manhattan_box_plot2 <- renderPlotly({if(!is.null(pplot2())){pplot2()}})
      
      column(width=12,
             style='margin-bottom:10px',
             column(width=6,plotlyOutput("manhattan_box_plot1", width = "480px", height = "400px")),
             column(width=6,plotlyOutput("manhattan_box_plot2", width = "500px", height = "400px"))
      
      )
    
    } else if(isolate(input$choose_graph)=='Volcano Plot with Box Plot' & !is.null(pplot1())){
    
      if (!is.null(id_interactive)){
        removeNotification(id_interactive)
        id_interactive <<- NULL
      }
      
      output$volcano_box_plot1 <- renderPlotly({if(!is.null(pplot1())){pplot1()}})
      output$volcano_box_plot2 <- renderPlotly({if(!is.null(pplot2())){pplot2()}})
      
      column(width=12,
             style='margin-bottom:10px',
             column(width=6,plotlyOutput("volcano_box_plot1", width = "480px", height = "400px")),
             column(width=6,plotlyOutput("volcano_box_plot2", width = "500px", height = "400px"))
             
      )
      
  }
    
    
  })
  
  
  output$downloadPlot1_manhattan_only <- downloadHandler(
    
    filename <- "type1_manhattan_plot.png",
    content = function(file) {
      ggsave(file, plot = plot1_single_plot(), device = "png", width=4, height=4)
    }
  )
  
  output$downloadPlot2_manhattan_only <- downloadHandler(
    
    filename <- "type2_manhattan_plot.png",
    content = function(file) {
      ggsave(file, plot = plot2_single_plot(), device = "png", width=4, height=4)
    }
  )
  
  output$downloadPlot1_volcano_vonly <- downloadHandler(
    
    filename <- "volcano_plot.png",
    content = function(file) {
      ggsave(file, plot = plot1_single_plot(), device = "png", width=4, height=4)
    }
  )
  
  output$downloadPlot1_manhattan_box <- downloadHandler(
    
    filename <- "type_manhattan_plot.png",
    content = function(file) {
      ggsave(file, plot = plot1_single_plot(), device = "png", width=4, height=4)
    }
  )
  
  output$downloadPlot2_1_manhattan_box <- downloadHandler(
    
    filename <- "box_plot.png",
    content = function(file) {
      ggsave(file, plot = plot4_single_plot(), device = "png", width=4, height=4)
    }
  )
  
  output$downloadPlot2_2_manhattan_box <- downloadHandler(
    
    filename <- "box_plot.png",
    content = function(file) {
      ggsave(file, plot = plot4_single_plot(), device = "png", width=4, height=4)
    }
  )
  
  output$downloadPlot1_volcano_box <- downloadHandler(
    
    filename <- "type_volcano_plot.png",
    content = function(file) {
      ggsave(file, plot = plot1_single_plot(), device = "png", width=4, height=4)
    }
  )
  
  output$downloadPlot2_1_volcano_box <- downloadHandler(
    
    filename <- "box_plot.png",
    content = function(file) {
      ggsave(file, plot = plot4_single_plot(), device = "png", width=4, height=4)
    }
  )
  
  output$downloadPlot2_2_volcano_box <- downloadHandler(
    
    filename <- "box_plot.png",
    content = function(file) {
      ggsave(file, plot = plot4_single_plot(), device = "png", width=4, height=4)
    }
  )
  
  #################### interactive plot end
  
} 


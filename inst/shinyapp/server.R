options(shiny.maxRequestSize=100*1024^2)
options(shiny.sanitize.errors=FALSE)
suppressPackageStartupMessages(library('xmsPANDA'))

suppressMessages(require(shiny))
                 suppressMessages(require(shinyjs))
                                  suppressMessages(require(shinyBS))
                                                   suppressMessages(require(DT))

                                                   suppressMessages(source("R/source_codes/xmsPANDA_v1.0.9.39.R"))

# Server logic

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
      sum(input$featselmethodi%in%c('limma','limma2way','limmarobust','limma2wayrobust','lm1wayanova','lm2wayanova','lmreg','logitreg','rfesvm','ttest','wilcox','RF','MARS','pamr'))>0
    }else{
      if(input$analysismode == 'classification' && input$pairedanalysis == 'TRUE'){
        sum(input$featselmethodii%in%c('limma1wayrepeat','limma2wayrepeat','limma1wayrepeatrobust','limma2wayrepeatrobust','lm1wayanovarepeat','lm2wayanovarepeat','ttestrepeat','wilcoxrepeat'))>0
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
  
  session_outloc <- eventReactive(input$go,{
    # if(input$go!=0 & check$count==1){
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
    #}else{
    #outloc<-paste('~/xmsPANDAout',sep="")
    #outloc
    #}
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
    
 
    
  
    
    checkellipse_conf_levelAlert <- TRUE
    
    all(checknumreplicateAlert, checksummarization_ratioAlert, checkall_missing_threshAlert, checkrsd_filt_listAlert, checkgroup_missing_threshAlert,
        checkpvalue_threshAlert, checkfdrthreshAlert, checkfoldchangethreshAlert, checkkfoldAlert, checkpls_vip_threshAlert,
        checkpls_ncompAlert, checkmax_comp_selAlert, checkpls_permut_countAlert, checkmax_varselAlert, checkabs_cor_threshAlert,checkcor_fdrthreshAlert,
        checkellipse_conf_levelAlert)
    
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
    req(check$count == 0)
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
                 
                 
                 # reset("nText2")
                 #reset("nText")
                 #reset("id1")
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
                 id1 <<- showNotification("Starting processing now. Your results will be available for download shortly. The processing time depends on the number of methods used.", duration=NULL)
                 
                 output$output_results <- renderUI({})
                 
               })
  
  
  
  #########################################
  
  
  featselmethod_check <-eventReactive(input$go,{
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
  
  
  
  
  featuretable <- eventReactive(input$go,{
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
  
  classlabel <- eventReactive(input$go,{
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
  
  
  
  
  
  ##########################################
  
  output$nText <- renderPrint({
    if(input$go!=0  & check$count==1 & !is.null(featselmethod_check()) & is.data.frame(featuretable()) & is.data.frame(classlabel()) & all_alert()==TRUE){
      
      
      
      
      
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
      
      if(input$hca.cex.legend==(-1)){
        input$hca.cex.legend=NA
        
      }
      
      if(input$ggplot.type1=="None"){
        
        ggplot.type1.val=NA
      }else{
        
        ggplot.type1.val=input$ggplot.type1
      }
      
      #"redblue","orangeblue","bluered","blueorange"
      if(input$heatmap_color_scheme=="redblue"){
        manhattanplot.col.opt=c("red3","darkblue")
        
      }else{
        if(input$heatmap_color_scheme=="bluered"){
          manhattanplot.col.opt=c("darkblue","red3")
          
        }else{
          
          if(input$heatmap_color_scheme=="blueorange"){
            manhattanplot.col.opt=c("darkblue","orange")
            
          }else{
            
            if(input$heatmap_color_scheme=="orangeblue"){
              manhattanplot.col.opt=c("orange","darkblue")
              
            }
          }
        }
        
      }
      
      # if(input$timeseries.lineplots==FALSE){
      
      #   timeseries.lineplots=FALSE
      #}
      
      ###################
      msg=""
      # if(input$workflow=='workflowI')
      # if(check$count==1)
      if(input$go)
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
          pca.ellipse=pca_ellipse,ellipse.conf.level=0.95,svm.acc.tolerance=5,pamr.threshold.select.max=FALSE,
          aggregation.method=aggregation_method,mars.gcv.thresh=1,pls.vip.selection=input$pls.vip.selection,limmadecideTests=input$limmadecideTests,lme.modeltype=input$lme.modeltype,
          
          #4) arguments for WGCNA and global clustering analysis (HCA and EM clustering)
          wgcnarsdthresh=input$rsd_filt_list,WGCNAmodules=WGCNAmodules,globalclustering=globalclustering,
          
          #5) arguments for correlation and network analysis using the selected features
          cor.method=cor_method, abs.cor.thresh = abs_cor_thresh, cor.fdrthresh=cor_fdrthresh,
          globalcor=globalcor,target.metab.file=NA,
          target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,
          degree_rank_method=degree_rank_method, #set to NA or "DiffRank"
          
          #6) arguments for graphical options: see manual for additional arguments
          output.device.type="png",pca.cex.val=6,legendlocation="bottomleft",
          net_node_colors=c("green","red"),net_legend=FALSE,
          manhattanplot.col.opt=manhattanplot.col.opt, #c("darkblue","red3"),
          heatmap.col.opt=input$heatmap_color_scheme,
          
          aggregation.max.iter=100,
          plots.width=input$plots.width,
          plots.height=input$plots.height,
          plots.type="cairo",
          cex.plots=input$cex.plots,
          
          add.pvalues=input$boxplot.pvalues,
          add.jitter=input$boxplot.jitter,
          timeseries.lineplots=input$timeseries.lineplots,
          alphabetical.order=input$alphabetical.order,
          ylab_text=input$ylabel.text,
          kegg_species_code="hsa",database="pathway",reference_set=NA,match_class_dist=TRUE,
          differential.network.analysis=differential.network.analysis,
          log2.transform.constant=input$log2.transform.constant,
          balance.classes=input$balance.classes,balance.classes.sizefactor=10,
          balance.classes.seed=1,cv.perm.count=100,multiple.figures.perpanel=FALSE,
          hca.labRow.value = input$hca.labRow.value, hca.labCol.value = input$hca.labCol.value,
          alpha.col=1,
          similarity.matrix=input$similarity.matrix,outlier.method=input$outlier.method,removeRda=TRUE,
          color.palette=input$color.palette,plot_DiNa_graph=FALSE,limma.contrasts.type=input$limma.contrasts.type,
          hca.cex.legend=input$hca.cex.legend,plot.boxplots.raw=input$plot.boxplots.raw,vcovHC.type=input$vcovHC.type,ggplot.type1=ggplot.type1.val,
          facet.nrow=input$facet.nrow
          
        ),silent=TRUE)
        
        check$count=0
        if(is(demetabs_res,"try-error")){
          done$count=1
          go <- reactiveValues(count = 0)
          # done <- reactiveValues(count = 0)
          #file.copy(paste(getwd(),'matrix_centrality.txt',sep='/'),session_outloc())
          #  print("Error processing the data. Error message:")
          #print(demetabs_res)
          #setwd(session_outloc())
          #zip(zipfile=paste(basename(session_outloc()),'zip',sep='.'), files='.')
          msg=paste("Error processing the data. Click on download button to save the partial results and review the Log.txt file. Error message:", demetabs_res)
          
          observeEvent(req(done$count==1),{
            if (!is.null(id1)){
              removeNotification(id1)
              id1 <<- showNotification(msg,duration=60)
            }
          })
          
        }else{
          #print(session_outloc())
          done$count=1
          go <- reactiveValues(count = 0)
          setwd(session_outloc())
          zip(zipfile=paste(basename(session_outloc()),'zip',sep='.'), files='.')
          print("Processing complete. Please click on download button to save the results.")
          #msg="Processing complete. Please click on download button to save the results."
          
          observeEvent(req(done$count==1),{
            if (!is.null(id1)){
              removeNotification(id1)
              id1 <<- showNotification("Processing complete. Please click on download button to save the results.",duration=60)
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
          
        }
        
      }
      
      #   input$go=0
      
      go <- reactiveValues(count = 0)
      # go <- reactiveValues(count = 0)
      #  suppressWarnings(reset("go"))
      
      #   done <- reactiveValues(count = 0)
      #go <- reactiveValues(count = 0)
      
      
      
      # }else{
      
      #NULL
    }
    
    #msg
  })
  
  ##########################################
  
  
  
  observeEvent(input$methodout,{
    
    if(input$methodout=="AggregatedResults"){
      l1 <- list.files(paste(session_outloc(),'AggregatedResults',sep="/"),".png",recursive=TRUE,full.names=FALSE)
      figurenum <- paste('Figure',seq(1:length(l1)))
    }else{
      folder <- grep(paste(input$methodout,"signalthresh",sep=""),list.dirs(paste(session_outloc(),'Stage2',sep='/'),recursive=FALSE,full.names=FALSE),value=TRUE)
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
  
  
  #observeEvent(input$normstart,
  eventReactive(input$normstart,
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
  
  #session_outloc2 <- reactive({
  session_outloc2 <- eventReactive(input$normstart,{
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
      normcheck2$count=0
      setwd(session_outloc2())
      zip(zipfile=paste(basename(session_outloc2()),'zip',sep='.'), files='.')
      print("Processing complete. Please click on download button to save the results.")
      
      # reset("normstart")
      
      normstart <- reactiveValues(count = 0)
      
    }else{
      
      NULL
    }
    
  })
  
  observeEvent(req(normdone2$count==1),{
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
  
  
  ###start pca
  pcacheck2 <- reactiveValues(count = 0)
  pcadone2 <- reactiveValues(count = 0)
  #normdone2 <- reactiveValues(res = list())
  #update2 <- reactiveValues(count = 0)
  pcaid3 <-NULL
  
  # observeEvent(input$pcastart,{pcacheck2$count=0})
  eventReactive(input$pcastart,{pcacheck2$count=0})
  
  #observeEvent
  #eventReactive
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
  
  session_outloc3 <- eventReactive(input$pcastart,{
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
  
  
  output$pcaText5 <- renderPrint({
    
    if(input$pcastart!=0  & pcacheck2$count==1){
      
      
      pca_ellipse_val<-(input$pca.ellipse)
     # showNotification(pca_ellipse_val)
      if(pca_ellipse_val=="TRUE"){
      pcares<-get_pcascoredistplots(X=pca_metab_data(),Y=pca_class_data(),feature_table_file=NA,parentoutput_dir=session_outloc3(),class_labels_file=NA,sample.col.opt=input$pca.sample.color.theme,
                                    plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=NA,pairedanalysis=input$pca.pairedanalysis,pca.cex.val=2,legendlocation="topright",
                                    pca.ellipse=TRUE,ellipse.conf.level=0.95,filename="PCA_results",paireddesign=NA,newdevice=TRUE,
                                    lineplot.col.opt=input$pca.sample.color.theme,lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                                    timeseries.lineplots=input$pca.timeseries.lineplots,pcacenter=TRUE,pcascale=TRUE,alphabetical.order=input$pca.alphabetical.order,analysistype=input$pca.analysistype,lme.modeltype=lme.modeltype) #,silent=TRUE)
      
      }else{
        
        pcares<-get_pcascoredistplots(X=pca_metab_data(),Y=pca_class_data(),feature_table_file=NA,parentoutput_dir=session_outloc3(),class_labels_file=NA,sample.col.opt=input$pca.sample.color.theme,
                                      plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=NA,pairedanalysis=input$pca.pairedanalysis,pca.cex.val=2,legendlocation="topright",
                                      pca.ellipse=FALSE,ellipse.conf.level=0.95,filename="PCA_results",paireddesign=NA,newdevice=TRUE,
                                      lineplot.col.opt=input$pca.sample.color.theme,lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                                      timeseries.lineplots=input$pca.timeseries.lineplots,pcacenter=TRUE,pcascale=TRUE,
                                      alphabetical.order=input$pca.alphabetical.order,analysistype=input$pca.analysistype,lme.modeltype=lme.modeltype) #,silent=TRUE)
        
      }
      
      pcadone2$count=1
      pcacheck2$count=0
      setwd(session_outloc3())
      
      
      zip(zipfile=paste(basename(session_outloc3()),'zip',sep='.'), files='.')
      
      print("Processing complete. Please click on download button to save the results.")
      
      pcacheck2$count=0
      
      #reset("pcastart")
      
      pcastart <- reactiveValues(count = 0)
      
      
    }
    
  })
  
  
  observeEvent(req(pcadone2$count==1),{
    if (!is.null(pcaid3)){
      removeNotification(pcaid3)
      pcaid3 <<- showNotification("Processing complete. Click on Download Results to save the results.", duration=NULL)
      
      
    }
    
    
    # column(12,style='padding-top:10px;padding-left:0;',tags$div(h4("Output")))
    output$pcaoutput_results <- renderText({
      
      pcal1 <- list.files(paste(session_outloc3()),".pdf",recursive=FALSE,full.names=FALSE)
      
      
      
      filename <- normalizePath(file.path(session_outloc3(),pcal1))
      
      
      PDFfile=filename
      print(paste("file exists:",file.exists(PDFfile)))
      list(src=PDFfile)
      
      
      
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
  
  
  #####Boxplots
  
  
  ###start boxplot
  boxplotcheck2 <- reactiveValues(count = 0)
  boxplotdone2 <- reactiveValues(count = 0)
  #normdone2 <- reactiveValues(res = list())
  #update2 <- reactiveValues(count = 0)
  boxplotid3 <-NULL
  
  # observeEvent(input$boxplotstart,{boxplotcheck2$count=0})
  eventReactive(input$boxplotstart,{boxplotcheck2$count=0})
  
 #eventReactive(input$boxplot.ylabel.text,{boxplot_ylabel=renderText(input$boxplot.ylabel.text)})
  #observeEvent
  #eventReactive
  observeEvent(input$boxplotstart,
               {
                 output$boxplotText4 <- renderText({shiny::validate(
                   need(input$boxplotinput1, "No data file was provided. Please upload your data file."),
                   need(input$boxplotinput1$type=="text/csv" || input$boxplotinput1$type=="text/plain", "The format of data file is not correct. Please upload the file with correct format.")
                 )})
                 shiny::validate(
                   need(input$boxplotinput1, "No data file was provided. Please upload your data file."),
                   need(input$boxplotinput1$type=="text/csv" || input$boxplotinput1$type=="text/plain", "The format of data file is not correct. Please upload the file with correct format.")
                 )
                 boxplotcheck2$count=1
                 boxplotid3 <<- showNotification("Data is processing now.", duration=NULL)
               })
  
  
  boxplot_metab_data <- reactive({
    
    if(input$boxplotstart!=0  & boxplotcheck2$count==1){
      
      if(input$boxplotinput1$type=="text/plain"){
        metab_data <- read.delim(input$boxplotinput1$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        
      }else{
        if(input$boxplotinput1$type=="text/csv"){
          metab_data <- read.csv(input$boxplotinput1$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
          
        }else{
          
          metab_data <- NULL
          
        }
      }
      metab_data
    }
    
  })
  
  boxplot_class_data <- reactive({
    
    if(input$boxplotstart!=0  & boxplotcheck2$count==1){
      
      if(input$boxplotinput2$type=="text/plain"){
        
        class_data <- read.delim(input$boxplotinput2$datapath,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
      }else{
        if(input$boxplotinput2$type=="text/csv"){
          
          class_data <- read.csv(input$boxplotinput2$datapath,header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)
        }else{
          
          
          class_data<-NULL
        }
      }
      class_data
    }
    
  })
  
  session_boxplotoutloc3 <- eventReactive(input$boxplotstart,{
    if(input$boxplotstart!=0  & boxplotcheck2$count==1){
      cur_date<-Sys.time()
      cur_date<-gsub(x=cur_date,pattern="-",replacement="")
      cur_date<-gsub(x=cur_date,pattern=":",replacement="")
      cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
      
      boxplotoutloc<-paste('~/xmsPANDAboxplots',cur_date,sep="")
      
      boxplotoutloc
    }else{
      cur_date<-Sys.time()
      cur_date<-gsub(x=cur_date,pattern="-",replacement="")
      cur_date<-gsub(x=cur_date,pattern=":",replacement="")
      cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
      boxplotoutloc<-paste('~/xmsPANDAboxplots',cur_date,sep="")
      boxplotoutloc
    }
    
    # normid3 <<- showNotification("Data is starting now.", duration=NULL)
  })
  
  
  output$boxplotText5 <- renderPrint({
    
    if(input$boxplotstart!=0  & boxplotcheck2$count==1){
      
      
      
      boxplot_ylab=(input$boxplot.ylabel.text)
      
      #showNotification(input$boxplot_jitter_1)
      #showNotification(input$boxplot_pvalues_1)
      
      if(input$boxplot.min.ylim==(-1) && input$boxplot.max.ylim==(-1)){
      boxplotres<-get_boxplots(X=boxplot_metab_data(),Y=boxplot_class_data(),feature_table_file=NA,parentoutput_dir=session_boxplotoutloc3(),class_labels_file=NA,
                             boxplot.col.opt=input$boxplot.sample.color.theme,
                             alphacol=1,newdevice=TRUE,cex.plots=0.8,replace.by.NA=FALSE,pairedanalysis=FALSE,
                             filename="boxplots",
                             ylabel=boxplot_ylab,
                             alphabetical.order=input$boxplot.alphabetical.order,name=NA,
                             add.jitter=input$boxplot_jitter_1,add.pvalues=input$boxplot_pvalues_1,class.levels=NA,fill.plots=TRUE,
                             connectpairedsamples=FALSE,boxplot.type=input$boxplot.type,
                             study.design=input$boxplot.analysistype,
                             multiple.figures.perpanel=FALSE,ggplot.type1=input$boxplot.ggplot.type1,replace.outliers=FALSE,
                             plot.height=input$boxplot.plots.height,plot.width=input$boxplot.plots.width,
                             extra_text=NA,group_by_mat=NA,position_dodge_width=0.75,
                             numnodes=2,hightlight.points=FALSE,ref.group.val=FALSE,facet.nrow=1)
      }else{
        boxplotres<-get_boxplots(X=boxplot_metab_data(),Y=boxplot_class_data(),feature_table_file=NA,parentoutput_dir=session_boxplotoutloc3(),class_labels_file=NA,
                                 boxplot.col.opt=input$boxplot.sample.color.theme,
                                 alphacol=1,newdevice=TRUE,cex.plots=0.8,replace.by.NA=FALSE,pairedanalysis=FALSE,
                                 filename="boxplots",
                                 ylabel=boxplot_ylab,
                                 alphabetical.order=input$boxplot.alphabetical.order,name=NA,
                                 add.jitter=input$boxplot_jitter_1,add.pvalues=input$boxplot_pvalues_1,class.levels=NA,fill.plots=TRUE,
                                 connectpairedsamples=FALSE,boxplot.type=input$boxplot.type,
                                 study.design=input$boxplot.analysistype,
                                 multiple.figures.perpanel=FALSE,ggplot.type1=input$boxplot.ggplot.type1,replace.outliers=FALSE,
                                 plot.height=input$boxplot.plots.height,plot.width=input$boxplot.plots.width,
                                 extra_text=NA,group_by_mat=NA,position_dodge_width=0.75,
                                 numnodes=2,hightlight.points=FALSE,ref.group.val=FALSE,facet.nrow=1,ylim.val=c(input$boxplot.min.ylim,input$boxplot.max.ylim))
        
      }
      boxplotdone2$count=1
      boxplotcheck2$count=0
      setwd(session_boxplotoutloc3())
      
      
      zip(zipfile=paste(basename(session_boxplotoutloc3()),'zip',sep='.'), files='.')
      
      print("Processing complete. Please click on download button to save the results.")
      
      boxplotcheck2$count=0
      
      #reset("boxplotstart")
      
      boxplotstart <- reactiveValues(count = 0)
      
      
    }
    
  })
  
  
  observeEvent(req(boxplotdone2$count==1),{
    if (!is.null(boxplotid3)){
      removeNotification(boxplotid3)
      boxplotid3 <<- showNotification("Processing complete. Click on Download Results to save the results.", duration=NULL)
      
      
    }
    
    
    # column(12,style='padding-top:10px;padding-left:0;',tags$div(h4("Output")))
    output$boxplotoutput_results <- renderText({
      
      boxplotl1 <- list.files(paste(session_boxplotoutloc3()),".pdf",recursive=FALSE,full.names=FALSE)
      
      
      
      filename <- normalizePath(file.path(session_boxplotoutloc3(),boxplotl1))
      
      
      PDFfile=filename
      print(paste("file exists:",file.exists(PDFfile)))
      list(src=PDFfile)
      
      
      
    }) #,deleteFile=FALSE)
    
  })
  
  
  
  output$downloadboxplotData <- downloadHandler(
    
    #if(input$go!=0 && input$featselmethod!="-" && input$feature_table_file!="" && input$class_labels_file!=""){
    
    filename <- function() {
      paste(basename(session_boxplotoutloc3()), "zip", sep=".")
    },
    content <- function(file) {
      fname1<-paste(session_boxplotoutloc3(),"/",basename(session_boxplotoutloc3()), ".zip", sep="")
      file.copy(fname1, file)
    },
    contentType = "application/zip"
    #}
  )
  
  
  
  
  ##########################
  
  
  ####Functional Class Scoring
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
  
  #eventReactive(input$start2,{check2$count=0})
  
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
                 id3 <<- showNotification(paste("Data is processing now using the ",input$fcs.database,"  database",sep=""), duration=NULL)
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
  
  output$nText5A <- renderPrint({
    
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
        done2$cluster_table <-get_fcs(target.data=cluster_metab_data(),target.data.annot=target.data.annot,
                                      kegg_species_code=kegg_species_code,database="custom",type.statistic=type.statistic,
                                      fcs.min.hits=input$fcs.min.hits,reference_set=cluster_metab_data2(),itrs=input$fcs.itrs)
        
      }else{
        done2$cluster_table <-get_fcs(target.data=cluster_metab_data(),target.data.annot=target.data.annot,kegg_species_code=kegg_species_code,database=input$fcs.database,type.statistic=type.statistic,fcs.min.hits=input$fcs.min.hits,reference_set=NA,itrs=input$fcs.itrs)
      }
      check2$count=0
      
      NULL
      #checktable1=1
      #print(dim(done2$cluster_table))
      #print("yes2")
      #  reset("start2")
      
      start2 <- reactiveValues(count = 0)
      
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
  
  
  observeEvent(req(dim(done2$cluster_table)[2]>1),{
    # eventReactive({if(dim(done2$cluster_table)[2]>1) TRUE else return()},{
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
        }else{
          
          if(input$fcs.database=="kegg_atlas"){
            
            map_cpd<-paste(table1[,1],"+",as.character(table1$XID),sep="")
            
            temp_pattern=paste(" ",input$path.bg.color,"+",sep="")
            map_cpd<-gsub(as.character(map_cpd),pattern=";",replacement=temp_pattern)
            map_cpd<-paste(map_cpd," ",input$path.bg.color,sep="")
            
            
            table1[,1] <- paste("<a target='_blank' href='https://www.genome.jp/kegg-bin/show_pathway?",map_cpd,"'>",table1[,1],"</a>",sep="")
          }
          
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
          if(isolate(input$fcs.database)=='kegg_atlas'){
            table[i,1]=paste("<a target='_blank' href='https://www.genome.jp/kegg-bin/show_pathway?",table[i,1],"'>",table[i,1],"</a>",sep="")
          }
          
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
               #eventReactive(input$start3,
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
    #if(input$start3!=0 & check3$count==1){
    cur_date<-Sys.time()
    cur_date<-gsub(x=cur_date,pattern="-",replacement="")
    cur_date<-gsub(x=cur_date,pattern=":",replacement="")
    cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
    outloc<-paste('~/metabolite_quantification_analysis_results',cur_date,sep="")
    outloc
    #  }else{
    #   NULL
    #  }
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
        check3$count=0
        #file.copy(paste(getwd(),'matrix_centrality.txt',sep='/'),session_outloc())
        setwd(session_outloc_quant())
        zip(zipfile=paste(basename(session_outloc_quant()),'zip',sep='.'), files='.')
        print("Processing complete. Please click on download button to save the results.")
        #reset("start3")
        
        start3 <- reactiveValues(count = 0)
        
      })
      
    }else{
      
      NULL
    }
    
  })
  
  
  observeEvent(req(done3$count==1),{
    if (!is.null(id4)){
      removeNotification(id4)
      id4 <<- showNotification(paste("Processing complete. Click on download button to save the results. Output location: ",session_outloc_quant(),sep=""), duration=NULL)
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
  
 # example_feat <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/feature_table_one_or_two_factor_analysis.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_feat <- read.delim("example_data/feature_table_one_or_two_factor_analysis.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  
  example_feat <- example_feat[1:5,1:7]
  colnames(example_feat) <- c(colnames(example_feat)[1:6],"...")
  output$example_feat <- renderTable({ example_feat }, striped = TRUE)
  
  #example_feat2 <- read.delim("https://raw.githubusercontent.com/kuppal2/xmsPANDA/master/inst/shinyapp/example_data/feature_table_one_or_two_factor_analysis_withNamesorIDs.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_feat2 <- read.delim("example_data/feature_table_one_or_two_factor_analysis_withNamesorIDs.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  
   example_feat2 <- example_feat2[1:5,1:7]
  colnames(example_feat2) <- c(colnames(example_feat2)[1:6],"...")
  output$example_feat2 <- renderTable({ example_feat2 }, striped = TRUE)
  
  example_multiclass_comparison <- read.delim("example_data/classlabels_multiclass_comparison.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_multiclass_comparison <- rbind(head(example_multiclass_comparison[example_multiclass_comparison$Factor1=="Group1",],8),head(example_multiclass_comparison[example_multiclass_comparison$Factor1=="Group2",],8))
  example_multiclass_comparison_covariates <- read.delim("example_data/classlabels_with_covariates.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_multiclass_comparison_covariates <- rbind(head(example_multiclass_comparison_covariates[example_multiclass_comparison_covariates$Class=="NonSmoker",],8),head(example_multiclass_comparison_covariates[example_multiclass_comparison_covariates$Class=="Smoker",],8))
  example_regression <- read.delim("example_data/classlabels_regression.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_regression <- example_regression[1:16,]
  example_two_way_anova <- read.delim("example_data/classlabels_two_way_anova.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_two_way_anova <- rbind(head(example_two_way_anova[example_two_way_anova$Factor1=="Group1",],8),head(example_two_way_anova[example_two_way_anova$Factor1=="Group2",],8))
  example_one_factor_repeatedmeasures <- read.delim("example_data/classlabels_one_factor_repeatedmeasures.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  example_one_factor_repeatedmeasures <- example_one_factor_repeatedmeasures[1:16,]
  example_two_factor_repeatedmeasures <- read.delim("example_data/classlabels_two_factor_repeatedmeasures.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
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
  
  
  
  #################### interactive plot start
  
  
  #################### interactive plot end
  
}


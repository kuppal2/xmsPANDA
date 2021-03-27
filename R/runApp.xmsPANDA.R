runApp.xmsPANDA <-
function() {
    
    suppressMessages(library(shiny))
    suppressMessages(library(shinyjs))
    suppressMessages(library(shinyBS))
    suppressMessages(library(DT))
    
    appDir <- system.file("shinyapp", package = "xmsPANDA")
    if (appDir == "") {
      stop("Could not find shinyapp directory. Try re-installing `xmsPANDA`.", call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal")
  }

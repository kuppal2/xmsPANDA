runApp.xmsPANDA <-
function() {
    
<<<<<<< HEAD
    library(shiny)
    library(shinyjs)
    library(shinyBS)
    library(DT)
    
    appDir <- system.file("shinyapp", package = "xmsPANDA")
    if (appDir == "") {
      stop("Could not find shinyapp directory. Try re-installing `xmsPANDA`.", call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal")
  }
=======
      library(shiny)
      library(shinyjs)
      library(shinyBS)
      library(DT)
    
    appDir <- system.file("shinyapp", package = "xmsPANDA")
    if (appDir == "") {
        stop("Could not find shinyapp directory. Try re-installing `xmsPANDA`.", call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal")
}
>>>>>>> 4ccdcb99e71707b6d2e6cfcfae418ec4bdb9aae3

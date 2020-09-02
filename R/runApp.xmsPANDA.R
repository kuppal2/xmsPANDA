runApp.xmsPANDA <-
function() {
<<<<<<< HEAD
    
      library(shiny)
      library(shinyjs)
      library(shinyBS)
      library(DT)
    
=======
 
    library(shiny)
    library(shinyjs)
    library(shinyBS)
    library(DT)
>>>>>>> c3ff66c6826817a36eed061db658de2fb3145900
    appDir <- system.file("shinyapp", package = "xmsPANDA")
    if (appDir == "") {
        stop("Could not find shinyapp directory. Try re-installing `xmsPANDA`.", call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal")
}

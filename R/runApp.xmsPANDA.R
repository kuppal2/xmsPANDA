runApp.xmsPANDA <-
function() {
<<<<<<< HEAD
    library(shiny)
    library(shinyBS)
=======
    library(shinyBS)
    library(shiny)
>>>>>>> 3c4f9d127a1a49d34d635ca6d7a0fdad51239eb6
    appDir <- system.file("shinyapp", package = "xmsPANDA")
    if (appDir == "") {
        stop("Could not find shinyapp directory. Try re-installing `xmsPANDA`.", call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal")
}

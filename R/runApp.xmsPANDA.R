runApp.xmsPANDA <-
function() {
    library(shiny)
    library(shinyBS)
    appDir <- system.file("shinyapp", package = "xmsPANDA")
    if (appDir == "") {
        stop("Could not find shinyapp directory. Try re-installing `xmsPANDA`.", call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal")
}

.get_data <-
function(fit, data = NULL, complain = FALSE) {
  if(is.null(data)){
    if (complain)
      warning ("The `data` argument is not provided. Data will be extracted from model fit.")
    data <- eval(fit$call$data)
  }
  data
}
.onAttach <-
function(libname, pkgname) {
  # to show a startup message
  #packageStartupMessage("xmsPANDA v1.2 successfully loaded.")
  
  suppressMessages(library(RColorBrewer))
  #suppressMessages(library(data.table))
  suppressMessages(library(plyr))
  suppressMessages(library(parallel))
  suppressMessages(library(CMA))
  suppressMessages(library(ggplot2))
  #suppressMessages(library(extrafont))
  
}
.onLoad <-
function(libname, pkgname) {
  # something to run
  packageStartupMessage("xmsPANDA v1.3.1 (03/29/2021) successfully loaded.")
  suppressMessages(library(RColorBrewer))
  #suppressMessages(library(data.table))
  suppressMessages(library(plyr))
  suppressMessages(library(parallel))
  suppressMessages(library(CMA))
  suppressMessages(library(ggplot2))
  # suppressMessages(library(extrafont))
}

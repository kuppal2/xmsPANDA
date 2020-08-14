\name{get_volcanoplots}
\alias{get_volcanoplots}
\title{get_volcanoplots}
\description{
Function to generate volcano plots
}
\usage{
get_volcanoplots(xvec, yvec, up_or_down, maintext = "", ythresh = 0.05,
 y2thresh = NA, ylab, xlab, colorvec = c("darkblue", "red3"), 
col_seq = c("brown", "chocolate3", "orange3", "coral", "pink",
 "skyblue", "blue", "darkblue", "purple", "violet"), 
xincrement = 1, yincrement = 1, xthresh = 1, pchvec = c(21, 21),
 background.points.col = "gray50", bad.feature.index = NA)
}
\arguments{
  \item{xvec}{
 Vector with log2 fold change values for
  the x-axis.
}
  \item{yvec}{
 Vector with values for the y-axis.
  eg: (-)log10 of p-values, VIP, loadings,
  regression coefficients.
}
  \item{up_or_down}{
  Same as xvec.
}
  \item{maintext}{
Text for the plot title
}
  \item{ythresh}{
Y-axis threshold for significance or differential
expression. eg: 3 for p=0.001; y=(-1)*log10(0.001)
or 2 for VIP from PLS
}
\item{y2thresh}{
Optional secondary Y-axis threshold for significance or differential
expression. eg: 1.3 for p=0.05; y=(-1)*log10(0.05)
or 1 for VIP from PLS
}
  \item{ylab}{
Y-axis label
}
  \item{xlab}{
X-axis label
}
  \item{colorvec}{
Vector of colors for representing up-regulation 
and down-regulation. eg: c("darkgreen", "firebrick1")
In this case, features that are up-regulated in class A
will have "darkgreen" color, and features that are 
up-regulated in class B will have "firebrick1" color.

}
  \item{xincrement}{
Window size for breaking the x-axis into different
segments for visualization purposes. eg: 150
}
  \item{yincrement}{
Window size for breaking the y-axis into different
segments for visualization purposes. eg: 1
}
  \item{xthresh}{
Absolute value of the threshold for log2 fold
change. e.g. 0, 1, 2
}

}
\details{
This function generates volcano plots. 
}
\value{
Volcano plot
}
\author{
Karan Uppal; kuppal2@emory.edu
}

\name{get_manhattanplots}
\alias{get_manhattanplots}
\title{get_manhattanplots}
\description{
Function to generate Manhattan plots.
}
\usage{
get_manhattanplots(xvec, yvec, up_or_down, maintext = "", ythresh = 0.05,
y2thresh=NA, ylab, xlab, colorvec = c("darkgreen", "firebrick1"),
col_seq = c("brown", "chocolate3", "orange3", "coral", "pink", "skyblue",
"blue", "darkblue", "purple", "violet"), xincrement = 150, yincrement = 1)
}
\arguments{
  \item{xvec}{
  Vector with values for the x-axis.
  eg: m/z or retention time values
}
  \item{yvec}{
 Vector with values for the y-axis.
  eg: (-)Log10 of p-values, VIP, loadings,
  regression coefficients, etc. 
}
  \item{up_or_down}{
Vector indicating directionality of change. 
eg: Fold change values
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
Secondary Y-axis threshold for significance or differential
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
  \item{col_seq}{
Vector of colors for plotting different segments of the
x-axis
}
  \item{xincrement}{
Window size for breaking the x-axis into different
segments for visualization purposes. eg: 150
}
  \item{yincrement}{
Window size for breaking the y-axis into different
segments for visualization purposes. eg: 1
}
}
\details{
This function can be used to generate
Type 1 Manhattan plots: signficance vs m/z
Type 2 Manhattan plots: significance vs retention time
Type 3 Manhattan plots: significance vs intensity
}
\value{
Manhattan plots
}
\author{
Karan Uppal
}
\note{
#Example
pdf("Manhattanplot.pdf")
get_manhattanplots(...) #pass arguments
dev.off()
}
plot.eigentrends <-
function(svdr, title1){
  v = svdr$v
  d = svdr$d
  ss = d^2
  Tk = signif(ss/sum(ss)* 100, 2)
  
  titles = paste("Trend ", 1:3, " (", Tk[1:3], "%)", sep = "")
  do.text = function(j) mtext(titles[j], cex=0.7, padj=-0.7, adj=1)
  range.y = range(as.numeric(v[,1:3]), na.rm=T)
  
  toplot1_1 = as.numeric(v[,1])
  toplot1_2 = as.numeric(v[,2])
  toplot1_3 = as.numeric(v[,3])
  
  plot(c(1:length(toplot1_1)), toplot1_1, type='b', ann=F, ylim=range.y)
  do.text(1)
  abline(h=0, lty=3)
  title(title1, cex.main = 1.2, font.main= 1, col.main= "purple", ylab=NULL)
  plot(c(1:length(toplot1_2)), toplot1_2, type='b', ann=F, ylim=range.y)
  do.text(2)
  abline(h=0, lty=3)
  plot(c(1:length(toplot1_3)), toplot1_3, type='b', ann=F, ylim=range.y)
  do.text(3)
  abline(h=0, lty=3)
  return(Tk)
}

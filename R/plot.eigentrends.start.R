plot.eigentrends.start <-
function(svdr, title1, pos1=1){
  # No check for valid range of pos1 is performed!!!
  v = svdr$v
  d = svdr$d
  ss = d^2
  Tk = signif(ss/sum(ss)* 100, 2)
  #  pe = signif(d/sum(d, na.rm=T)*100, 2)
  titles = paste("Trend ", pos1:(pos1+3), " (", Tk[pos1:(pos1+3)], "%)", sep = "")
  do.text = function(j) mtext(titles[j], cex=0.7, padj=-0.7, adj=1)
  range.y = range(as.numeric(v[,pos1:(pos1+3)]), na.rm=T)
  
  toplot1_1 = as.numeric(v[,pos1])
  toplot1_2 = as.numeric(v[,(pos1+1)])
  toplot1_3 = as.numeric(v[,(pos1+2)])
  
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

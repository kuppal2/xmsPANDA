get_pcascoredistplots <-
function(X=NA,Y=NA,feature_table_file,parentoutput_dir,class_labels_file,
                                sample.col.opt=c("journal", "npg", "nejm", "jco", "lancet", "custom1", "brewer.RdYlBu", "brewer.RdBu", "brewer.PuOr", 
                                                 "brewer.PRGn", "brewer.PiYG", "brewer.BrBG", "brewer.Set2", "brewer.Paired", "brewer.Dark2", "brewer.YlGnBu", "brewer.YlGn",
                                                 "brewer.YlOrRd", "brewer.YlOrBr", "brewer.PuBuGn", "brewer.PuRd", "brewer.PuBu", "brewer.OrRd", "brewer.GnBu", "brewer.BuPu",
                                                 "brewer.BuGn", "brewer.blues", "black", "grey65", "terrain", "rainbow", "heat", "topo"),
                                plots.width=2000,plots.height=2000,plots.res=300,
                                alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,pca.cex.val=3,
                                legendlocation="topright",pca.ellipse=FALSE,ellipse.conf.level=0.95,
                                filename="all",paireddesign=NA,error.bar=TRUE,lineplot.col.opt="black",
                                lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                                newdevice=FALSE,timeseries.lineplots=FALSE,alphabetical.order=FALSE,pcascale=TRUE,
                                pcacenter=TRUE,study.design="oneway",lme.modeltype="RI",cex.plots=0.8,
                                ypos.adj.factor=0.5,...)
{
  if(FALSE){
  if(length(sample.col.opt)==1){
  sample.col.opt=tolower(sample.col.opt)
  sample.col.opt=get_hexcolors_for_palettes(color.palette=sample.col.opt[1],alpha.col=alphacol[1])
  }
  
  if(length(lineplot.col.opt)==1){
    lineplot.col.opt=tolower(lineplot.col.opt)
    lineplot.col.opt=get_hexcolors_for_palettes(color.palette=lineplot.col.opt[1],alpha.col=alphacol[1])
  }
  }  
  
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  
  match_col.opt=match(sample.col.opt,c("journal","npg","nejm","jco","lancet","custom1","brewer.RdYlBu","brewer.RdBu","brewer.PuOr","brewer.PRGn","brewer.PiYG","brewer.BrBG",
                                       "brewer.Set2","brewer.Paired","brewer.Dark2","brewer.YlGnBu","brewer.YlGn","brewer.YlOrRd","brewer.YlOrBr","brewer.PuBuGn",
                                       "brewer.PuRd","brewer.PuBu",
                                       "brewer.OrRd","brewer.GnBu","brewer.BuPu","brewer.BuGn","brewer.blues","black","grey65","terrain","rainbow","heat","topo"))
  
  match_col.opt=length(which(is.na(match_col.opt)==TRUE))
    
  if(length(match_col.opt)<1){
    
    sample.col.opt=sample.col.opt[1]
  }else{
    
    if(length(grep(sample.col.opt,pattern="brewer."))>1){
      
      sample.col.opt=sample.col.opt[1]
    }
  }
  
  if(length(sample.col.opt)==1){
    sample.col.opt=tolower(sample.col.opt)
    sample.col.opt=get_hexcolors_for_palettes(color.palette=sample.col.opt[1],alpha.col=alphacol[1])
  }
  
  lineplot.col.opt=sample.col.opt
  
  res<-get_pcascoredistplots_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,
                                   class_labels_file=class_labels_file,sample.col.opt=sample.col.opt,plots.width=plots.width,
                                   plots.height=plots.height,plots.res=plots.res, alphacol=alphacol,col_vec=col_vec,
                                   pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,
                                   pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename=filename,
                                   paireddesign=paireddesign,error.bar=error.bar,lineplot.col.opt=lineplot.col.opt,
                                   lineplot.lty.option=lineplot.lty.option,newdevice=newdevice,timeseries.lineplots=timeseries.lineplots,
                                   alphabetical.order=alphabetical.order,pcascale=pcascale,pcacenter=pcacenter,study.design=study.design,lme.modeltype=lme.modeltype,cex.plots=cex.plots,ypos.adj.factor=ypos.adj.factor,...)
  
  
  return(res)
}

get_boxplots <-
function(X=NA,Y=NA,feature_table_file,parentoutput_dir,class_labels_file,
                       boxplot.col.opt=c("journal", "npg", "nejm", "jco", "lancet", "custom1", "brewer.RdYlBu", "brewer.RdBu", "brewer.PuOr", 
                                         "brewer.PRGn", "brewer.PiYG", "brewer.BrBG", "brewer.Set2", "brewer.Paired", "brewer.Dark2", "brewer.YlGnBu", "brewer.YlGn",
                                         "brewer.YlOrRd", "brewer.YlOrBr", "brewer.PuBuGn", "brewer.PuRd", "brewer.PuBu", "brewer.OrRd", "brewer.GnBu", "brewer.BuPu",
                                         "brewer.BuGn", "brewer.blues", "black", "grey65", "terrain", "rainbow", "heat", "topo"),
                       alphacol=1,newdevice=TRUE,cex.plots=0.8,replace.by.NA=FALSE,pairedanalysis=FALSE,
                       filename="",ylabel="Intensity",alphabetical.order=FALSE,name=NA,
                       add.jitter=FALSE,add.pvalues=FALSE,class.levels=NA,fill.plots=TRUE,
                       connectpairedsamples=FALSE,boxplot.type="ggplot",
                       study.design=c("multiclass","onewayanova","twowayanova","onewayanovarepeat",
                                      "twowayanovarepeat"),
                       multiple.figures.perpanel=FALSE,ggplot.type1=TRUE,replace.outliers=FALSE,
                       plot.height=8,plot.width=8,
                       extra_text=NA,group_by_mat=NA,position_dodge_width=0.75,
                       numnodes=2,hightlight.points=FALSE,ref.group.val=FALSE,facet.nrow=1,facet.ncol=NULL,...)
{
  
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  match_col.opt=match(boxplot.col.opt,c("journal","npg","nejm","jco","lancet","custom1","brewer.RdYlBu","brewer.RdBu","brewer.PuOr","brewer.PRGn","brewer.PiYG","brewer.BrBG",
                                       "brewer.Set2","brewer.Paired","brewer.Dark2","brewer.YlGnBu","brewer.YlGn","brewer.YlOrRd","brewer.YlOrBr","brewer.PuBuGn",
                                       "brewer.PuRd","brewer.PuBu",
                                       "brewer.OrRd","brewer.GnBu","brewer.BuPu","brewer.BuGn","brewer.blues","black","grey65","terrain","rainbow","heat","topo"))
  
  match_col.opt=length(which(is.na(match_col.opt)==TRUE))
  
  if(length(match_col.opt)<1){
    
    boxplot.col.opt=boxplot.col.opt[1]
  }else{
    
    if(length(grep(boxplot.col.opt,pattern="brewer."))>1){
      
      boxplot.col.opt=boxplot.col.opt[1]
    }
  }
  
  #color.palette=get_hexcolors_for_palettes(color.palette=color.palette,alpha.col=alpha.col[1])
  
  if(length(boxplot.col.opt)==1){
    boxplot.col.opt=tolower(boxplot.col.opt)
    boxplot.col.opt=get_hexcolors_for_palettes(color.palette=boxplot.col.opt[1],alpha.col=alphacol[1])
  }

  res<-suppressWarnings(suppressMessages(get_boxplots_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,boxplot.col.opt,alphacol=alphacol,newdevice=newdevice,cex.plots=cex.plots,
                     replace.by.NA=replace.by.NA,pairedanalysis=pairedanalysis,filename=filename,ylabel=ylabel,alphabetical.order=alphabetical.order,name=name,add.jitter=add.jitter,add.pvalues=add.pvalues,class.levels=class.levels,fill.plots=fill.plots,
                     connectpairedsamples=connectpairedsamples,boxplot.type=boxplot.type,study.design=study.design,
                     multiple.figures.perpanel=multiple.figures.perpanel,ggplot.type1=ggplot.type1,
                     replace.outliers=replace.outliers,plot.height=plot.height,
                     plot.width=plot.width,extra_text=extra_text,group_by_mat=group_by_mat,
                     position_dodge_width=position_dodge_width,numnodes=numnodes,
                     hightlight.points=hightlight.points,ref.group.val=ref.group.val,facet.nrow=facet.nrow,facet.ncol=facet.ncol,...)))
  
  if(newdevice==TRUE){
    
    try(dev.off(),silent=TRUE)
  }
  
  try(unlink("Rplots.pdf"),silent=TRUE)
  return(res)
}

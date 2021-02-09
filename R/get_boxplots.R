get_boxplots <-
function(X=NA,Y=NA,feature_table_file,parentoutput_dir,class_labels_file,
                       boxplot.col.opt="journal",alphacol=0.3,newdevice=TRUE,cex.plots=0.8,replace.by.NA=FALSE,pairedanalysis=FALSE,
                       filename="",ylabel="Intensity",alphabetical.order=FALSE,name=NA,
                       add.jitter=FALSE,add.pvalues=FALSE,class.levels=NA,fill.plots=TRUE,
                       connectpairedsamples=FALSE,boxplot.type="ggplot",
                       study.design=c("multiclass","onewayanova","twowayanova","onewayanovarepeat",
                                      "twowayanovarepeat"),
                       multiple.figures.perpanel=TRUE,ggplot.type1=TRUE,replace.outliers=FALSE,
                       plot.height=7,plot.width=7,
                       extra_text=NA,group_by_mat=NA,position_dodge_width=0.75,
                       numnodes=2,hightlight.points=FALSE,...)
{
  
  get_boxplots_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,boxplot.col.opt,alphacol=alphacol,newdevice=newdevice,cex.plots=cex.plots,
                     replace.by.NA=replace.by.NA,pairedanalysis=pairedanalysis,filename=filename,ylabel=ylabel,alphabetical.order=alphabetical.order,name=name,add.jitter=add.jitter,add.pvalues=add.pvalues,class.levels=class.levels,fill.plots=fill.plots,
                     connectpairedsamples=connectpairedsamples,boxplot.type=boxplot.type,study.design=study.design,
                     multiple.figures.perpanel=multiple.figures.perpanel,ggplot.type1=ggplot.type1,
                     replace.outliers=replace.outliers,plot.height=plot.height,
                     plot.width=plot.width,extra_text=extra_text,group_by_mat=group_by_mat,
                     position_dodge_width=position_dodge_width,numnodes=numnodes,hightlight.points=hightlight.points,...)
  
  if(newdevice==TRUE){
    
    try(dev.off(filename),silent=TRUE)
  }
}

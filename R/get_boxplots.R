get_boxplots <-
function(X=NA,Y=NA,feature_table_file,parentoutput_dir,class_labels_file,boxplot.col.opt="white",sample.col.opt="rainbow",alphacol=0.3,newdevice=TRUE,cex.plots=0.8,replace.by.NA=FALSE,pairedanalysis=FALSE,
filename="",ylabel="Intensity",alphabetical.order=FALSE,name=NA,add.jitter=FALSE,add.pvalues=FALSE,class.levels=NA,fill.plots=FALSE,connectpairedsamples=FALSE,boxplot.type="ggplot",study.design=c("multiclass","onewayanova","twowayanova","onewayanovarepeat","twowayanovarepeat"))
{

get_boxplots_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,boxplot.col.opt,sample.col.opt=sample.col.opt,alphacol=alphacol,newdevice=newdevice,cex.plots=cex.plots,
replace.by.NA=replace.by.NA,pairedanalysis=pairedanalysis,filename=filename,ylabel=ylabel,alphabetical.order=alphabetical.order,name=name,add.jitter=add.jitter,add.pvalues=add.pvalues,class.levels=class.levels,fill.plots=fill.plots,
connectpairedsamples=connectpairedsamples,boxplot.type=boxplot.type,study.design=study.design)

}

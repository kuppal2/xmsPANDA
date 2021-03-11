get_pcascoredistplots <-
function(X=NA,Y=NA,feature_table_file,parentoutput_dir,class_labels_file,
                                sample.col.opt="journal",plots.width=2000,plots.height=2000,plots.res=300,
                                alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,pca.cex.val=3,
                                legendlocation="topright",pca.ellipse=FALSE,ellipse.conf.level=0.95,
                                filename="all",paireddesign=NA,error.bar=TRUE,lineplot.col.opt="black",
                                lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                                newdevice=FALSE,timeseries.lineplots=FALSE,alphabetical.order=FALSE,pcascale=TRUE,
                                pcacenter=TRUE,study.design="oneway",lme.modeltype="RI",cex.plots=0.8,
                                ypos.adj.factor=0.5,...)
{
  
  
  res<-get_pcascoredistplots_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,sample.col.opt=sample.col.opt,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, alphacol=alphacol,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename=filename,paireddesign=paireddesign,error.bar=error.bar,lineplot.col.opt=lineplot.col.opt,lineplot.lty.option=lineplot.lty.option,newdevice=newdevice,timeseries.lineplots=timeseries.lineplots,alphabetical.order=alphabetical.order,pcascale=pcascale,pcacenter=pcacenter,study.design=study.design,lme.modeltype=lme.modeltype,cex.plots=cex.plots,ypos.adj.factor=ypos.adj.factor,...)
  
  
  return(res)
}

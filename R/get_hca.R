get_hca <-
function(feature_table_file=NA,parentoutput_dir,class_labels_file=NA,X=NA,Y=NA,heatmap.col.opt="RdBu",cor.method="spearman",
                  is.data.znorm=FALSE,analysismode="classification",
                  sample.col.opt=c("journal", "npg", "nejm", "jco", "lancet", "custom1", "brewer.RdYlBu", "brewer.RdBu", "brewer.PuOr", 
                                   "brewer.PRGn", "brewer.PiYG", "brewer.BrBG", "brewer.Set2", "brewer.Paired", "brewer.Dark2", "brewer.YlGnBu", "brewer.YlGn",
                                   "brewer.YlOrRd", "brewer.YlOrBr", "brewer.PuBuGn", "brewer.PuRd", "brewer.PuBu", "brewer.OrRd", "brewer.GnBu", "brewer.BuPu",
                                   "brewer.BuGn", "brewer.blues", "black", "grey65", "terrain", "rainbow", "heat", "topo"),
                  plots.width=8,plots.height=8,plots.res=600, plots.type="cairo", alphacol=1, hca_type="two-way",
                  newdevice=FALSE,input.type="intensity",mainlab="",cexRow=0.5, cexCol=0.5,plot.bycluster=FALSE,color.rows=TRUE,
                  similarity.matrix="correlation",deepsplit=4,minclustsize=2,mergeCutHeight=0.05,num_nodes=2,alphabetical.order=FALSE,
                  pairedanalysis=FALSE,cutree.method=c("dynamic","default"),study.design=c("multiclass","onewayanova","twowayanova","onewayanovarepeat",
                                                                                 "twowayanovarepeat"),labRow.value = FALSE,
                  labCol.value = FALSE,power_val=6,row.col.opt="journal",show.silhouette=FALSE,cexLegend=0.7,ylab_text="",xlab_text="")
{
  match_col.opt=match(sample.col.opt,c("journal","npg","nejm","jco","lancet","custom1","brewer.RdYlBu","brewer.RdBu","brewer.PuOr","brewer.PRGn","brewer.PiYG","brewer.BrBG",
                                       "brewer.Set2","brewer.Paired","brewer.Dark2","brewer.YlGnBu","brewer.YlGn","brewer.YlOrRd","brewer.YlOrBr","brewer.PuBuGn",
                                       "brewer.PuRd","brewer.PuBu",
                                       "brewer.OrRd","brewer.GnBu","brewer.BuPu","brewer.BuGn","brewer.blues","black","grey65","terrain","rainbow","heat","topo"))
  
  match_col.opt=length(which(is.na(match_col.opt)==TRUE))
  
  cutree.method=cutree.method[1]
  
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
  
  
  h73<-get_hca_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=is.data.znorm,analysismode=analysismode,
                     sample.col.opt=sample.col.opt,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, plots.type=plots.type, 
                     alphacol=alphacol, hca_type=hca_type,newdevice=newdevice,input.type=input.type,mainlab=mainlab,cexRow=cexRow, 
                     cexCol=cexCol,plot.bycluster=plot.bycluster,color.rows=color.rows,similarity.matrix,deepsplit,minclustsize,
                     mergeCutHeight,num_nodes,alphabetical.order,pairedanalysis,cutree.method,study.design,labRow.value,labCol.value,power_val,row.col.opt,
                     show.silhouette,cexLegend,ylab_text,xlab_text)
  return(h73)
}

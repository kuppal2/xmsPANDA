plot_graph <-
function(df,net_node_colors=c("orange", "green","blue","pink"),graphmethod="radial",label.cex=0.3,filename="",
                     net_node_shape=c("rectangle","circle","triangle","star","square","csquare","crectangle","vrectangle"),maxnodesperclass=100,
                     seednum=100,vertex.size=6,interactive=TRUE,Xname=NA,Yname=NA,Zname=NA,Wname=NA,classname=NA,layout.type="fr1",mtext.val=NA){
  
  df<-as.data.frame(df)
  
 
  df$from<-as.character(df$from)
  df$to<-as.character(df$to)
  
  df$weight<-as.numeric(df$weight)
 #save(df,file="df.Rda")
  
  if(length(which(df$weight==0))>0){
    
    df<-df[-which(df$weight==0),]
  }
  
  if(nrow(df)>0){
    
  df<-df[order(df$weight,df$to),] #decreasing=FALSE),]
  
  nodes_vec_class<-gsub(df$to,pattern="[0-9|_|\\.]*",replacement="") #c("X","Y","Z","W") #
  
  df$from<-gsub(df$from,pattern="X",replacement="")
  df$to<-gsub(df$to,pattern="Y",replacement="")
  
  
  dfall<-df
  
  
  t1<-c("X","Y","Z","W") #c(Xname,Yname,Zname,Wname) #levels(as.factor(nodes_vec_class))

  
  
  
  if(graphmethod=="BAM"){
    sg<-graphBAM(df)
    nodes_vec<-nodes(sg)
  }else{
    sg <- graph.data.frame(d = df, directed = FALSE)
    #sg <- graph.adjacency(df, mode="undirected", weighted=TRUE)
    nodes_vec<-V(sg)$name
  }
  nodes_vec_class<-gsub(nodes_vec,pattern="[0-9]*",replacement="")
  nodes_col_vec<-rep("#11BCFF99",length(nodes_vec))
  nodes_shape_vec<-rep("circle",length(nodes_vec))
  nodes_size_vec<-rep(vertex.size,length(nodes_vec))
  
  if(graphmethod=="BAM"){
    shapes_sel_vec<-c("circle","ellipse","box","rectangle")
  }else{
    
    # clips as a circle
    add_shape("triangle",plot=mytriangle) #, clip=shapes("circle")$clip,plot=mytriangle)
    # no clipping, edges will be below the vertices anyway
    add_shape("star", clip=shape_noclip,plot=mystar, parameters=list(vertex.norays=5))
    
    shapes_sel_vec<-net_node_shape  #c("triangle","sphere","rectangle","star")
  }
  
  for(t1_ind in 1:length(t1)){
    
    nodes_col_vec[which(nodes_vec_class==t1[t1_ind])]<-net_node_colors[t1_ind]
    
    nodes_shape_vec[which(nodes_vec_class==t1[t1_ind])]<-shapes_sel_vec[t1_ind]
    
    #nodes_size_vec[which(nodes_vec_class==t1[t1_ind])]<-3/t1_ind
    
  }
  
  nAttrs<-list()
  nAttrs$fillcolor<-nodes_col_vec
  names(nAttrs$fillcolor)<-nodes_vec #(sg)
  nAttrs$shape<-nodes_shape_vec
  names(nAttrs$shape)<-nodes_vec #(sg)
  #nAttrs$size<-nodes_size_vec
  #names(nAttrs$size)<-nodes_vec #(sg)
  # clips as a circle
  add_shape("triangle",plot=mytriangle) #, clip=shapes("circle")$clip,plot=mytriangle)
  # no clipping, edges will be below the vertices anyway
  add_shape("star", clip=shape_noclip,plot=mystar, parameters=list(vertex.norays=5))
  
  
  if(graphmethod=="radial"){
    edge_colors<-rep("blue",length(E(sg)$weight))
    edge_colors[which(E(sg)$weight>=0)]<-"red"
    
    V(sg)$color<-nodes_col_vec
    V(sg)$shape<-nodes_shape_vec
    V(sg)$label.cex<-label.cex
    
    V(sg)$vertex.size<-vertex.size
    E(sg)$color=edge_colors
    
    set.seed(seednum)
    
    if(FALSE){
      if(is.na(maxnodesperclass)==TRUE){
        cytoscape_fname<-paste(filename,"cytoscapeall.gml",sep="")
        rda_fname<-paste(filename,"all.Rda",sep="")
      }else{
        cytoscape_fname<-paste(filename,"cytoscapetop",maxnodesperclass,".gml",sep="")
        rda_fname<-paste(filename,"top",maxnodesperclass,".Rda",sep="")
      }
      fname1<-paste(filename,"_linkmatrix.Rda",sep="")
      #  ##save(df,file=fname1)
      ###save(sg,file=rda_fname)
      write.graph(sg, file =cytoscape_fname, format = "gml")
    }
    
    fname1<-paste("Tables/",filename,"_linkmatrix.txt",sep="")
    write.table(df,file=fname1,sep="\t",row.names=FALSE)
    
    {
      
      set.seed(seednum)
      
      l = get_layout(sg,layout.type=layout.type) #layout.fruchterman.reingold(sg, weights = (1 -abs(E(sg)$weight)))
      
      
      set.seed(seednum)
      
      check_plot<-plot.igraph(sg,layout=l,vertex.size=vertex.size,vertex.label=V(sg)$name,vertex.label.cex=label.cex,edge.color=edge_colors)
      
      
      if(is(check_plot,"try-error")){
        
        set.seed(seednum)
        sg$layout <- l
        set.seed(seednum)
        
        plot.igraph(sg,layout=l,vertex.size=vertex.size,vertex.label=V(sg)$name,vertex.label.cex=label.cex,edge.color=edge_colors) #,silent=TRUE)
        
        
      }
      #mtext("Red: +ve correlation; Blue: -ve correlation",side=1)
      
      if(is.na(classname)==TRUE){
        
        if(is.na(maxnodesperclass)==TRUE){
          
          if(is.na(mtext.val)==FALSE){
            mtext(mtext.val,side=3,line=3,cex=1,adj=NA)
          }else{
            mtext("Showing all pairwise correlations",side=3,line=3,cex=1,adj=NA)
          }
        }else{
          
          if(is.na(mtext.val)==FALSE){
            mtext(mtext.val,side=3,line=3,cex=1,adj=NA)
          }else{
          mtext(paste("Showing top ",maxnodesperclass," correlations",sep=""),side=3,line=3,cex=1,adj=NA)
          }
        }
        #mtext("Using all samples",side=3,line=3,cex=0.6,adj=NA)
      }else{
        
        mtext(paste("Using samples in class ",classname,sep=""),side=3,line=3,cex=0.6,adj=NA)
      }
      
      mtext("(Edges) Red: +ve correlation; Blue: -ve correlation",line=0,side=1,cex=0.8,adj=0)
      
      mtext_community<-paste("(Nodes) ",net_node_shape[1],": ",Xname,"; ",net_node_shape[2],": ",Yname,sep="")
      if(is.na(Zname)==FALSE){
        mtext_community<-paste(mtext_community,"; ",net_node_shape[3],": ",Zname,sep="")
        
      }
      if(is.na(Wname)==FALSE){
        mtext_community<-paste(mtext_community,"; ",net_node_shape[4],": ",Wname,sep="")
        
      }
      
      mtext(mtext_community,side=1,cex=0.8,line=1,adj=0)
      
      
      filename_text<-paste(filename,".png",sep="")
      try(mtext(filename_text,line=3,cex=0.6,col="brown",side=1,adj=0),silent=TRUE)
      
    }
    ###save(list=ls(),file="plot.Rda")
    nodes<-V(sg)$name
    nodes<-as.data.frame(nodes)
    links<-df[,c(1:2)]
    links<-as.data.frame(links)
    
    if(is.na(maxnodesperclass)==FALSE){
      if(interactive==TRUE){
        
        tkid <- tkplot(sg,canvas.width=1050, canvas.height=750)
        
        
        colnames(links)<-c("from","to")
        colnames(nodes)<-c("id")
        
        nodes$color.background<-"green"
        nodes_vec_class<-gsub(nodes[,1],pattern="[0-9]+",replacement="")
        nodes$shape<-"dot"
        nodes$type=nodes_vec_class
        nodes_vec_class<-unique(nodes_vec_class)
        
        for(i in 1:length(unique(nodes_vec_class))){
          
          nodes$shape[which(nodes$type==nodes_vec_class[i])] <- net_node_shape[i]
          nodes$color.background[which(nodes$type==nodes_vec_class[i])] <- net_node_colors[i]
        }
        
        nodes$shadow <- TRUE # Nodes will drop shadow
        nodes$title <- nodes$media # Text on click
        nodes$label <- nodes$type.label # Node label
        nodes$size <- vertex.size #nodes$audience.size # Node size
        nodes$borderWidth <- 2 # Node border width
        
        #nodes$color.background <- c("slategrey", "tomato", "gold")[nodes$type]
        nodes$color.border <- "black"
        #nodes$color.highlight.background <- "orange"
        #nodes$color.highlight.border <- "darkred"
        links$color <- edge_colors
        
        #network<-visNetwork(nodes, links)
        
        #vis##save(network, file = "network.html")
        
        
      }
    }
  }else{
    set.seed(seednum)
    plot(sg,nodeAttrs=nAttrs)
  }
  return(list(nodes=nodes,links=links,sg=sg,link_matrix=df,layout_matrix=l))
  }else{
    
    return(NULL)
  }
}

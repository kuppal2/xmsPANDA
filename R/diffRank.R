diffRank <-
function(adjMtrxSample1,adjMtrxSample2,degree.centrality.method="eigenvector",node_names=NA,class_labels=c(1,2),
                   plot_graph_bool=FALSE,label.cex=0.5,vertex.size=0.7,ignore.edge.weights=FALSE)
  {
  
  suppressMessages(library(igraph))
  
 #save(adjMtrxSample1,adjMtrxSample2,degree.centrality.method,node_names,ignore.edge.weights,class_labels,file="adj.Rda")
  #N is the number if genes in the sample.
  N=dim(adjMtrxSample1)[1];
  
  
  diag(adjMtrxSample1)<-0
  diag(adjMtrxSample2)<-0
  
  if(is.na(node_names)==FALSE){
    
    rownames(adjMtrxSample1)<-node_names
    colnames(adjMtrxSample1)<-node_names
    rownames(adjMtrxSample2)<-node_names
    colnames(adjMtrxSample2)<-node_names
    
  }
  
  if(ignore.edge.weights==FALSE){
  # Create a graph from the dataset; one graph for each condition.
  graph1=graph.adjacency(adjMtrxSample1, mode="undirected",weighted=TRUE)#, diag=FALSE )
  graph2=graph.adjacency(adjMtrxSample2, mode="undirected",weighted=TRUE)#, diag=FALSE )
  
  }else{
    graph1=graph.adjacency(adjMtrxSample1, mode="undirected")#, diag=FALSE )
    graph2=graph.adjacency(adjMtrxSample2, mode="undirected")#, diag=FALSE )
    
  }
  
  if(plot_graph_bool==TRUE){
    
    edge_colors<-rep("blue",length(E(graph1)$weight))
    edge_colors[which(E(graph1)$weight>=0)]<-"red"
   # V(sg)$color<-nodes_col_vec
    #V(sg)$shape<-nodes_shape_vec
    V(graph1)$label.cex<-label.cex
    
    V(graph1)$vertex.size<-vertex.size
    E(graph1)$color=edge_colors
    
    #E(graph1)$weight <- edge.betweenness(graph1)
    
    if(ignore.edge.weights==FALSE){
     
      
    c1 = cluster_fast_greedy(graph1,weights=abs(E(graph1)$weight)) #cluster_louvain(graph1,weights = NA)
    c2 = cluster_fast_greedy(graph2,weights=abs(E(graph2)$weight)) #cluster_louvain(graph2,weights=NA)
    }else{
      c1 = cluster_fast_greedy(graph1,weights=NA) #cluster_louvain(graph1,weights = NA)
      c2 = cluster_fast_greedy(graph2,weights=NA) 
      
    }
    #clustering coefficient
    cc1<-transitivity(graph1)
    
    cc2<-transitivity(graph2)
    plot.layout.type="circular"
    
    if(plot.layout.type=="circular"){
    #layout circular
    coords1<-layout_in_circle(graph1,order=V(graph1))
    coords2<-layout_in_circle(graph2,order=V(graph2))
    
      
    main_graph1=paste(class_labels[1],"\nclustering coefficient: ",round(cc1,3),sep="")
    main_graph2=paste(class_labels[2],"\nclustering coefficient: ",round(cc2,3),sep="")
    
    set.seed(555)
    plot.igraph(graph1,main=main_graph1,vertex.label=V(graph1)$name,vertex.label.cex=label.cex,edge.color=edge_colors,
                vertex.size=vertex.size,layout=coords1,vertex.color=membership(c1))
    
    edge_colors<-rep("blue",length(E(graph2)$weight))
    edge_colors[which(E(graph2)$weight>=0)]<-"red"
    
    # V(sg)$color<-nodes_col_vec
    #V(sg)$shape<-nodes_shape_vec
    V(graph2)$label.cex<-label.cex
    V(graph2)$vertex.size<-vertex.size
      E(graph2)$color=edge_colors
      set.seed(555)
      plot.igraph(graph2,main=main_graph2,vertex.label=V(graph2)$name,vertex.label.cex=label.cex,edge.color=edge_colors,
                  vertex.size=vertex.size,layout=coords2,vertex.color=membership(c2))
    }else{
      state_col = c("TssA" = "#E41A1C",    "TssAFlnk" = "#E41A1C",
                    "TxFlnk" = "#E41A1C",  "Tx" = "#E41A1C",
                    "TxWk" = "#E41A1C",    "EnhG" = "#E41A1C",
                    "Enh" = "#E41A1C",     "ZNF/Rpts" = "#E41A1C",
                    "Het" = "#377EB8",     "TssBiv" = "#377EB8",
                    "BivFlnk" = "#377EB8", "EnhBiv" = "#377EB8",
                    "ReprPC" = "#377EB8",  "ReprPCWk" = "#377EB8",
                    "Quies" = "black")
      
      # one for rows and one for columns
      state_col2 = c(state_col, state_col)
      names(state_col2) = c(rownames(mat), colnames(mat))
      
      colmat = rep(state_col2[rownames(mat)], n_states)
      colmat = rgb(t(col2rgb(colmat)), maxColorValue = 255)
      
      #qati = quantile(mat, 0.7)
      colmat[mat > 0.7] = paste0(colmat[mat > 0.7], "A0")
      colmat[mat <= (-0.7)] = paste0(colmat[mat <= (-0.7)], "20")
      dim(colmat) = dim(mat)
      gm1<-as_edgelist(graph1, names = FALSE)
      gm1<-cbind(gm1,as.matrix(E(graph1)$weight))
      colnames(gm1)<-c("X","Y","weight")
      gm1<-as.data.frame(gm1)
      set.seed(555)
      chordDiagram(gm1,annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))
      
      gm2<-as_edgelist(graph2, names = FALSE)
      gm2<-cbind(gm2,as.matrix(E(graph2)$weight))
      colnames(gm2)<-c("X","Y","weight")
      gm2<-as.data.frame(gm2)
      set.seed(555)
      chordDiagram(gm2,annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))
      
    } 
      set.seed(555)
      try(plot_dendrogram(c1,cex=label.cex),silent=TRUE)
      
      set.seed(555)
      try(plot_dendrogram(c2,cex=label.cex),silent=TRUE)
      
      #save(graph1,graph2,c1,c2,file="graphdebug.Rda")
      
  }
  
  if(degree.centrality.method=="betweenness"){
    
    betweennessGraph1=betweenness(graph1,directed = FALSE,weights=abs(E(graph1)$weight),normalized=TRUE);
    betweennessGraph2=betweenness(graph2,directed = FALSE,weights=abs(E(graph2)$weight),normalized=TRUE);
    solution=abs(betweennessGraph1-betweennessGraph2)
    
  }else{
    
    if(degree.centrality.method=="eigenvector"){
          #eigenvector
          eigenGraph1=eigen_centrality(graph1,directed = FALSE,weights=abs(E(graph1)$weight))$vector #,normalized=TRUE);
          eigenGraph2=eigen_centrality(graph2,directed = FALSE,weights=abs(E(graph2)$weight))$vector #,normalized=TRUE);
          solution=abs(eigenGraph1-eigenGraph2)
    }else{
        if(degree.centrality.method=="hybrid.DBC" | degree.centrality.method=="hybrid.DEC"){
          
          if(degree.centrality.method=="hybrid.DBC"){
            betweennessGraph1=betweenness(graph1,directed = FALSE,weights=abs(E(graph1)$weight),normalized=TRUE);
            betweennessGraph2=betweenness(graph2,directed = FALSE,weights=abs(E(graph2)$weight),normalized=TRUE);
            
            
            
            DBC=abs(betweennessGraph1-betweennessGraph2)
            
          }else{
            eigenGraph1=eigen_centrality(graph1,directed = FALSE,weights=abs(E(graph1)$weight))$vector #,normalized=TRUE);
            eigenGraph2=eigen_centrality(graph2,directed = FALSE,weights=abs(E(graph2)$weight))$vector #,normalized=TRUE);
            
            DBC=abs(eigenGraph1-eigenGraph2)
            
          }
          #Create and initialize Delta_C_i
          Delta_C_i=abs(adjMtrxSample1-adjMtrxSample2);
          
          
          Delta_C_i<-apply(Delta_C_i,1,function(x)
          {
            
            if(sum(x,na.rm=TRUE)==0)
            {
              x=rep(1/N,length(x))
            }
            x=x/sum(x)
            return(x)
            #})
          })
          
          #Delta_C_i<-do.call(rbind,Delta_C_i)
          Delta_C_i<-t(Delta_C_i)
          error=100;
          count=1;
          eps=0.001
          lambda=0.5
          
          
          # Create and initialize solution array to 1/N.
          # solution contains the ranks for all genes.
          solution=array(1/N,dim=c(N,1))
          solutionEachIteration=solution;
          
          # eps (EPSILON) is a value of the difference between 2 consecutive solutions to stop the iterations.
          # Do iterate while the difference between 2 consecutive solutions is not that much big (> eps).
          while(error > eps)
          {
            count=count+1;
            
            #Keep the previous solution
            formerSoulution=solution;
            
            #Find a new solution.
            for (v in 1:N){
              
              #solution<-lapply(1:N,function(v){
              s = sum(Delta_C_i[,v] * solution); # Solution, is the (Pi ) in this case
              solution[v]=(1-lambda) * DBC[v] + lambda * s;
              #solution[v]=(lambda) * DBC[v] + lambda * s + lambda*DEC[v];
              #res=(1-lambda) * DBC[v] + lambda * s;
              #return(res)
              #})
            }
            solution<-unlist(solution)
            solutionEachIteration=array(c(solutionEachIteration,solution),dim=c(N,count))
            
            #find the error to stop the iteration, the error in this case that there no significant difference between the old and the new solution
            error=sum((formerSoulution - solution)^2)
            
            
          }
        }
      
    }
  }

  
  
  
  
  #return(list(solution=solution,DBC=DBC,DeltaK=DeltaK,DeltaK1=Delta_C_i))
  ##save(solution,DBC,DeltaK,Delta_C_i,DEC,adjMtrxSample1,file="diffrank.Rda")
  return(solution)
}

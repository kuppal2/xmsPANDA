diffRank <-
function(adjMtrxSample1,adjMtrxSample2,degree.centrality.method="eigenvector"){
  
  suppressMessages(library(igraph))
  
  #N is the number if genes in the sample.
  N=dim(adjMtrxSample1)[1];
  
  # Create a graph from the dataset; one graph for each condition.
  graph1=graph.adjacency(adjMtrxSample1, mode="undirected",weighted=TRUE)#, diag=FALSE )
  graph2=graph.adjacency(adjMtrxSample2, mode="undirected",weighted=TRUE)#, diag=FALSE )
  
  if(degree.centrality.method=="betweenness"){
    
    betweennessGraph1=betweenness(graph1,directed = FALSE,weights=abs(E(graph1)$weight),normalized=TRUE);
    betweennessGraph2=betweenness(graph2,directed = FALSE,weights=abs(E(graph2)$weight),normalized=TRUE);
    
    
    
    DBC=abs(betweennessGraph1-betweennessGraph2)
    
  }else{
    #eigenvector
    eigenGraph1=eigen_centrality(graph1,directed = FALSE,weights=abs(E(graph1)$weight),normalized=TRUE)$vector #,normalized=TRUE);
    eigenGraph2=eigen_centrality(graph2,directed = FALSE,weights=abs(E(graph2)$weight),normalized=TRUE)$vector #,normalized=TRUE);
    
    DBC=abs(eigenGraph1-eigenGraph2)
  }
  
  
  
  
  #Create and initialize Delta_C_i
  Delta_C_i=abs(adjMtrxSample1-adjMtrxSample2);
  
  
  Delta_C_i<-apply(Delta_C_i,1,function(x)
  {
    
    if(sum(x)==0)
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
    #for (v in 1:N){
    
    solution<-lapply(1:N,function(v){
      s = sum(Delta_C_i[,v] * solution); # Solution, is the (Pi ) in this case
      solution[v]=(1-lambda) * DBC[v] + lambda * s;
      #solution[v]=(lambda) * DBC[v] + lambda * s + lambda*DEC[v];
      #res=(1-lambda) * DBC[v] + lambda * s;
      #return(res)
    })
    # }
    solution<-unlist(solution)
    solutionEachIteration=array(c(solutionEachIteration,solution),dim=c(N,count))
    
    #find the error to stop the iteration, the error in this case that there no significant difference between the old and the new solution
    error=sum((formerSoulution - solution)^2)
    
    #print(error)
  }
  #print(head(solution))
  
  #rnk = rank(1-solution)
  DeltaK=apply(Delta_C_i,1,mean)
  
  #return(list(solution=solution,DBC=DBC,DeltaK=DeltaK,DeltaK1=Delta_C_i))
  ##save(solution,DBC,DeltaK,Delta_C_i,DEC,adjMtrxSample1,file="diffrank.Rda")
  return(solution)
}

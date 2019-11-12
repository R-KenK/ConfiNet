#' Convert an adjacency matrix to an edge list with repeated lines
#'
#' Wrapper for igraph functions in order to obtain an adjacency matrix where the cell a pair represent the edge weight. Useful to obtain the desired association index after bootstrapping the weighted edge list.
#'
#' @param edge data.table or two-columns matrix. Each line represent an edge (order matters if the graph is directed), repetition of a given (ordered) pair represent the weight of an edge.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return adjacency matrix
#' @export
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph E
#' @importFrom igraph simplify
#' @importFrom igraph as.undirected
#' @importFrom igraph as_adjacency_matrix
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[upper.tri(Adj)]<- sample(0:2,n*(n-1)/2,replace = TRUE)
#' Adj[lower.tri(Adj)]<- sample(0:1,n*(n-1)/2,replace = TRUE)
#' Adj
#'
#' edge<- KenNet::adj.to.edge(Adj,mode = "directed")
#' edge.to.adj(edge,mode="directed")
#' edge.up<- KenNet::adj.to.edge(Adj,mode = "upper")
#' edge.to.adj(edge.up,mode="upper")
#'
#' edge.plus<- KenNet::adj.to.edge(Adj,mode = "plus")
#' edge.to.adj(edge.plus,mode="plus")
#'
#' edge.low<- KenNet::adj.to.edge(Adj,mode = "lower")
#' edge.to.adj(edge.low,mode="lower")
edge.to.adj<- function(edge,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus")){
  if(length(mode)>1) {mode<- "directed"}
  if(mode %in% c("directed", "undirected", "max","min", "upper", "lower", "plus")){
    if(mode=="directed"){
      directed<- TRUE;
      type<- "both";
    }else{
      directed<- FALSE;
      switch (mode,"upper" =,"lower" =type<- mode,type<- "both")
    }
  }else{
    stop("Argument mode not recognized. Should be \"directed\", \"undirected\", \"max\",\"min\", \"upper\", \"lower\", \"plus\"")
  }

  graph<- igraph::graph_from_edgelist(as.matrix(edge),directed = directed);igraph::E(graph)$weight<-1;

  if(directed){
    graph<- igraph::simplify(graph,remove.multiple = TRUE,remove.loops = FALSE,edge.attr.comb = "sum")
  }else{
    graph<- igraph::as.undirected(igraph::simplify(graph,remove.multiple = TRUE,remove.loops = FALSE,edge.attr.comb = "sum"),mode = "collapse")
  }


  Adj<- as.matrix(igraph::as_adjacency_matrix(graph,attr = "weight",type=type))
  Adj<- Adj[order(row.names(Adj)),order(colnames(Adj))]

  switch (mode,
          "upper" ={Adj<- Adj+t(Adj); Adj[lower.tri(Adj)]<-0;Adj},
          "lower" ={Adj<- Adj+t(Adj); Adj[upper.tri(Adj)]<-0;Adj},
          Adj)
}

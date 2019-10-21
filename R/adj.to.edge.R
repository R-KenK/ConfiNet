#' Convert an adjacency matrix to an edge list with repeated lines
#'
#' Wrapper for igraph functions in order to obtain an edge list where the number of line for a pair represent the edge weight. Useful for bootstrapping dyads when the observation protocol is based on dyadic observations.
#'
#' @param Adj adjacency matrix. Can be directed or undirected. Should probably be weighted if the network is to be bootstrapped...
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return edge list (as a data table) with column id and tar for the two nodes involved in a dyad, with repeated lines representing the weight of the edge.
#' @export
#' @importFrom igraph graph.adjacency
#' @importFrom igraph get.edgelist
#' @importFrom data.table data.table
#'
#' @examples
#' set.seed(42)
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[upper.tri(Adj)]<- sample(0:3,n*(n-1)/2,replace = TRUE)
#' Adj[lower.tri(Adj)]<- sample(0:3,n*(n-1)/2,replace = TRUE)
#' Adj
#'
#' adj.to.edge(Adj,mode = "directed")
#' adj.to.edge(Adj,mode = "plus")

adj.to.edge<- function(Adj,mode){
  graph<- igraph::graph_from_adjacency_matrix(Adj,mode = mode)
  edgelist<- data.table::data.table(igraph::get.edgelist(graph))
  colnames(edgelist)<-c("id","tar")
  edgelist
}

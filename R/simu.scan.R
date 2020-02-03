#' Simulate a group scan
#'
#' Simulate a list of instantaneous adjacency matrices
#'
#' @param n_scans integer, number of scans to simulate.
#' @param n_nodes integer, number of nodes to simulate.
#' @param nodes_names optional character vector of the nodes' names, otherwise simply numbered.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return list of adjacency matrices
#' @export
#'
#' @examples
#' set.seed(42)
#' n<- 10L;nodes<- letters[1:n];n_scans<- 50L;
#'
#' simu.scan(n_scans,n_nodes = n,mode = "undirected")
#' simu.scan(n_scans,n_nodes = n,mode = "directed")
#' simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "plus")
#'
simu.scan<- function(n_scans,n_nodes, nodes_names=as.character(1:n_nodes),mode = c("directed", "undirected", "max","min", "upper", "lower", "plus")){
  if(length(mode)>1) {mode<- "directed"}
  switch(mode,
         "directed" ={
           lapply(1:n_scans,
                  function(scan) {
                    Scan<- matrix(data = 0,nrow = n_nodes,ncol = n_nodes,dimnames = list(nodes_names,nodes_names))
                    Scan[upper.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    Scan[lower.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    Scan
                  }
           )
         },
         "undirected" ={
           lapply(1:n_scans,
                  function(scan) {
                    Scan<- matrix(data = 0,nrow = n_nodes,ncol = n_nodes,dimnames = list(nodes_names,nodes_names))
                    Scan[upper.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    Scan+t(Scan)
                  }
           )
         },
         "max" ={
           lapply(1:n_scans,
                  function(scan) {
                    Scan<- matrix(data = 0,nrow = n_nodes,ncol = n_nodes,dimnames = list(nodes_names,nodes_names))
                    Scan[upper.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    Scan[lower.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    ifelse(Scan+t(Scan)>=1,1,0) #conserve a connection between nodes if there's one in either directions (either adjacency triangle)
                  }
           )
         },
         "min" ={
           lapply(1:n_scans,
                  function(scan) {
                    Scan<- matrix(data = 0,nrow = n_nodes,ncol = n_nodes,dimnames = list(nodes_names,nodes_names))
                    Scan[upper.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    Scan[lower.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    ifelse(Scan+t(Scan)==2,1,0) #only conserve a connection between nodes who have one in both directions (each adjacency triangle)
                  }
           )
         },
         "upper" ={
           lapply(1:n_scans,
                  function(scan) {
                    Scan<- matrix(data = 0,nrow = n_nodes,ncol = n_nodes,dimnames = list(nodes_names,nodes_names))
                    Scan[upper.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    Scan
                  }
           )
         },
         "lower" ={
           lapply(1:n_scans,
                  function(scan) {
                    Scan<- matrix(data = 0,nrow = n_nodes,ncol = n_nodes,dimnames = list(nodes_names,nodes_names))
                    Scan[lower.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    Scan
                  }
           )
         },
         "plus" = {
           lapply(1:n_scans,
                  function(scan) {
                    Scan<- matrix(data = 0,nrow = n_nodes,ncol = n_nodes,dimnames = list(nodes_names,nodes_names))
                    Scan[upper.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    Scan[lower.tri(Scan)]<- sample(0:1,n_nodes*(n_nodes-1)/2,replace = TRUE)
                    Scan<- Scan+t(Scan)
                    Scan
                  }
           )
         }
  )
}

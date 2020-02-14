#' Single group scan
#' Generate a random group scan based on provided reference adjacency matrix and sampling effort.
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return a square binary matrix representing the whle group scan
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:30,n*(n-1),replace = TRUE)
#' Adj
#'
#' do.scan(Adj,42)
#' Reduce("+",lapply(1:42,function(s) do.scan(Adj,42)))

do.scan<- function(Adj,total_scan,
                   mode = c("directed", "undirected", "max","min", "upper", "lower", "plus")){
  if(nrow(Adj)==ncol(Adj)) {n<- nrow(Adj);nodes_names<- row.names(Adj)} else {stop("Adj is not a square matrix")}
  mode<- match.arg(mode)
  Scan<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes_names,nodes_names))
  Adj.subfun<- switch(mode,
                      "directed" = ,
                      "undirected" = ,
                      "max" = ,
                      "min" = ,
                      "plus" = non.diagonal,
                      "upper" = upper.tri,
                      "lower" =  lower.tri
  )

  prob<- scale.to.binary.prob(Adj=Adj,total_scan=total_scan,mode = mode)
  Scan[Adj.subfun(Scan)]<- sapply(1:length(Scan[Adj.subfun(Scan)]),
                                  function(dyad) {
                                    Scan[Adj.subfun(Scan)][dyad]<- sample(c(1,0),1,replace = TRUE,prob=prob[dyad,c("present","absent")])
                                  }
  )

  switch(mode,
         "undirected" = ,
         "max" = ifelse(Scan+t(Scan)>=1,1,0), #conserve a connection between nodes if there's one in either directions (either adjacency triangle)
         "min" = ifelse(Scan+t(Scan)==2,1,0), #only conserve a connection between nodes who have one in both directions (each adjacency triangle)
         "plus" = Scan+t(Scan),
         "directed" = ,
         "upper" = ,
         "lower" =  Scan
  )
}

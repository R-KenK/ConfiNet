#' Calculate network metric
#'
#' Wrapper for internal use.
#'
#' @param graph graph (igraph) or list of graphsp for which the metric will be calculated
#' @param FUN function to calculate the desired metric. Should take at least
#' @param args.list arguments of FUN. graph will be internally used as a mandatory argument (hopefully being called graph is consistent throughout functions...)
#'
#' @return output of FUN, or list of outputs
#' @export
#' @importFrom igraph is.igraph
#'
#' @examples
#' # set.seed(42)
#' n<- 10;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[upper.tri(Adj)]<- sample(0:8,n*(n-1)/2,replace = TRUE)
#' Adj[lower.tri(Adj)]<- sample(0:8,n*(n-1)/2,replace = TRUE)
#' Adj
#'
#' effort<- rowSums(Adj)+sample(5:8,n,replace = TRUE)
#' effort
#'
#' boot<- KenNet::KenBoot(Adj,effort,10,replacement = TRUE,proportion = 1.0,
#'         index = "focal",type = "total",mode = "plus",output = "graph")
#' metric.boot(graph = boot[[2]],igraph::transitivity,args.list = list(type="global"))
#' metric.boot(boot,igraph::strength,args.list = list(mode="total",loops=FALSE))

metric.boot<- function(graph,FUN,args.list){
  if(igraph::is.igraph(graph)){
    do.call(what = FUN,args = c(list(graph=graph),args.list))
  }
  else{
    if(is.list(graph)&igraph::is.igraph(graph[[1]])){
      lapply(graph,
             function(g) do.call(what = FUN,args = c(list(graph=g),args.list)))
    }
    else{
      stop("Unrecognized structure of graph argument")
    }
  }
}

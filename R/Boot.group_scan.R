#' Bootstrap Adjacency
#' Perform a bootstrap of group scans with probabilities derived from a provided adjancecy matrix, to produce a new adjancecy matrix.
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
#' @param n.cores number of threads to use while performingh the bootstrap
#'
#' @return according to output a list of scans, or the bootstrapped adjacency matrix
#' @export
#' @importFrom parallel detectCores
#' @importFrom snow makeCluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach `%dopar%`
#' @importFrom foreach foreach
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:30,n*(n-1),replace = TRUE)
#' Adj
#'
#' Boot.group_scan(Adj,42,output="adj")
#' Boots<- lapply(
#'   1:5,
#'     function(b){
#'         Boot.group_scan(Adj,42,mode = "directed",output = "adj",n.core=1)
#'     }
#' )

Boot.group_scan<- function(Adj,total_scan,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),
                           output=c("list","adjacency"),n.cores=(parallel::detectCores()-1)){
  output<- match.arg(output)
  mode<- match.arg(mode)

  cl<- snow::makeCluster(n.cores)
  doSNOW::registerDoSNOW(cl)
  boot.list<- foreach::`%dopar%`(foreach::foreach(b=1:total_scan,.export = c("do.scan","non.diagonal","Binary.prob")),
                                 do.scan(Adj = Adj,total_scan = total_scan,mode = mode)
  )
  snow::stopCluster(cl)

  switch(output,
         "list" = boot.list,
         "adjacency" = Reduce("+",boot.list))
}

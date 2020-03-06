#' Bootstrap scans
#' Bootstrap iterated binary group or focal scans with probabilities derived from a provided adjancecy matrix.
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param n.boot integer, number of bootstrap to perform.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#' @param focal.list Character vector, indicate the list of focals to consider throughout the scans.
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
#' @param n.cores number of threads to use while performingh the bootstrap
#' @param cl Optional cluster object (cf snow package), otherwise created according to n.cores.
#'
#' @return according to output and method: a list of iterated scans, or of adjacency matrix
#'
#' @export
#' @importFrom parallel detectCores
#' @importFrom snow makeCluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom pbapply pblapply
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
#' Adj
#'
#' focal.list<- sample(nodes,42,replace = TRUE)
#' table(focal.list)
#'
#' Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,scaled = TRUE,
#'            method = "group",mode = "directed",output = "list")
#' Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,scaled = TRUE,
#'            method = "focal",mode = "directed",output = "adj")
#' # system.time( #single threaded
#' #   Boot_scans(Adj,10,total_scan = 42,focal.list = focal.list,scaled = TRUE,
#' #              method = "both",mode = "directed",output = "adj",n.cores = 1)
#' # )
#' # system.time( #multi threaded (with number of thread existing - 1)
#' #   Boot_scans(Adj,10,total_scan = 42,focal.list = focal.list,scaled = TRUE,
#' #             method = "both",mode = "directed",output = "adj")
#' # )
Boot_scans<- function(Adj,n.boot,total_scan,
                      method=c("group","focal","both"),focal.list=NULL,scaled=FALSE,
                      mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),
                      output=c("list","adjacency"),n.cores=(parallel::detectCores()-1),cl=NULL){
  b<-NULL; #irrelevant bit of code, only to remove annoying note in R CMD Check...
  method<- match.arg(method)
  output<- match.arg(output)
  mode<- match.arg(mode)

  if(is.null(focal.list)){
    focal.list<- sample(rownames(Adj),total_scan,replace=TRUE)
  }else{
    if(length(focal.list)!=total_scan) {stop("Provided focal.list's length doesn't match total number of scans to perform.")}
  }

  if(is.null(cl)){
    cl<- snow::makeCluster(n.cores)
    doSNOW::registerDoSNOW(cl);on.exit(snow::stopCluster(cl))
  }

  pbapply::pblapply(
    1:n.boot,
    function(b){
      iterate_scans(Adj = Adj,total_scan = total_scan,focal.list = focal.list,scaled = scaled,
                    method = method,mode = mode,output = output,n.cores = n.cores,cl=cl)
    }
  )
}

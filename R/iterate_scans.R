#' Iterate scans
#' Internal use in Boot_scan. Iterate several binary group or focal scans with probabilities derived from a provided adjancecy matrix, to produce a new adjancecy matrix.
#'
#' @param Adj square integers matrix of occurences of dyads. Optional if using presence.prob. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param presence.prob. square probability matrix of presence (as in Bernouilli trial) of dyads. Optional if using Adj and total_scan.
#' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
#' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
##' @param ... additional argument to be used, to use produce a scan in a desired way.
#' @param focal.list Character vector, indicate the list of focals to consider throughout the scans.
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument
#'
#' @return according to output and method: a list of iterated scans, or an adjacency matrix
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
#' n<- 5;nodes<- as.character(1:n);total_scan<- 42;total_scan.rare<- 4200;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
#' Adj
#' Adj.rare<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj.rare[non.diagonal(Adj.rare)]<- sample(0:10,n*(n-1),replace = TRUE)
#' Adj.rare
#'
#' presence.prob<- Binary.prob(Adj,50)
#' obs.prob<- matrix(runif(n*n,0,1),n,n);diag(obs.prob)<- 0
#'
#' focal.list<- sample(nodes,total_scan,replace = TRUE)
#' table(focal.list)
#' focal.list.rare<- sample(nodes,total_scan.rare,replace = TRUE)
#' table(focal.list.rare)
#'
#' iterate_scans(Adj,total_scan,scaled = FALSE,method = "both",output = "list",obs.prob = 0.8)
#' iterate_scans(Adj,total_scan,scaled = TRUE,method = "group",output = "adjacency",obs.prob = obs.prob)
#' iterate_scans(Adj,total_scan,scaled = FALSE,method = "focal",output = "adjacency")
#' iterate_scans(Adj,total_scan,focal.list = focal.list,scaled = TRUE,obs.prob = 0.7,method = "both",mode = "directed",output = "all")
#' iterate_scans(Adj.rare,total_scan.rare,scaled = FALSE,method = "both",output = "list",obs.prob = 0.8,use.rare.opti = TRUE)
#' iterate_scans(Adj.rare,total_scan.rare,scaled = TRUE,method = "group",output = "adjacency",obs.prob = obs.prob,use.rare.opti = TRUE)
#' iterate_scans(Adj.rare,total_scan.rare,scaled = FALSE,method = "focal",output = "adjacency",use.rare.opti = TRUE)
#' iterate_scans(Adj.rare,total_scan.rare,focal.list = focal.list.rare,scaled = TRUE,obs.prob = 0.7,method = "both",mode = "directed",output = "all",use.rare.opti = TRUE)

iterate_scans<- function(Adj=NULL,total_scan,method=c("theoretical","group","focal","both"),
                         output=c("list","adjacency","all"),scaled=FALSE,...,
                         use.rare.opti = FALSE){
  scan.default.args(Adj = Adj,total_scan = total_scan,method = method,...)
  if(!use.rare.opti){
    scan_list<- lapply(
      1:total_scan,
      function(i){
        do.scan(presence.prob = presence.prob,method = method,focal = focal.list[i],obs.prob = obs.prob,Adj.subfun = Adj.subfun)
      }
    )
  }else{
    scan_list<- simulate_zeros.non.zeros(total_scan,presence.prob)
    zero<- attr(scan_list,"n.zeros")
    scan_list<- lapply(
      seq_along(scan_list),
      function(i){
        do.non.zero.scan(presence.prob = presence.prob,method = method,focal = focal.list[i],obs.prob = obs.prob,Adj.subfun = Adj.subfun)
      }
    )
    attr(scan_list,"n.zeros")<- zero
  }
  switch(output,
         "list" = scan_list,
         "adjacency" = sum_up.scans(scan_list = scan_list,scaled = scaled,method = method,mode = mode),
         "all" = list(
           list = scan_list,
           adjacency = sum_up.scans(scan_list = scan_list,scaled = scaled,method = method,mode = mode)
         )
  )
}

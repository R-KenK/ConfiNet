#' Bootstrap scans
#' Bootstrap iterated binary group or focal scans with probabilities derived from a provided adjancecy matrix.
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param n.boot integer, number of bootstrap to perform.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.
#' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
#' @param focal.list Character vector, indicate the list of focals to consider throughout the scans.
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param n.cores number of threads to use while performingh the bootstrap
#' @param cl Optional cluster object (cf snow package), experimentally set to put the makeCluster and stopCluster out of the bootable function. (WIP, next implementation should rethink this).
#' @param ... additional argument to be used, to use produce a scan in a desired way.
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation obs.probability matrix"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a general dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
#' @param n.cores number of threads to use while performingh the bootstrap
#' @param cl Optional cluster object (cf snow package), otherwise created according to n.cores.
#' @param use.rare.opti logical: should the optimization for rare event be used? If left NULL, choice is made automatically by decide_use.rare.opti().
#'
#' @return according to output and method: a list of iterated scans, or of adjacency matrix, with attributes to keep track of certain data
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
#'            method = "group",use.rare.opti=FALSE,mode = "directed",output = "list")
#' Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,scaled = TRUE,
#'            method = "focal",mode = "directed",output = "adj")
#' # system.time( #single threaded
#' #   Boot_scans(Adj,10,total_scan = 42,focal.list = focal.list,scaled = TRUE,
#' #              method = "both",mode = "directed",output = "adj",n.cores = 1)
#' # )
#' # system.time( #multi threaded (with number of thread existing - 1)
#' #   Boot_scans(Adj,10,total_scan = 42,focal.list = focal.list,
#' #             scaled = TRUE,obs.prob=0.7,keep=TRUE,
#' #             method = "both",mode = "directed",output = "adj")
#' # )
Boot_scans<- function(Adj,total_scan,method=c("theoretical","group","focal","both"),focal.list=NULL,n.boot,...,
                      scaled=FALSE,obs.prob=1,pbapply=FALSE,
                      mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector"),
                      output=c("list","adjacency","all"),use.rare.opti=NULL,n.cores=(parallel::detectCores()-1),cl){
  b<-NULL; #irrelevant bit of code, only to remove annoying note in R CMD Check...
  method<- match.arg(method)
  output<- match.arg(output)
  mode<- match.arg(mode)

  scan.default.args(Adj = Adj,total_scan = total_scan,method = method,...)

  Adj.subfun<- switch(mode,
                      "directed" = ,
                      "undirected" = ,
                      "max" = ,
                      "min" = ,
                      "plus" = non.diagonal,
                      "upper" = upper.tri,
                      "lower" =  lower.tri,
                      "vector" = function(V){rep(TRUE,length(V))}
  )
  presence.prob<- Binary.prob(Adj=Adj,total_scan=total_scan,mode = mode)

  if(is.null(use.rare.opti)){
    use.rare.opti<- FALSE
    # n<- nrow(Adj)
    # use.rare.opti<- decide_use.rare.opti(n = n*(n-1),total_scan = total_scan,max.obs = max(Adj))
  }
  if(pbapply){
    if(missing(cl)) {
      .export<- c("iterate_scans","do.non.zero.scan","do.scan","non.diagonal",
                  "sum_up.scans","sum_up.scan","matrix_sum_na.rm","quick.sample",
                  "Binary.prob","binary_adjacency_mode","scan.default.args","adjacency_mode",
                  "n.observed_edges","observable_edges","adjust.conditional.prob","focal.scan")
      cl<- snow::makeCluster(n.cores);
      snow::clusterExport(cl,list = .export);snow::clusterExport(cl,list = ls(sys.frame(which = 1)),envir = sys.frame(which = 1));
      on.exit(snow::stopCluster(cl));
    }
    Bootstrap<- pbapply::pblapply(
      1:n.boot,
      function(b){
        iterate_scans(Adj = Adj,total_scan = total_scan,presence.prob = presence.prob,
                      focal.list = focal.list,scaled = scaled,obs.prob = obs.prob,
                      method = method,mode = mode,output = output,use.rare.opti = use.rare.opti,Adj.subfun = Adj.subfun)
      },cl = cl
    )
  }else{
    Bootstrap<- lapply(
      1:n.boot,
      function(b){
        iterate_scans(Adj = Adj,total_scan = total_scan,presence.prob = presence.prob,
                      focal.list = focal.list,scaled = scaled,obs.prob = obs.prob,
                      method = method,mode = mode,output = output,use.rare.opti = use.rare.opti,Adj.subfun = Adj.subfun)
      }
    )
  }
  Bootstrap_add.attributes(Bootstrap = Bootstrap,method = method,scaled = scaled,
                           mode = mode,output = output,total_scan = total_scan)
}

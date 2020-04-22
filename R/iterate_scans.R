#' Iterate scans
#' Internal use in Boot_scan. Iterate several binary group or focal scans with probabilities derived from a provided adjancecy matrix, to produce a new adjancecy matrix.
#'
#' @param Adj square integers matrix of occurences of dyads. Optional if using presence.prob. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param presence.prob. square probability matrix of presence (as in Bernouilli trial) of dyads. Optional if using Adj and total_scan.
#' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
#' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param n.cores number of threads to use while performingh the bootstrap
#' @param cl Optional cluster object (cf snow package), experimentally set to put the makeCluster and stopCluster out of the bootable function. (WIP, next implementation should rethink this).
#' @param ... additional argument to be used, to use produce a scan in a desired way.
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
#' n<- 5;nodes<- as.character(1:n);
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
#' Adj
#'
#' presence.prob<- Binary.prob(Adj,50)
#' obs.prob<- matrix(runif(n*n,0,1),n,n);diag(obs.prob)<- 0
#'
#' focal.list<- sample(nodes,42,replace = TRUE)
#' table(focal.list)
#'
#' iterate_scans(Adj,42,scaled = FALSE,method = "group",output = "list",obs.prob = 0.8,n.cores = 1)
#' iterate_scans(Adj,42,scaled = TRUE,method = "group",output = "adjacency",obs.prob = obs.prob,n.cores = 1)
#' iterate_scans(Adj,42,scaled = FALSE,method = "focal",output = "adjacency",n.cores = 1)
#' iterate_scans(Adj,42,focal.list = focal.list,scaled = TRUE,obs.prob = 0.7,method = "both",mode = "directed",output = "all")

iterate_scans<- function(Adj=NULL,total_scan,method=c("theoretical","group","focal","both"),
                         output=c("list","adjacency","all"),scaled=FALSE,...,
                         n.cores=(parallel::detectCores()-1),cl){
  scan.default.args(Adj = Adj,total_scan = total_scan,method = method,...)

  if(missing(cl)) {
    .export<- c("do.scan","non.diagonal","quick.sample","Binary.prob","binary_adjacency_mode","scan.default.args","observable_edges");
    cl<- snow::makeCluster(n.cores);
    snow::clusterExport(cl,list = .export);snow::clusterExport(cl,list = ls(sys.frame(which = 1)),envir = sys.frame(which = 1));
    doSNOW::registerDoSNOW(cl);on.exit(snow::stopCluster(cl))
  }

  switch(method,
         "theoretical" = ,
         "group" = {
           scan_list<- pbapply::pblapply(
             1:total_scan,
             function(i){
               do.scan(presence.prob = presence.prob,method = method,obs.prob = obs.prob,Adj.subfun = Adj.subfun,check.defaults = FALSE)
             },cl = cl
           )
         },
         "focal" = {
           scan_list<- pbapply::pblapply(
             1:total_scan,
             function(i){
               do.scan(presence.prob = presence.prob,method = "focal",focal = focal.list[i],Adj.subfun = Adj.subfun)
             },cl = cl
           )
         },
         "both" = {
           scan_list<- pbapply::pblapply(
             1:total_scan,
             function(i){
               do.scan(presence.prob = presence.prob,method = "both",obs.prob = obs.prob,focal = focal.list[i],Adj.subfun = Adj.subfun)
             },cl = cl
           )
         }
  )

  switch(output,
         "list" = scan_list,
         "adjacency" = sum_up.scans(scan_list = scan_list,scaled = scaled,method = method,mode = mode),
         "all" = list(
           list = scan_list,
           adjacency = sum_up.scans(scan_list = scan_list,scaled = scaled,method = method,mode = mode)
         )
  )
}

# #' Title
# #'
##' @param P
##' @param N
##' @param ...
##' @param n
##'
##' @return
##' @export
##'
##' @examples
# iterate_rare.scans<- function(P,N,...,n = length(P)){
#   scan[non.zero,]<- data.table(rbind_lapply(seq_len(nrow(scan[non.zero,])),function(k) do.non.zero.scan(P)))
#   scan
# }
#
# iterate_rare.scans<- function(Adj,total_scan,method=c("group","focal","both"),
#                               focal.list=NULL,scaled=FALSE,obs.prob=NULL,keep=FALSE,
#                               mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),
#                               output=c("list","adjacency","all"),n.cores=(parallel::detectCores()-1),cl=NULL,Adj.subfun=NULL,prob=NULL){
#   b<-NULL; #irrelevant bit of code, only to remove annoying note in R CMD Check...
#   if(is.null(cl)) {cl<- snow::makeCluster(n.cores);doSNOW::registerDoSNOW(cl);on.exit(snow::stopCluster(cl))} # left to avoid error if the function is used alone, but should probably be used internally from Boot_scans() now.
#
#   switch(method,
#          "group" = {
#            scan<- data.table(matrix(0,N,n));
#            non.zero<- rbinom(N,1,1-prod(1-P))==1;
#
#            scan_list<- foreach::`%dopar%`(
#              foreach::foreach(b=1:total_scan,
#                               .export = c("do.scan","non.diagonal","Binary.prob","observable_edges","binary_adjacency_mode")),
#              do.non.zero.scan(Adj = Adj,total_scan = total_scan,
#                               focal = NULL,obs.prob=obs.prob,keep=keep,
#                               mode = mode, method = "group",Adj.subfun = Adj.subfun,prob = prob)
#            )
#          },
#          "focal" = {
#            scan_list<- foreach::`%dopar%`(
#              foreach::foreach(b=1:total_scan,
#                               .export = c("do.scan","non.diagonal","Binary.prob","binary_adjacency_mode")),
#              do.scan(Adj = Adj,total_scan = total_scan,
#                      focal = focal.list[b],obs.prob=NULL,keep=FALSE,
#                      mode = mode,method = "focal",Adj.subfun = Adj.subfun,prob = prob)
#            )
#          },
#          "both" = {
#            scan_list<- foreach::`%dopar%`(
#              foreach::foreach(b=1:total_scan,
#                               .export = c("do.scan","non.diagonal","Binary.prob","observable_edges","binary_adjacency_mode")),
#              do.scan(Adj = Adj,total_scan = total_scan,
#                      focal = focal.list[b],obs.prob=obs.prob,keep=keep,
#                      mode = mode,method = "both",Adj.subfun = Adj.subfun,prob = prob)
#            )
#          }
#   )
#
#   switch(output,
#          "list" = return(scan_list),
#          "adjacency" = return(sum_up.scans(Adj = Adj,scan_list = scan_list,keep = keep,scaled = scaled,method = method,mode = mode)),
#          "all" = return(
#            list(
#              list = scan_list,
#              adjacency = sum_up.scans(Adj = Adj,scan_list = scan_list,keep = keep,scaled = scaled,method = method,mode = mode)
#            )
#          )
#   )
# }

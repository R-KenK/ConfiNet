#' Iterate scans
#' Internal use in Boot_scan. Iterate several binary group or focal scans with probabilities derived from a provided adjancecy matrix, to produce a new adjancecy matrix.
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#' @param focal.list Character vector, indicate the list of focals to consider throughout the scans.
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation obs.probability matrix"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a general dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param keep logical. Relevant if group scans are performed. Indicate if the original "theoretical" group scan should be kept track of.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
#' @param n.cores number of threads to use while performingh the bootstrap
#' @param cl Optional cluster object (cf snow package), experimentally set to put the makeCluster and stopCluster out of the bootable function. (WIP, next implementation should rethink this).
#'
#' @return according to output and method: a list of iterated scans, or an adjacency matrix
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
#' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
#' Adj
#'
#' focal.list<- sample(nodes,42,replace = TRUE)
#' table(focal.list)
#' iterate_scans(Adj,42,scaled = FALSE,method = "group",output = "adjacency",n.cores = 1)
#' iterate_scans(Adj,42,focal.list,scaled = TRUE,obs.prob = 0.7,keep=TRUE,
#'               method = "focal",mode = "directed",output = "list")

iterate_scans<- function(Adj,total_scan,method=c("group","focal","both"),
                         focal.list=NULL,scaled=FALSE,obs.prob=NULL,keep=FALSE,
                         mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),
                         output=c("list","adjacency","all"),n.cores=(parallel::detectCores()-1),cl=NULL){
  b<-NULL; #irrelevant bit of code, only to remove annoying note in R CMD Check...
  if(is.null(cl)) {cl<- snow::makeCluster(n.cores);doSNOW::registerDoSNOW(cl);on.exit(snow::stopCluster(cl))} # left to avoid error if the function is used alone, but should probably be used internally from Boot_scans() now.

  switch(method,
         "group" = {
           scan_list<- foreach::`%dopar%`(
             foreach::foreach(b=1:total_scan,
                              .export = c("do.scan","non.diagonal","Binary.prob","observable_edges")),
             do.scan(Adj = Adj,total_scan = total_scan,
                     focal = NULL,obs.prob=obs.prob,keep=keep,
                     mode = mode, output = "group")
           )
         },
         "focal" = {
           scan_list<- foreach::`%dopar%`(
             foreach::foreach(b=1:total_scan,
                              .export = c("do.scan","non.diagonal","Binary.prob")),
             do.scan(Adj = Adj,total_scan = total_scan,
                     focal = focal.list[b],obs.prob=NULL,keep=FALSE,
                     mode = mode,output = "focal")
           )
         },
         "both" = {
           scan_list<- foreach::`%dopar%`(
             foreach::foreach(b=1:total_scan,
                              .export = c("do.scan","non.diagonal","Binary.prob","observable_edges")),
             do.scan(Adj = Adj,total_scan = total_scan,
                     focal = focal.list[b],obs.prob=obs.prob,keep=keep,
                     mode = mode,output = "both")
           )
         }
  )

  switch(output,
         "list" = return(scan_list),
         "adjacency" = return(sum_up.scans(Adj = Adj,scan_list = scan_list,keep = keep,scaled = scaled,method = method)),
         "all" = return(
           list(
             list = scan_list,
             adjacency = sum_up.scans(Adj = Adj,scan_list = scan_list,keep = keep,scaled = scaled,method = method)
           )
         )
  )
}

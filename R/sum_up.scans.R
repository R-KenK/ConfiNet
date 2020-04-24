#' Sum up and flatten list of binary scans
#' Internal use. Sum up binary adjacency matrices or focal vectors according to the method used. sum_up.scan*s*() is a wrapper for sum_up.scan() that returns summed_up scans for the differnet requested methods0
#'
#' @param scan_list list of binary adjacency matrices, binary focal vectors, or both, outputed by iterating do.scan()
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param use.rare.opti logical: should the optimization for rare event be used?
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return An integer adjacency matrix, or a list of two if method = "both". Now returns NAs for dyad scaled with no observation at all.
#' @export
#'
#' @examples
#' #Internal use for readability
sum_up.scans<- function(scan_list,scaled=FALSE,use.rare.opti=FALSE,obs.prob = NULL,
                        method = c("theoretical","group","focal","both"),mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
  method<- match.arg(method);
  switch(method,
         "theoretical" = {
           list(
             theoretical = sum_up.scan(scan_list,"theoretical",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob)
           )
         },
         "group" = {
           list(
             theoretical = sum_up.scan(scan_list,"theoretical",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob),
             group = sum_up.scan(scan_list,"group",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob)
           )
         },
         "focal" = {
           list(
             theoretical = sum_up.scan(scan_list,"theoretical",scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob),
             focal = sum_up.scan(scan_list,"focal",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob)
           )
         },
         "both" = {
           list(
             theoretical = sum_up.scan(scan_list,"theoretical",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob),
             group = sum_up.scan(scan_list,"group",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob),
             focal = sum_up.scan(scan_list,"focal",mode = mode,scaled = scaled,use.rare.opti = use.rare.opti,obs.prob = obs.prob)
           )
         }
  )
}

#' Sum up and flatten list of binary scans
#' Internal use. Sum up binary adjacency matrices or focal vectors according to the method used. sum_up.scan*s*() is a wrapper for sum_up.scan() that returns summed_up scans for the differnet requested methods0
#'
#' @param scan_list list of binary adjacency matrices, binary focal vectors, or both, outputed by iterating do.scan()
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param use.rare.opti logical: should the optimization for rare event be used?
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#'
#' @return An integer adjacency matrix, or a list of two if method = "both". Now returns NAs for dyad scaled with no observation at all.
#' @export
#'
#' @examples
#' #Internal use for readability
sum_up.scan<- function(scan_list,method,mode = NULL,scaled,use.rare.opti = FALSE,obs.prob = NULL){
  if(use.rare.opti) {n.zeros<- attr(scan_list,"n.zeros")}
  scan_list.method<- lapply(scan_list,function(scan) scan[[method]]);
  summed_up.method<- Reduce(matrix_sum_na.rm,scan_list.method)
  if(scaled){
    if(method=="theoretical"){
      summed_up.method<- summed_up.method/length(scan_list)
    }else{
      summed_up.method<- summed_up.method/n.observed_edges(scan_list.method,diag = 1,use.rare.opti = use.rare.opti,obs.prob = obs.prob,n.zeros = n.zeros) # here diag = 1 because while the count of the diagonal is irrelevant, it shouldn't be 0/0.
      summed_up.method<- ifelse(!is.nan(summed_up.method),summed_up.method,NA)
    }
  }
  adjacency_mode(summed_up.method,mode = mode)
}

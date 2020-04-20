#' Sum up and flatten list of binary scans
#' Internal use. Sum up binary adjacency matrices or focal vectors according to the method used
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param scan_list list of binary adjacency matrices, binary focal vectors, or both, outputed by iterating do.scan()
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return An integer adjacency matrix, or a list of two if method = "both"
#' @export
#'
#' @examples
#' #Internal use for readability

sum_up.scans<- function(Adj,scan_list,scaled=FALSE,
                        method = c("theoretical","group","focal","both"),mode = c("directed", "undirected", "max","min", "upper", "lower", "plus","vector")){
  method<- match.arg(method);total_scan<- length(scan_list)

  scan_list.theoretical<- lapply(scan_list,function(l) l$theoretical);
  summed_up.theoretical<- Reduce(matrix_sum_na.rm,scan_list.theoretical)

  if(method %in% c("group","both")){
    scan_list.group<- lapply(scan_list,function(l) l$group);
    summed_up.group<- Reduce(matrix_sum_na.rm,scan_list.group)

    if(scaled){
      summed_up.theoretical<- summed_up.theoretical/total_scan
      summed_up.group<- summed_up.group/n.observed_edges(scan_list.group,diag = 1) # here and hereafter diag = 1 because while the count of the diagonal is irrelevant, it shouldn't be 0/0.
    }
    adj.group<- adjacency_mode(summed_up.group,mode = mode)
  }

  if(method %in% c("focal","both")){
    scan_list.focal<- lapply(scan_list,function(l) l$focal);
    focal.list<- sapply(scan_list.focal,function(scan) attr(scan,"focal"))
    summed_up.focal<- rbind_lapply(rownames(Adj),
                                   function(node){
                                     scan_list.focal.node<- scan_list.focal[focal.list==node]
                                     if(length(scan_list.focal.node)==0) {
                                       scan_list.focal.node<- rep(0,ncol(Adj))
                                     }else{
                                       Reduce("+",scan_list.focal.node)
                                     }
                                   }
    )
    if(scaled){
      summed_up.theoretical<- summed_up.theoretical/total_scan
      summed_up.focal<- summed_up.focal/rbind_lapply(rownames(Adj),function(node) rep(length(scan_list.focal[focal.list==node]),ncol(Adj)))
    }
    adj.focal<- adjacency_mode(summed_up.focal,mode = mode)
  }

  switch(method,
         "theoretical" = {
           list(
             theoretical = summed_up.theoretical,
           )
         },
         "group" = {
           list(
             theoretical = summed_up.theoretical,
             group = adj.group
           )
         },
         "focal" = {
           list(
             theoretical = summed_up.theoretical,
             focal = adj.focal
           )
         },
         "both" = {
           list(
             theoretical = summed_up.theoretical,
             group = adj.group,
             focal = adj.focal
           )
         }
  )
}

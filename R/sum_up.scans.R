#' Sum up and flatten list of binary scans
#' Internal use. Sum up binary adjacency matrices or focal vectors according to the method used
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param scan_list list of binary adjacency matrices, binary focal vectors, or both, outputed by iterating do.scan()
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param keep logical. Relevant if group scans are performed. Indicate if the original "theoretical" group scan should be kept track of.
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#'
#' @return A non-binary adjacency matrix, or a list of two if method = "both"
#' @export
#'
#' @examples
#' #Internal use for readability

sum_up.scans<- function(Adj,scan_list,scaled=FALSE,keep=FALSE,
                        method = c("group","focal","both")){
  method<- match.arg(method)
  switch(method,
         "group" = {
           if(keep){
             scan_list.theoretical<- lapply(scan_list,function(l) l$theoretical);
             scan_list.observed<- lapply(scan_list,function(l) l$observed);
             list(
               theoretical = Reduce(matrix_sum_na.rm,scan_list.theoretical)/ifelse(scaled,n.observed_edges(scan_list.theoretical,diag = 1),1),
               group = Reduce(matrix_sum_na.rm,scan_list.observed)/ifelse(scaled,n.observed_edges(scan_list.observed,diag = 1),1)
             )
           }else{
             Reduce(matrix_sum_na.rm,scan_list)/ifelse(scaled,n.observed_edges(scan_list,diag = 1),1)
           }
         },
         "focal" = {
           summed_up<- do.call(rbind,
                               lapply(rownames(Adj),
                                      function(node) {
                                        scan_list.focal<- scan_list[sapply(scan_list,function(b) attr(b,"focal"))==node]
                                        Reduce("+",scan_list.focal)/ifelse(scaled,length(scan_list.focal),1)
                                      }
                               )
           )
           row.names(summed_up)<- rownames(Adj)
           summed_up
         },
         "both" = {
           list(
             group = {
               scan_list.group<- lapply(scan_list,function(l) l$group);
               if(keep){
                 scan_list.theoretical<- lapply(scan_list.group,function(l) l$theoretical);
                 scan_list.observed<- lapply(scan_list.group,function(l) l$observed);
                 list(
                   theoretical = Reduce(matrix_sum_na.rm,scan_list.theoretical)/ifelse(scaled,n.observed_edges(scan_list.theoretical,diag = 1),1),
                   observed = Reduce(matrix_sum_na.rm,scan_list.observed)/ifelse(scaled,n.observed_edges(scan_list.observed,diag = 1),1)
                 )
               }else{
                 Reduce(matrix_sum_na.rm,scan_list)/ifelse(scaled,n.observed_edges(scan_list,diag = 1),1)
               }
             },
             focal = {
               summed_up<- do.call(rbind,
                                  lapply(rownames(Adj),
                                         function(node) {
                                           scan_list.focal<- scan_list[sapply(scan_list,function(l) attr(l$focal,"focal"))==node]
                                           scan_list.focal<- lapply(scan_list.focal,function(l) l$focal)
                                           Reduce("+",scan_list.focal)/ifelse(scaled,length(scan_list.focal),1)
                                         }
                                  )
               );
               row.names(summed_up)<- rownames(Adj)
               summed_up
             }
           )
         }
  )
}

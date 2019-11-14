#' Calculate network metric
#'
#' Wrapper for internal use.
#'
#' @param graph graph (igraph) for which the metric will be calculated
#' @param FUN function to calculate the desired metric
#' @param ... arguments of FUN
#'
#' @return output of FUN
#' @export
#'
#' @examples
#' print("internal use. Look around.")
metric.boot<- function(graph,FUN,...){
  FUN(graph,...)
}

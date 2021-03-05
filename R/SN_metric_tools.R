#' TO WRITE
#'
#' @param Adj TO WRITE
#' @param mode TO WRITE
#'
#' @return TO WRITE
#' @export

#' @importFrom DirectedClustering ClustF
#'
#' @examples
#' # TO WRITE
GlobalCC <- function(Adj,mode) {
  if (inherits(Adj,"igraph")) {
    Adj <- igraph::get.adjacency(Adj,type = "both",sparse = FALSE,attr = "weight")
  }
  switch(mode,
         "upper" = ,
         "lower" = {
           Adj <- Adj + t(Adj)
           type <- "undirected"
         },
         "directed" = type <- "directed",
         type <- "undirected"
  )
  DirectedClustering::ClustF(Adj,type = type)[["GlobalCC"]]
}

#' Simulate a focal scan
#'
#' Simulate a list of instantaneous adjacency matrices
#'
#' @param n_focals integer, number of focals to simulate.
#' @param n_scans integer, number of scans to simulate.
#' @param n_nodes integer, number of nodes to simulate.
#' @param nodes_names optional character vector of the nodes' names, otherwise simply numbered.
#' @param focals optional character vector of the nodes' names, otherwise simply numbered.
#' @param prob optional character vector of the nodes' names, otherwise simply numbered.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param output Character scalar, specifies if the function should return a list of lists of scans (n_focals-list of n_scans vectors), a list of egocentric network (n_focals-list, can have several per node), or an adjacency matrix. Only adjacency matrices can be undirected (???).
#'
#' @return depends on output
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' simu.focal(n_focals = 15,n_scans = 5,n_nodes = 10,mode = "directed",output = "scan")
#' simu.focal(n_focals = 15,n_scans = 5,n_nodes = 10,mode = "directed",output = "ego")
#' simu.focal(n_focals = 500,n_scans = 5,n_nodes = 10,mode = "directed",output = "adj")
#' simu.focal(n_focals = 15,n_scans = 5,n_nodes = 10,mode = "plus",output = "adj")
simu.focal<- function(n_focals,n_scans,n_nodes, nodes_names=as.character(1:n_nodes),focals=nodes_names,prob=NULL,
                      mode = c("directed", "undirected", "plus"),
                      output = c("scans","egos","adjacency")){
  if(length(mode)>1) {mode<- "directed"}
  if(length(output)>1) {output<- "scans"}

  if((mode %in% c("undirected","plus")) & !(output %in% c("adj","Adj","adjacency"))) {stop("Incompatible mode and output provided")}

  if(identical(focals,nodes_names)) {focals<- sample(nodes_names,n_focals,replace = TRUE)}
  if(length(focals)!=n_focals) {
    warning("focal list provided not the same length of provided n_focals. n_focals argument ignored")
    n_focals<- length(focals)
  }

  Focal.scans<- lapply(focals,
                       function(focal) {
                         lapply(1:n_scans,
                                function(scan) {
                                  ego<- rep(NA,n_nodes);names(ego)<- nodes_names;
                                  ego[names(ego)!=focal]<- sample(0:1,n_nodes-1,replace = TRUE,prob = prob)
                                  ego
                                }
                         )
                       }
  )
  names(Focal.scans)<- focals
  if(output %in% c("scans","scan")) return(Focal.scans)

  Combined.focals<- lapply(1:n_focals,
                           function(foc) {
                             ego.scans<- ifelse(nodes_names==focals[foc],NA,0);names(ego.scans)<- nodes_names;
                             for(scan in 1:n_scans){
                               ego.scans[nodes_names!=focals[foc]]<- ego.scans[nodes_names!=focals[foc]]+Focal.scans[[foc]][[scan]][nodes_names!=focals[foc]]
                             }
                             ego.scans
                           }
  )
  if(output %in% c("egos","ego","egonetwork")) {
    names(Combined.focals)<- focals;
    return(Combined.focals)
    }


  Adj.focals<- do.call(rbind,
                       lapply(nodes_names,
                              function(node){
                                Combined.ego<- Combined.focals[sapply(Combined.focals,function(foc) which(is.na(foc)))==node];
                                Combined.ego<- Reduce("+",Combined.ego)
                                if(!is.null(Combined.ego)) Combined.ego
                                else ifelse(nodes_names==node,NA,0)
                              }
                       )
  )
  rownames(Adj.focals)<- nodes_names;
  diag(Adj.focals)<- 0;

  if(output %in% c("adj","Adj","adjacency")) {
    switch(mode,
           "directed" = return(Adj.focals),
           "undirected" =,
           "plus" = return(Adj.focals+t(Adj.focals)), # CAN UNDIRECTED BE DIFFERENTLY DEFINED THAN PLUS?
           stop("wtf with the mode?")
    )
  }
  stop("wtf with the output?")
}

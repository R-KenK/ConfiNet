#' Simulate a group scan
#'
#' Simulate a list of instantaneous adjacency matrices
#'
#' @param n_scans integer, number of scans to simulate.
#' @param n_nodes integer, number of nodes to simulate.
#' @param nodes_names optional character vector of the nodes' names, otherwise simply numbered.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return list of adjacency matrices
#'
#' @examples
#' set.seed(42)
#' n<- 10L;nodes<- letters[1:n];n_scans<- 50L;
#'
#' simu.scan(n_scans,n_nodes = n,mode = "undirected")
#' simu.scan(n_scans,n_nodes = n,mode = "directed",output="scan")
#' simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "plus",output="adj")
#'
simu.scan<- function(n_scans,n_nodes, nodes_names=as.character(1:n_nodes),prob=NULL,
                     mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),
                     output=c("scan","adjacency")){
  if(length(mode)>1) {mode<- "directed"}
  if(length(output)>1) {output<- "scan"}
  if(is.null(prob)) {
    prob<- matrix(data = 1/n_nodes,nrow = n_nodes,ncol = n_nodes,dimnames = list(nodes_names,nodes_names));
    diag(prob)<- 0
  }else{
    if(!is.matrix(prob)) stop("dyadic probability should be provided as a square matrix")
    if(nrow(prob)!=ncol(prob)|nrow(prob)!=n_nodes) stop("dimensions of dyadic probability matrix incorrect or incompatible with specified n_nodes")
    if(is.null(row.names(prob))) {row.names(prob)<- nodes_names;colnames(prob)<- nodes_names;}
  }

  Scans<- lapply(1:n_scans,
                 function(scan) {
                  # do.scan(Adj=prob,mode=mode) # WIP
                 }
  )

  if(output %in% c("scan","scans","list")) return(Scans)

  if(output %in% c("adj","adjacency")) return(Reduce("+",Scans)) else stop("invalid output provided. Should be either \"scan\" or \"adjacency\".")
}

# set.seed(42)
# n<- 10L;nodes<- as.character(1:n);n_scans<- 10L;
#
#
# plot(sort(dnorm(1:n,mean = n/2,sd=n/5)))
#
# prob<- dnorm(1:n,mean = n/2,sd=n/5)
#
# prob.test.lin<- sapply(1:n,
#                    function(i){
#                      sapply(1:n,
#                             function(j){
#                               ifelse(i!=j,i*j,0)
#                             }
#                      )
#                    }
# )
#
# ggplot2::ggplot(lin.dt,ggplot2::aes(y,-x,colour=prob,size=prob))+ggplot2::geom_point()
# ggplot2::ggplot(norm.dt,ggplot2::aes(y,-x,colour=prob,size=prob))+ggplot2::geom_point()
#
#
#
# scale.to.binary.prob(prob.test.lin,n*(n-1),mode = "upper")
#
# prob<- dnorm(1:n,mean = n/2,sd=n/3)
# prob.norm<- sapply(1:n,
#                    function(i){
#                      sapply(1:n,
#                             function(j){
#                               ifelse(i!=j,prob[i]*prob[j],0)
#                             }
#                      )
#                    }
# )
#
# scale.to.binary.prob(prob.norm,max(prob.norm),mode = "max")
# scale.to.binary.prob(test.Adj,1138,mode = "max")
#
#
#
# identical(do.scan(prob.test.lin,n*(n-1),mode = "max"),do.scan(prob.test.lin,n*(n-1),mode = "plus"))
#
# prob.scaled<- prob.test.lin[upper.tri(prob.test.lin)]/(min(prob.test.lin[upper.tri(prob.test.lin)])+max(prob.test.lin[upper.tri(prob.test.lin)]))
# bin<- data.table::data.table(yes=prob.scaled,no=1-prob.scaled)
# test<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(as.character(1:n),as.character(1:n)))
#
# sapply(1:(n*(n-1)/2),
#        function(d) {
#          test[upper.tri(test)][d]<- sample(c(1,0),1,replace = TRUE,prob=bin[d,c("yes","no")])
#        }
# )
#
#
# test[upper.tri(test)]
#
# simu<- simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "undirected",output="adj")
#
# library(igraph)
# set.seed(42)
# n<- 25L;nodes<- as.character(1:n);n_scans<- 10L;
#
# prob.lin<- sapply(1:n,
#                        function(i){
#                          sapply(1:n,
#                                 function(j){
#                                   ifelse(i!=j,i*j,0)
#                                 }
#                          )
#                        }
# )
#
# prob<- dnorm(1:n,mean = n/2,sd=n/3)
# prob.norm<- sapply(1:n,
#                         function(i){
#                           sapply(1:n,
#                                  function(j){
#                                    ifelse(i!=j,prob[i]*prob[j],0)
#                                  }
#                           )
#                         }
# )
#
# simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "plus",output="scan")
# simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "directed",output="scan")
# simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "max",output="adj")
# simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "max",output="adj",prob = prob.lin)
# simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "max",output="adj",prob = prob.norm)
#
#
# SN.norm<- igraph::graph_from_adjacency_matrix(simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "max",output="adj",prob = prob.norm),
#                                               weighted = TRUE,mode = "max")
# plot(SN.norm,edge.width=1*(E(SN.norm)$weight)/4,
#      edge.alpha=E(SN.norm)$weight,
#      edge.color="grey90",
#      vertex.size=eigen_centrality(SN.norm)$vector*30)
#
#
#
# SN.lin<- igraph::graph_from_adjacency_matrix(simu.scan(n_scans,n_nodes = n,nodes_names = nodes,mode = "max",output="adj",prob = prob.lin),
#                                               weighted = TRUE,mode = "max")
# plot(SN.lin,edge.width=E(SN.lin)$weight,
#      edge.alpha=E(SN.lin)$weight,
#      edge.color="grey90",
#      vertex.size=eigen_centrality(SN.lin)$vector*25)
#
# plot(igraph::graph_from_adjacency_matrix(prob.test.norm,weighted = TRUE,mode = "max"))
#

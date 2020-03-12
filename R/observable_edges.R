#' Hide unobservable edges
#' Simulate that some dyads might reasonable not be observable during a group scan
#'
#' @param Scan a square matrix
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation obs.probability matrix"}{of same dimension as Scan}
#'  \item{"a dyad observation vector"}{subsetted similarly as Scan (through the non.diagonal() function for instance)}
#'  \item{"a general dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param keep logical. Relevant if group scans are performed. Indicate if the original "theoretical" group scan should be kept track of.
#'
#' @return a similar matrix as Scan where some non diagonal edges have a obs.probability to be NAs. If keep is TRUE, returns a list with theoretical and observed scan.
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 6;
#' Scan<- matrix(sample(c(1,0),6*6,replace=TRUE),n,n);diag(Scan)<- 0
#' obs.prob<- matrix(runif(6*6,0,1),n,n);diag(obs.prob)<- 0
#' obs.prob.single<- runif(1,0,1)
#'
#' observable_edges(Scan,obs.prob)
#' observable_edges(Scan,obs.prob.single)
#'

observable_edges<- function(Scan,obs.prob=NULL,keep=FALSE){
  if(is.matrix(obs.prob)) {obs.prob<- non.diagonal(obs.prob,"vector")}
  if(length(obs.prob)==1) {ifelse(obs.prob<=1 & obs.prob>=0,obs.prob<- rep(obs.prob,length(non.diagonal(Scan,"vector"))),stop("Single observation obs.probability provided should be within [0,1]"))}

  observable<- sapply(obs.prob,function(p) sample(c(TRUE,FALSE),1,prob = c(p,1-p)))

  observed<- Scan;
  observed[non.diagonal(observed)]<- ifelse(observable,non.diagonal(Scan,"vector"),NA)
  if(keep){
    list("theoretical" = Scan,
         "observed" = observed)
  }else{
    observed
  }
}
#
# set.seed(42)
#
# n<- 10;nodes<- as.character(1:n);
# total_scan<- 100; #from original paper
# n.boot<- 5;
#
# Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
# Adj[non.diagonal(Adj)]<- sample((0:round(total_scan*.50)),n*(n-1),replace = TRUE)
#
# dirty_EV(Adj,"undirected")
# dirty_EV(Adj,"undirected")
#
# EV<- dirty_EV(Adj,"undirected")
#
# obs.prob.trait.plus<- matrix(0,length(EV),length(EV),dimnames = list(row.names(Adj),row.names(Adj)))
# for(i in seq_along(row.names(Adj))){
#   for(j in seq_along(row.names(Adj))){
#     if(i!=j) {obs.prob.trait.plus[i,j]<- i+j}
#   }
# }
# obs.prob.trait.plus
#
# obs.prob.trait.prod<- matrix(0,length(EV),length(EV),dimnames = list(row.names(Adj),row.names(Adj)))
# for(i in seq_along(EV)){
#   for(j in seq_along(EV)){
#     if(i!=j) {obs.prob.trait.prod[i,j]<- i*j}
#   }
# }
# obs.prob.trait.prod
#
# obs.prob.net.plus<- matrix(0,length(EV),length(EV),dimnames = list(row.names(Adj),row.names(Adj)))
# for(i in seq_along(row.names(Adj))){
#   for(j in seq_along(row.names(Adj))){
#     if(i!=j) {obs.prob.net.plus[i,j]<- EV[i]+EV[j]}
#   }
# }
# obs.prob.net.plus
#
# obs.prob.net.prod<- matrix(0,length(EV),length(EV),dimnames = list(row.names(Adj),row.names(Adj)))
# for(i in seq_along(EV)){
#   for(j in seq_along(EV)){
#     if(i!=j) {obs.prob.net.prod[i,j]<- EV[i]*EV[j]}
#   }
# }
# obs.prob.net.prod
#
# G<- igraph::read_graph("C:/R/Git/asnr/Networks/Mammalia/rhesusmacaque_association_weighted/weighted_Contact_Sits_Macaque.graphml",format = "graphml")
# total_scan<- 1138;
# Adj<- as.matrix(igraph::as_adj(G,attr = "weight"))
# row.names(Adj)<- as.character(1:nrow(Adj));colnames(Adj)<- row.names(Adj)
# Adj
# sorta.default.plot(Adj,edge.with.mul = 1/100,vertex.size.mul = 20,centrality.fun = "EV")
#
# focal.list<- sample(nodes,total_scan,replace = TRUE)
#
# Scan<- do.scan(Adj,total_scan)
# obs.prob<- matrix(runif(10*10,0,1),10,10);diag(obs.prob)<- 0
#
# observable_edges(Scan,obs.prob = obs.prob.net.plus/2,keep = T)
#
# obs.test<- obs.prob.net.plus/2
#
# do.scan(Adj,total_scan,obs.prob = obs.test,keep = TRUE)
#
# iterate_scans(Adj,total_scan = 1138,focal.list = focal.list,
#               scaled = FALSE,keep = TRUE,obs.prob = 0.5,
#               method = "both",mode = "max",output = "all",n.cores = 7)
#
#
# dirty_EV<- function(Adj,mode = mode){
#   graph<- igraph::graph.adjacency(Adj,weighted = TRUE,diag = TRUE,mode = mode)
#   igraph::eigen_centrality(graph)$vector
# }

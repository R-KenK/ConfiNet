# default_explorative_net.plot<- function(graph,edge.with.mul=NULL,vertex.size.mul=NULL,centrality.fun=c("degree","strength","eigen_vector"),...){
#   if(!igraph::is.igraph(graph)) {
#     if(is.matrix(graph) & dim(graph)[1]==dim(graph)[2]){
#       graph<- igraph::graph_from_adjacency_matrix(graph,weighted = TRUE,mode = "plus")
#       message("Provided graph argument is not an igraph object, but is compatible with adjacency matrix structure and has been considered as such.")
#     } else stop("Invalid graph argument provided.")
#   }
#   if(centrality.fun %in% c("EV","ev")) {centrality.fun<- "eigen_vector"} else {centrality.fun<- match.arg(centrality.fun)}
#
#   Adj<- igraph::get.adjacency(graph,type = "both",attr = ifelse(igraph::is.weighted(graph),"weight",NULL),names = TRUE,sparse = FALSE)
#
#   if(is.null(edge.with.mul)) {edge.with.mul<-5*1/max(Adj)}
#   if(is.null(vertex.size.mul)) {vertex.size.mul<-3*1/max(Adj)}
#
#   centrality.fun<- switch(centrality.fun,
#                           "degree" = igraph::degree,
#                           "strength" = igraph::strength,
#                           "eigen_vector" = function(graph){igraph::eigen_centrality(graph)$vector}
#   )
#   plot(graph,edge.width=igraph::E(graph)$weight*edge.with.mul,
#        edge.alpha=igraph::E(graph)$weight*edge.with.mul,
#        edge.color="grey70",vertex.label=NA,
#        vertex.size=centrality.fun(graph)*vertex.size.mul,...)
# }

# Title
#
# @param old.path
# @param new.path
# @param correcting.factor
#
# @return
# @export
#
# @examples
# export_to_graphml<- function(old.path,new.path,correcting.factor){
#   Adj.corrected<- ceiling(import_from_graphml(old.path,"adjacency")*correcting.factor)
#   G.corrected<- igraph::graph.adjacency(Adj.corrected,weighted = TRUE,mode = "max")
#   igraph::write_graph(graph = G.corrected,file = new.path,format = "graphml")
# }

# Title
#
# @param dir_path
#
# @return
# @export
#
# @examples
# get.total_scan<- function(dir_path){
#   if(file.exists(paste0(dir_path,"/total_scan.R"))){
#     source(paste0(dir_path,"/total_scan.R"))
#   } else{
#     total_scan<- NULL;
#   }
#   total_scan
# }

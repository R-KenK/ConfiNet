#' Network Bootstrap of edge-lists
#'
#' Bootstrap edge lists where repeated lines represent dyads' weights
#'
#' @param Adj raw adjacency matrix
#' @param effort vector of sampling effort per individual. WIP
#' @param type type of sampling effort. cf \link[KenNet]{adj.to.assoc}
#' @param boot number of bootstrap iteration. WIP
#' @param replacement resampling with replacement or not. WIP
#' @param proportion proportion to sub-sample, within [0,1[. WIP
#' @param index Type of desired association index. At the molment only focal method implemented. WIP
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param output desired output of the bootstrap: edge lists, adjacency matrices, association index matrices, igraph graphs, or network metrics
#' @param FUN if output is "metric", indicate which function to use to calculate said metric.
#' @param args.list if output is "metric", arguments for function FUN.
#'
#' @return according to output, returns a list of edge lists, adjacency matrices, association index matrices, or igraph graphs
#' @export
#'
#' @examples
#' print("WIP")
KenBoot<- function(Adj,effort,type = c("total","yab"),boot=100,replacement=TRUE,proportion=1.0,
                   index=c("joint","sri","hw","tw","sqrt","socaff","both","focal"),
                   mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),
                   output=c("adjacency","edge","assoc","graph","metric"),
                   FUN=NULL, args.list){

  if(length(index)>1) {index<- "sri"}
  if(length(type)>1) {type<- "total"}
  if(length(mode)>1) {mode<- "plus"}

  if(nrow(Adj)!=ncol(Adj)) {stop("Adj is not a square matrix.")}
  if(length(effort)!=nrow(Adj)) {stop("effort is not the same length as the number of individuals in Adj")}
  if(!replacement&proportion>=1) stop("If resampling is done without replacement, please provide a proportion of data in [0,1[ to be sub-sampled.")

  if(!is.null(effort)) {
    if(length(type)>1) {
      if(all(rowSums(Adj)<=effort)) {
        type<- "total";
        message("provided effort assumed to be total effort. If it is not the case, please explicit the type argument (cf. ?adj.to.assoc).")
      }else{
        type<- "yab";
        message("provided effort assumed to be number of time an individual has been observed not associated (cf. ?adj.to.assoc).")
      }
    }

    switch(type,
           # "total" = if(all(rowSums(Adj)<=effort)) diag(Adj)<- effort-rowSums(Adj) else stop("effort provided smaller than number of associations. This shouldn't be possible."),
           "total" = diag(Adj)<- effort-rowSums(Adj),
           "yab" = diag(Adj)<- effort,
           "diag" = stop("\"diag\" type shouldn't be used for KenBoot, only for adj.to.assoc."),
           stop("unrecognized effort type.")

    )
  }

  edge<- KenNet::adj.to.edge(Adj,mode = "directed")

  if(output %in% c("edge","adjacency","assoc","graph","metric")){
    e.boot<- lapply(1:boot,function(b) {
      if(replacement)
        edge[sample(1:nrow(edge),nrow(edge),replace = TRUE)]
      else
        edge[sample(1:nrow(edge),round(nrow(edge)*proportion,0),replace = FALSE)]
    }
    )
    if(output %in% c("adjacency","assoc","graph","metric")){
      a.boot<- lapply(e.boot,function(edge) {
        KenNet::edge.to.adj(edge,mode = "directed")
      }
      )
      if(output %in% c("assoc","graph","metric")){
        assoc.boot<- lapply(a.boot,function(Adj) {
          KenNet::adj.to.assoc(Adj,effort = diag(Adj),index = "focal",type = "diag")
        }
        )
        if(output %in% c("graph","metric")){
          graph.boot<- lapply(assoc.boot,function(assoc) {
            igraph::graph_from_adjacency_matrix(assoc,mode = "plus",weighted = TRUE)
          }
          )
          if(output %in% c("metric")){
            netmet.boot<- lapply(graph.boot,function(graph) {
              KenNet::metric.boot(graph,FUN = FUN,args.list=args.list)
            }
            )
          }
        }
      }
    }
  }

  switch (output,
          "edge" = e.boot,
          "adjacency" = a.boot,
          "assoc" = assoc.boot,
          "graph" = graph.boot,
          "metric" = netmet.boot,
          stop("output type not recognized.")
  )
}

# set.seed(42)
#
# n<- 10;nodes<- letters[1:n];
# Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
# Adj[upper.tri(Adj)]<- sample(0:8,n*(n-1)/2,replace = TRUE)
# Adj[lower.tri(Adj)]<- sample(0:8,n*(n-1)/2,replace = TRUE)
# Adj
#
# effort<- rowSums(Adj)+sample(5:8,n,replace = TRUE)
# effort
#
# boot<- KenNet::KenBoot(Adj,effort,10,replacement = TRUE,proportion = 1.0,
#         index = "focal",type = "total",mode = "plus",output = "graph")
#
#
# edge<- KenNet::adj.to.edge(Adj)
# el<- KenBoot(edge,boot = 10,output = "edge")
#
# al<- lapply(el,
#             function(edge) {
#               ADJ<- matrix(0,nrow(Adj),ncol(Adj),dimnames = list(row.names(Adj.tmp),row.names(Adj.tmp)))
#               Adj.tmp<- edge.to.adj(edge,mode = "directed")
#               ADJ[rownames(Adj.tmp),colnames(Adj.tmp)]<- Adj.tmp
#               ADJ[sort(rownames(ADJ)),sort(colnames(ADJ))]
#             })
# al
#
# a<- al[[2]]
#
# for(x in a>0){
#   a[[x]]
# }
#
# pos<- non.zero.non.diag(a)
#
# pos[2,]["row"]
#
# a[pos]<- sapply(1:nrow(pos),
#                 function(ij){
#                   i<- pos[ij,]["row"];j<- pos[ij,]["col"];
#                   a[i,j]/(a[i,i]+a[j,j])
#                 })
# a
#
# str(pos)
# Matrix::
#   a[a>0]
#
# ef<- rowSums(a)
#
#
# adj.mat.ij<- matrix(0,nrow = length(IDs),ncol = length(IDs));                                       #make matrices with the right dimensions as place-holders
# row.names(adj.mat.ij)<- IDs;colnames(adj.mat.ij)<- IDs;                                           #name the matrices' row and columns right
#
# Adj.tmp<- edge.to.adj(edge.test,mode = "directed")
# if(!identical(row.names(Adj.tmp),row.names(Adj))){
#   outer
# }
#
# edge<- KenNet::adj.to.edge(Adj,mode = "plus")
# Adj<- KenNet::edge.to.adj(edge,mode="plus")
#
# round(do.call(rbind,
#               lapply(1:nrow(Adj),function(r) {
#                 Adj[r,]<-Adj[r,]/sum(Adj[r,])
#               })
# ),2)
#
#
#
# for(i in ID){                                                                                                           #for each cell
#   for(j in ID){
#     IJ<- sum(Adj[i,i],Adj[j,j])                                                       #check how many observations has been made for focal i AND for focal j from the first non transposed matrix
#     Adj[i,j]<- ifelse(IJ!=0,Adj[i,j]/IJ, Adj[i,j])                                    # divide each cell by this value: taking into account observation effort of focal sampling i AND j
#   }
# }

# Adjacency matrix to association indices matrix
#
# Calculate the association indices of choice from an adjacency matrix.
#
# @param Adj Adjacency matrix of edge weights
# @param effort facultative vector of sampling effort per individual
# @param type character string. Type of effort per individual provided. Either "yab", indicates the number of times an individual have been observed associated with no-one, or "total" for its total sampling effort (cf. Whitehead, 2008). Logically, the number or association of an individual cannot be smaller than the total number of observation for this individual.
# @param index character string. Type of association index (cf. Whitehead, 2008).
# @param mode Not implemented yet. WIP
#
# @return A matrix of association index.
# @export
#
# @examples
# set.seed(42)
#
# n<- 5;nodes<- letters[1:n];
# Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
# Adj[upper.tri(Adj)]<- sample(0:8,n*(n-1)/2,replace = TRUE)
# Adj[lower.tri(Adj)]<- sample(0:8,n*(n-1)/2,replace = TRUE)
# Adj
#
# effort<- rowSums(Adj)+sample(5:8,n,replace = TRUE)
#
# adj.to.assoc(Adj,effort,type="total",index="focal")
#
# yab<- sample(5:8,n,replace = TRUE)
#
# adj.to.assoc(Adj,yab,type = "yab",index="focal")
# adj.to.assoc(Adj,yab,index="focal")
# adj.to.assoc<- function(Adj,effort=NULL,index=c("joint","sri","hw","tw","sqrt","socaff","both","focal"),
#                         type = c("total","yab"),
#                         mode = c("directed", "undirected", "max","min", "upper", "lower", "plus")){
#   if(length(mode)>1) {mode<- "plus"}
#   if(length(index)>1) {index<- "sri"}
#
#   if(nrow(Adj)!=ncol(Adj)) {stop("Adj is not a square matrix.")}
#   if(length(effort)!=nrow(Adj)) {stop("effort is not the same length as the number of individuals in Adj")}
#
#   ID<- row.names(Adj)
#
#   if(!is.null(effort)) {
#     if(length(type)>1) {
#       if(all(rowSums(Adj)<=effort)) {
#         type<- "total";
#         message("provided effort assumed to be total effort. If it is not the case, please explicit the type argument (cf. ?adj.to.assoc).")
#       }else{
#         type<- "yab";
#         message("provided effort assumed to be number of time an individual has been observed not associated (cf. ?adj.to.assoc).")
#       }
#     }
#
#     switch(type,
#            "total" = if(all(rowSums(Adj)<=effort)) diag(Adj)<- effort else stop("effort provided smaller than number of associations. This shouldn't be possible."),
#            "yab" = diag(Adj)<- effort+rowSums(Adj),
#            stop("unrecognized effort type.")
#
#     )
#   }
#
#   switch(index,
#          "sri"={
#           lapply(nrow(Adj),function(r) {
#             Adj[r,]<-Adj[r,]/sum(Adj[r,])
#           })
#          },
#          "focal"={
#            IJ<- non.zero.non.diag(Adj)
#            Adj[Adj>0&!diagonal(Adj)]<- sapply(1:nrow(IJ),
#                                               function(ij){
#                                                 i<- IJ[ij,]["row"];j<- IJ[ij,]["col"];
#                                                 Adj[i,j]/(Adj[i,i]+Adj[j,j])
#                                               })
#            diag(Adj)<-0
#            Adj
#          },
#          stop("WIP. Other methods not implemented yet. Sorry (owo)"))
# }


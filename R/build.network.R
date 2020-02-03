#' Network building function wrapper
#'
#' Wrapper for building and bootstrapping networks
#'
#' @param Obs WIP
#' @param id WIP
#' @param tar WIP
#' @param mode WIP
#' @param boot WIP
#' @param rdraw WIP
#' @param replacement WIP
#' @param directed WIP
#' @param IDS WIP
#'
#' @return WIP
#' @export
#' @importFrom igraph graph.data.frame
#'
#' @examples
#' print("WIP")
build.network<- function(Obs,id="id",tar="tar",mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),
                         boot = 10,rdraw=FALSE,replacement=TRUE, directed=FALSE,IDS=NULL){
  if(length(mode)>1) {mode<- "directed"}
  interaction.list<- ConfiBoot(Obs = Obs,id = id,tar = tar,boot = boot,replacement = replacement,proportion = 1.0,output = "for igraph")
  interaction.list
  lapply(interaction.list,function(l) igraph::graph.data.frame(l))
  # if(is.null(IDS)) IDs<- unique(factor(c(as.character(interaction.list$id),                                          #get the names of the individual present, among all individuals (error if not the same levels)
  #                                        as.character(interaction.list$tar))));
  # adj.mat.ij<- matrix(0,nrow = length(IDs),ncol = length(IDs));                                       #make matrices with the right dimensions as place-holders
  # row.names(adj.mat.ij)<- IDs;colnames(adj.mat.ij)<- IDs;                                           #name the matrices' row and columns right
  # for(l in 1:nrow(interaction.list)){
  #   adj.mat.ij[as.character(interaction.list$id[l]),as.character(interaction.list$tar[l])]<- interaction.list$N[l]
  # }
  # ifelse(directed,no = adj.mat<- adj.mat.ij+t(adj.mat.ij))
  # for(i in IDs){                                                                                                           #for each cell
  #   for(j in IDs){
  #     IJ<- sum(adj.mat.ij[i,],adj.mat.ij[j,])                                                       #check how many observations has been made for focal i AND for focal j from the first non transposed matrix
  #     adj.mat[i,j]<- ifelse(IJ!=0,adj.mat[i,j]/IJ, adj.mat[i,j])                                    # divide each cell by this value: taking into account observation effort of focal sampling i AND j
  #   }
  # }
  # SN.igraph<- graph.adjacency(adj.mat,mode = "undirected",weighted = TRUE)                          #get igraph network from adjacency matrix
  # return(SN.igraph)                                                                                 #get network
}

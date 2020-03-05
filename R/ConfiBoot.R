#' Network Bootstrap of edge-lists
#'
#' Bootstrap edge lists where repeated lines represent dyads' weights
#'
#' @param Obs WIP
#' @param id WIP
#' @param tar  WIP
#' @param boot  WIP
#' @param replacement  WIP
#' @param proportion  WIP
#' @param output  WIP
#'
#' @return  WIP
#' @export
#' @importFrom data.table data.table
#' @importFrom data.table data.table
#'
#'
#' @examples
#' print("WIP")

ConfiBoot<- function(Obs,id="id",tar="tar",boot=100,replacement=TRUE,proportion=1.0,output="list"){
  .N<-NULL;.<- NULL; #irrelevant bit of code, only to remove annoying note in R CMD Check...
  switch(class(Obs)[1],
         "data.table" = Obs,
         "data.frame" = Obs<- data.table::data.table(Obs),
         "igraph" = stop("Please provide observations as a data table containing columns id and tar (or else provide them as unquoted arguments). You can use the function adj.to.edge() for instance."),
         stop("Class of Obs not recognized. Supported classes are data.table and data.frame containing columns id and tar.")
  )
  if(is.null(Obs[[id]])|is.null(Obs[[tar]])) stop(paste0("Wrong Obs structure: either column ",id," or ",tar," is missing."))
  if(!replacement&proportion>=1) stop("If resampling is done without replacement, please provide a proportion of data in [0,1[ to be sub-sampled.")

  switch (output,
          "list" ={
            lapply(1:boot,function(b)
              if(replacement)
                Obs[sample(1:nrow(Obs),nrow(Obs),replace = TRUE)][,.N,by=c(id,tar)]
              else
                Obs[sample(1:nrow(Obs),round(nrow(Obs)*proportion,0),replace = FALSE)][,.N,by=c(id,tar)]
            )
          },
          "edge" = {
            lapply(1:boot,function(b) {
              if(replacement)
                Obs[sample(1:nrow(Obs),nrow(Obs),replace = TRUE)]
              else
                Obs[sample(1:nrow(Obs),round(nrow(Obs)*proportion,0),replace = FALSE)]
            }
            )
          },
          "for igraph" = {
            lapply(1:boot,function(b) {
              if(replacement)
                Obs<-Obs[sample(1:nrow(Obs),nrow(Obs),replace = TRUE)][,.(weight=.N),by=c(id,tar)]
              else
                Obs<-Obs[sample(1:nrow(Obs),round(nrow(Obs)*proportion,0),replace = FALSE)][,.(weight=.N),by=c(id,tar)]
              Obs<-Obs[,c(id,tar,"weight"),with=FALSE]
              colnames(Obs)<- c("from","to","weight")
              Obs
            }
            )
          }
  )
}


# edge<- ConfiNet::adj.to.edge(Adj)
# el<- ConfiBoot(edge,boot = 10,output = "edge")
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
# edge<- ConfiNet::adj.to.edge(Adj,mode = "plus")
# Adj<- ConfiNet::edge.to.adj(edge,mode="plus")
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

adj.to.assoc<- function(Adj,effort=NULL,index=c("joint","sri","hw","tw","sqrt","socaff","both"),
                        mode = c("directed", "undirected", "max","min", "upper", "lower", "plus")){
  if(length(mode)>1) {mode<- "plus"}
  if(length(index)>1) {index<- "sri"}

  ID<- row.names(Adj)
  if(!is.null(effort)) {diag(Adj)<- effort}

  switch(index,
         "sri"={
          lapply(nrow(Adj),function(r) {
            Adj[r,]<-Adj[r,]/sum(Adj[r,])
          })
         },
         stop("WIP. Other methods not implemented yet. Sorry ¯\_(ツ)_/¯"))
}

set.seed(42)

n<- 5;nodes<- letters[1:n];
Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
Adj[upper.tri(Adj)]<- sample(0:2,n*(n-1)/2,replace = TRUE)
Adj[lower.tri(Adj)]<- sample(0:1,n*(n-1)/2,replace = TRUE)
Adj

edge<- KenNet::adj.to.edge(Adj,mode = "plus")
Adj<- KenNet::edge.to.adj(edge,mode="plus")

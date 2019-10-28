KenBoot<- function(Obs,id="id",tar="tar",boot=100,replacement=TRUE,proportion=1.0,output="list"){
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


KenBoot(Obs)
KenBoot(Obs.2,"atata","ototo")

KenBoot(Obs,output = "for igraph")


build.network(Obs)
test<- build.network(Obs.2,"atata","ototo")
lapply(test,function(g) igraph::plot.igraph(g))

set.seed(42)
n<- 5;nodes<- letters[1:n];
Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
Adj[upper.tri(Adj)]<- sample(0:3,n*(n-1)/2,replace = TRUE)
Adj[lower.tri(Adj)]<- sample(0:3,n*(n-1)/2,replace = TRUE)
Adj
Obs<- adj.to.edge(Adj,mode = "directed")

Obs.2<- Obs
colnames(Obs.2)<- c("atata","ototo")

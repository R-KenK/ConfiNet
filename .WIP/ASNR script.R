asnr.weighted<- list.files("C:/R/Git/asnr/Networks/Mammalia/",pattern = "_weighted",full.names = TRUE)

asnr.list<- lapply(asnr.weighted,
                   function(path){
                     list.files(path,pattern = ".graphml",full.names = TRUE)
                   }
)
asnr.weighted[which(sapply(asnr.list,length)>1)]

asnr.Adj<- lapply(asnr.list[which(sapply(asnr.list,length)==1)],
                  function(path){
                    Adj<- import_from_graphml(path,"adjacency")
                    attr(Adj,"path")<- path
                    Adj
                  }
)

asnr.graph<- lapply(asnr.list[which(sapply(asnr.list,length)==1)],
                  function(path){
                    G<- import_from_graphml(path,"graph")
                    attr(G,"path")<- path
                    G
                  }
)


asnr.Adj
asnr.graph
sapply(asnr.Adj,max)
sapply(asnr.Adj,sum)/2


plot(import_from_graphml(asnr.list[which(sapply(asnr.list,length)==1)][[10]],"graph"),sparse = F)

igraph::edge.attributes(import_from_graphml(asnr.list[which(sapply(asnr.list,length)==1)][[10]],"graph"))

default_explorative_net.plot(asnr.graph[[2]],edge.with.mul = 9*1/max(asnr.Adj[[2]]),vertex.size.mul = 3*1/max(asnr.Adj[[2]]),centrality.fun = "strength")

default_explorative_net.plot(asnr.graph[[8]],2*1/max(asnr.Adj[[8]]),5*1/max(asnr.Adj[[8]]),"strength")

lapply(seq_along(asnr.graph),
       function(G){
         default_explorative_net.plot(asnr.graph[[G]],3*1/max(asnr.Adj[[G]]),5*1/max(asnr.Adj[[G]]),"strength")
       }
)

lapply(seq_along(asnr.graph),
       function(G){
         hist(igraph::strength(graph = asnr.graph[[G]]))
       }
)


igraph::vertex.attributes(asnr.graph[[8]])

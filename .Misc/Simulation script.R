set.seed(42)
G<- igraph::read_graph("C:/R/Git/asnr/Networks/Mammalia/rhesusmacaque_association_weighted/weighted_Contact_Sits_Macaque.graphml",format = "graphml")

total_scan<- 1138; #from original paper

Adj<- as.matrix(igraph::as_adj(G,attr = "weight"))
Adj<- as.matrix(igraph::as_adj(G,attr = "weight"))
row.names(Adj)<- as.character(1:nrow(Adj));colnames(Adj)<- row.names(Adj)
Adj

focal.list<- sample(row.names(Adj),1138,replace = TRUE)
table(focal.list)

test<- ConfiNet::Boot_scans(Adj = Adj,n.boot = 100,total_scan = total_scan,
                     method = "both",focal.list = focal.list,scaled = TRUE,mode = "max",output = "adj",n.cores = 7)
head(test)





# progress bar testing ----------------------------------------------------

source("R/matrix.tools.R")
source("R/Binary.prob.R")
source("R/do.scan.R")
source("R/sum_up.scans.R")
source("R/iterate_scans.R")
source("R/Boot_scans.R")
source("R/observable_edges.R")
source("R/Bootstrap_tools.R")

G<- igraph::read_graph("C:/R/Git/asnr/Networks/Mammalia/rhesusmacaque_association_weighted/weighted_Contact_Sits_Macaque.graphml",format = "graphml")
total_scan<- 1138;n.boot<- 100;
Adj<- as.matrix(igraph::as_adj(G,attr = "weight"))
row.names(Adj)<- as.character(1:nrow(Adj));colnames(Adj)<- row.names(Adj)
do.scan(Adj=Adj,total_scan=total_scan,mode = "max")

focal.list<- sample(row.names(Adj),total_scan,replace = TRUE)
table(focal.list)

Bootstrap<- Boot_scans(Adj = Adj,n.boot = n.boot,total_scan = total_scan,obs.prob = 0.6,keep = TRUE,
                  method = "both",focal.list = focal.list,scaled = TRUE,mode = "directed",output = "all",n.cores = 7)

sapply(1:n.boot,
       function(b) {
         cor(c(Boot_get.list(test,"theoretical")$adjacency[[b]]),   # c() flattens the matrix to consider it like a vector
             c(Boot_get.list(test,"observed")$adjacency[[b]]))
       }
)

sapply(1:n.boot,
       function(b) {
         cor(c(Boot_get.list(test,"theoretical")$adjacency[[b]]),
             c(Boot_get.list(test,"focal")$adjacency[[b]]))
       }
)

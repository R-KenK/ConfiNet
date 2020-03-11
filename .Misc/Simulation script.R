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



# Simulation script backbone ----------------------------------------------

## import needed function, but ultimately should use library(ConfiNet) --------
source("R/matrix.tools.R")
source("R/Binary.prob.R")
source("R/do.scan.R")
source("R/sum_up.scans.R")
source("R/iterate_scans.R")
source("R/Boot_scans.R")
source("R/observable_edges.R")
source("R/Bootstrap_tools.R")


# import and generate objects ---------------------------------------------
set.seed(42)

n<- 10;nodes<- letters[1:n];
total_scan<- 100; #from original paper
n.boot<- 50;

Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
Adj[non.diagonal(Adj)]<- sample((0:round(total_scan*.50)),n*(n-1),replace = TRUE)
Adj

focal.list<- sample(nodes,total_scan,replace = TRUE)
table(focal.list)


# Parameter choices -------------------------------------------------------
OBS.PROB<- seq(0.1,0.9,by = 0.2)
MODE<- c("directed","max")
FOCAL.LIST<- list(random = sample(nodes,total_scan,replace = TRUE),
                  even = rep_len(nodes,length.out = total_scan),
                  biased = "TO IMPLEMENT")

# Boot_scans() wrapper to loop through parameters -------------------------
Boot_with.parameter.list<- function(Adj,n.obs,total_scan,obs.prob,mode,focal.list){
  Boot_scans(Adj = Adj,n.boot = n.boot,total_scan = total_scan,obs.prob = 0.6,keep = TRUE,
             method = "both",focal.list = focal.list,scaled = TRUE,mode = "directed",output = "all",n.cores = 7)
}
start<- Sys.time()
Bootstrap.list<- lapply(seq_along(OBS.PROB),
                        function(param.obs){
                          lapply(seq_along(FOCAL.LIST[1:2]),
                                 function(param.foc){
                                   lapply(seq_along(MODE),
                                          function(param.mode){
                                            Boot_with.parameter.list(Adj,n.obs,total_scan,
                                                                     obs.prob = OBS.PROB[[param.obs]],
                                                                     mode = MODE[[param.mode]],
                                                                     focal.list = FOCAL.LIST )
                                          }
                                   )
                                 }
                          )
                        }
)
stop<- Sys.time()

# sapply(1:n.boot,
#        function(b) {
#          cor(c(Boot_get.list(test,"theoretical")$adjacency[[b]]),   # c() flattens the matrix to consider it like a vector
#              c(Boot_get.list(test,"observed")$adjacency[[b]]))
#        }
# )
#
# sapply(1:n.boot,
#        function(b) {
#          cor(c(Boot_get.list(test,"theoretical")$adjacency[[b]]),
#              c(Boot_get.list(test,"focal")$adjacency[[b]]))
#        }
# )

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


# Quick plot for presentation ---------------------------------------------

# Simulation script backbone ----------------------------------------------

## import needed function, but ultimately should use library(ConfiNet) --------
source("R/matrix.tools.R")
source('R/rare.event.optimization_tools.R')
source("R/Misc_tools.R")
source("R/Binary.prob.R")
source("R/do.scan.R")
source("R/sum_up.scans.R")
source("R/iterate_scans.R")
source("R/Boot_scans.R")
source("R/obs.prob_tools.R")
source("R/focal.list.R")
source("R/Bootstrap_tools.R")
source(".WIP/ASNR.tools.R")
load("C:/R/Git/ConfiNet/R/sysdata.rda")

# import and generate objects ---------------------------------------------
# Here preferably should be implemented as automatic import from ASNR/networkdata

set.seed(42)
n.boot<- 10;

asnr.weighted.dir<- list.files("C:/R/Git/asnr/Networks/Mammalia/",pattern = "_weighted",full.names = TRUE)

asnr.list<- lapply(asnr.weighted.dir,
                   function(path){
                     list.files(path,pattern = ".graphml",full.names = TRUE)
                   }
)

asnr.Adj<- lapply(asnr.list[which(sapply(asnr.list,length)==1)],
                  function(path){
                    Adj<- import_from_graphml(path,"adjacency")
                    attr(Adj,"path")<- path
                    Adj
                  }
)
ADJ<- asnr.Adj
sapply(ADJ,function(adj) attr(adj,"path"))

TOTAL_SCAN<- lapply(asnr.weighted.dir[which(sapply(asnr.list,length)==1)],get.total_scan)

with.total_scan<- !sapply(TOTAL_SCAN,is.null)
ADJ<- ADJ[with.total_scan]
TOTAL_SCAN<- TOTAL_SCAN[with.total_scan]

ADJ<- lapply(seq_along(ADJ),
             function(a) {
               Adj<- ADJ[[a]]
               n<- nrow(Adj)
               attr(Adj,"total_scan")<- TOTAL_SCAN[[a]]
               attr(Adj,"use.rare.opti")<- decide_use.rare.opti(n = n,total_scan = total_scan,max.obs = max(Adj),alpha = 0.05)
               Adj
             }
)


with.n.inf100<- sapply(ADJ,nrow)<100
ADJ<- ADJ[with.n.inf100]
TOTAL_SCAN<- TOTAL_SCAN[with.n.inf100]


path.list<- lapply(ADJ,function(Adj) attr(Adj,"path"))
TOTAL_SCAN<- sapply(ADJ,function(Adj) attr(Adj,"total_scan"))
net.list<- lapply(path.list,function(path) import_from_graphml(path,"graph"))

species.name<- lapply(path.list,
                      function(path) {
                        gsub("[_].*$",
                             "",
                             gsub("^C:/R/Git/asnr/Networks/Mammalia/",
                                  "",
                                  path)
                        )
                      }
)


library(ggplot2)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Make a color palette
colour.palette<- gg_color_hue(length(ADJ))

# Manual loop lol
net<- 1

graph<- net.list[[net]]
name<- species.name[[net]]
Adj<- ADJ[[net]]
colour<- colour.palette[[net]]
png(filename = paste0(".WIP/SN plots/",name,".png"),width = 10,height = 10,units = "in",res = 300)
default_explorative_net.plot(graph,edge.with.mul = 7.5*1/max(Adj),vertex.size.mul = 5*1/max(Adj),centrality.fun = "strength",main = name, vertex.color=colour)
dev.off()
net<- net+1


# Cor plot theoretical/original -------------------------------------------
OG.theo<- lapply(seq_along(ADJ),
                   function(a){
                     cat(paste0("Network ",a,"/",length(ADJ)," @ ",Sys.time(),"\n"))
                     Adj<- ADJ[[a]]
                     total_scan<- attr(Adj,"total_scan")
                     iterate_scans(Adj = Adj,total_scan = total_scan,obs.prob = 1,use.rare.opti = FALSE,
                                   method = "theoretical",scaled = FALSE,mode = "max",output = "adjacency")
                   }
)

library(data.table)
OG.theo.dt<- rbind_lapply(seq_along(ADJ),
                          function(a){
                            cat(paste0("Network ",a,"/",length(ADJ)," @ ",Sys.time(),"\n"))
                            Adj<- non.diagonal(ADJ[[a]],"vec")
                            Theo<- non.diagonal(OG.theo[[a]]$theoretical,"vec")
                            data.table(Network = as.character(a),OG=Adj,Theoretical=Theo)
                          }
)
plot<- ggplot(OG.theo.dt,aes(OG,Theoretical,colour=Network))+
  facet_wrap(.~Network,scales = "free",ncol = 4)+
  geom_abline(slope = 1,0,colour="black",lty="dashed")+
  geom_smooth(method = "lm",colour="darkred",alpha=0.3,se=FALSE)+
  geom_point(alpha=0.5)+
  guides(colour=FALSE)+theme_bw()
ggsave(".WIP/SN plots/OG.Theo.png",plot,width = 16,height = 10,units = "in",dpi = 200)

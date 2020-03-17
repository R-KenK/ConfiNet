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
# Here preferably should be impleemnted as automatic import from ASNR

set.seed(42)

n<- 15;nodes<- as.character(1:n);
n.boot<- 100;

Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
Adj[non.diagonal(Adj)]<- sample((0:round(total_scan*.50)),n*(n-1),replace = TRUE)
Adj

# Parameter choices -------------------------------------------------------
# for each variable type, there should be only a non-nested list of parameters. I'll figure out later how to group similar values through factor ifelse() and substring I guess...

OBS.PROB<- list(trait = obs.prob_bias(Adj = Adj,obs.prob_fun = prod,bias_fun = NULL,reverse = FALSE),
                network = obs.prob_bias(Adj = Adj,obs.prob_fun = prod,
                                        bias_fun = function(node) igraph::strength(igraph::graph.adjacency(Adj,weighted = TRUE))[node],
                                        reverse = FALSE)
)
OBS.PROB<- c({unb<- seq(0.1,0.9,by = 0.2);names(unb)<- paste0("unbiased_",unb);as.list(unb)},OBS.PROB) # c() over two lists makes them flat while allowing for shorter calls
MODE<- as.list(c(directed = "directed",max = "max",min = "min",plus = "plus"))
FOCAL.LIST<- list(random = sample(nodes,total_scan,replace = TRUE),
                  even = rep_len(nodes,length.out = total_scan),
                  biased = "TO IMPLEMENT")

parameters.comb<- expand.grid(
  list(mode = 1:length(MODE),
       focal.list = 1:length(FOCAL.LIST[1:2]),
       obs.prob = 1:length(OBS.PROB)
  )
)

parameters.list<- lapply(1:nrow(parameters.comb),
                         function(p){
                           list(
                             obs.prob = {
                               obs.prob<- OBS.PROB[[parameters.comb[p,"obs.prob"]]];
                               attr(obs.prob,"name")<- names(OBS.PROB)[parameters.comb[p,"obs.prob"]];
                               obs.prob
                             },
                             focal.list = {
                               focal.list<- FOCAL.LIST[[parameters.comb[p,"focal.list"]]];
                               attr(focal.list,"name")<- names(FOCAL.LIST)[parameters.comb[p,"focal.list"]];
                               focal.list
                             },
                             mode = {
                               mode<- MODE[[parameters.comb[p,"mode"]]];
                               attr(mode,"name")<- names(MODE)[parameters.comb[p,"mode"]];
                               mode
                             }
                           )
                         }
)

# Iterated Boot_scans() through parameters.list -------------------------
start<- Sys.time()
Bootstrap.list<- lapply(seq_along(parameters.list),
                        function(p){
                          obs.prob<- parameters.list[[p]]$obs.prob;
                          mode<- parameters.list[[p]]$mode;
                          focal.list<- parameters.list[[p]]$focal.list
                          boot_progress.param(p)
                          Boot_scans(Adj = Adj,n.boot = n.boot,total_scan = total_scan,obs.prob = obs.prob,keep = TRUE,
                                     method = "both",focal.list = focal.list,scaled = TRUE,mode = mode,output = "all",n.cores = 7)
                        }
)
stop<- Sys.time()
stop-start

data.long<- rbind_lapply(seq_along(Bootstrap.list),
                         function(B){
                           rbind(data.frame(boot = B,method = "group",
                                            Boot_get.param(parameters.list[[B]]),
                                            cor = adjacency_cor(Bootstrap = Bootstrap.list[[B]],what = "observed")),
                                 data.frame(boot = B,method = "focal",
                                            Boot_get.param(parameters.list[[B]]),
                                            cor = adjacency_cor(Bootstrap = Bootstrap.list[[B]],what = "focal")))
                         }
)

library(ggplot2)
library(data.table)
data.long<- data.table(data.long)
data.summary<- data.long[,.(cor=median(cor),sd=sd(cor)),by = .(obs.prob,focal.list,mode,method)]

ggplot(data.long[obs.prob %in% c("network","trait")],aes(interaction(method,obs.prob),cor,colour = method))+facet_grid(mode~focal.list)+geom_jitter(alpha=0.2)+geom_boxplot()+theme_bw()
ggplot(data.long[!(obs.prob %in% c("network","trait"))],aes(interaction(method,obs.prob),cor,colour = method))+facet_grid(mode~focal.list)+geom_jitter(alpha=0.2)+geom_boxplot()+theme_bw()

ggplot(data.summary,aes(interaction(method,obs.prob),cor,fill = method))+facet_grid(mode~focal.list)+
  geom_errorbar(aes(ymin = cor-sd,ymax=cor+sd),width = 0.2)+geom_bar(stat = "identity",alpha=0.8)+theme_bw()


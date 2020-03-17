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

n<- 10;nodes<- as.character(1:n);
total_scan<- 100; #from original paper
n.boot<- 5;

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
attr(parameters.list[[1]]$obs.prob,"name")
# Boot_scans() wrapper to loop through parameters -------------------------
#' Title
#'
#' @param Adj
#' @param n.obs
#' @param total_scan
#' @param obs.prob
#' @param mode
#' @param focal.list
#'
#' @return
#' @export
#'
#' @examples
boot_progress.param<- function(obs.prob,mode,focal.list){
  cat(paste0("obs.prob = ",attr(obs.prob,"name")," - focal.list = ",attr(focal.list,"name")," - mode = ",attr(mode,"name"),"\n"))
}

#' Title
#'
#' @param Bootstrap
#' @param n.boot
#' @param what
#'
#' @return
#' @export
#'
#' @examples
adjacency_cor<- function(Bootstrap,what = c("observed","focal"),n.boot = length(Bootstrap)){
  what<- match.arg(what)
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           cor(c(Boot_get.list(Bootstrap,"theoretical")$adjacency[[b]]),   # c() flattens the matrix to consider it like a vector
               c(Boot_get.list(Bootstrap,what)$adjacency[[b]]))
         }
  )
}

Bootstrap.list<- lapply(parameters.list,
                        function(p){
                          obs.prob<- p$obs.prob;mode<- p$mode;focal.list<- p$focal.list
                          boot_progress.param(obs.prob = obs.prob,mode = mode,focal.list = focal.list)
                          Boot_scans(Adj = Adj,n.boot = n.boot,total_scan = total_scan,obs.prob = obs.prob,keep = TRUE,
                                     method = "both",focal.list = focal.list,scaled = scaled,mode = mode,output = "all",n.cores = 7)
                        }
)

data<- rbind_lapply(Bootstrap.list,
                    function(param.obs){
                      rbind_lapply(param.obs,
                                   function(param.foc){
                                     rbind_lapply(param.foc,
                                                  function(L){
                                                    L$"stats"
                                                  }
                                     )
                                   }
                      )
                    }
)

library(ggplot2)
library(data.table)
data.long<- data.table(reshape2::melt(data,id.vars = 1:3,measure.vars = 4:5,variable.name = "method",value.name = "cor"))
data.summary<- data.long[,.(cor=median(cor),sd=sd(cor)),by = .(obs.prob,focal.list,mode,method)]
ggplot(data.long,aes(interaction(method,obs.prob),cor,colour = method))+facet_grid(mode~focal.list)+geom_jitter(alpha=0.2)+geom_boxplot()+theme_bw()
ggplot(data.summary,aes(interaction(method,obs.prob),cor,fill = method))+facet_grid(mode~focal.list)+
  geom_errorbar(aes(ymin = cor-sd,ymax=cor+sd),width = 0.2)+geom_bar(stat = "identity")+theme_bw()

str(Bootstrap.list)

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

n<- 10;nodes<- letters[1:n];
total_scan<- 100; #from original paper
n.boot<- 5;

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
Boot_with.parameter.list<- function(Adj,n.obs,total_scan,obs.prob,mode,focal.list){
  cat(paste0("obs.prob = ",obs.prob," - focal.list = ",names(focal.list)," - mode = ",mode,"\n"))
  Boot_scans(Adj = Adj,n.boot = n.boot,total_scan = total_scan,obs.prob = 0.6,keep = TRUE,
             method = "both",focal.list = focal.list[[1]],scaled = TRUE,mode = "directed",output = "all",n.cores = 7)
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

#' Title
#'
#' @param X
#' @param FUN
#'
#' @return
#' @export
#'
#' @examples
rbind_lapply<- function(X,FUN){
  do.call(rbind,lapply(X = X,FUN = FUN))
}

start<- Sys.time()
Bootstrap.list<- lapply(seq_along(OBS.PROB),
                        function(param.obs){
                          lapply(seq_along(FOCAL.LIST[1:2]),
                                 function(param.foc){
                                   lapply(seq_along(MODE),
                                          function(param.mode){
                                            Bootstrap<- Boot_with.parameter.list(Adj,n.obs,total_scan,
                                                                                 obs.prob = OBS.PROB[[param.obs]],
                                                                                 mode = MODE[[param.mode]],
                                                                                 focal.list = FOCAL.LIST[param.foc])
                                            list(Bootstrap = Bootstrap,
                                                 stats = data.table::data.table(obs.prob = OBS.PROB[[param.obs]],
                                                                                focal.list = as.factor(names(FOCAL.LIST[param.foc])),
                                                                                mode = as.factor(MODE[[param.mode]]),
                                                                                group.cor = adjacency_cor(Bootstrap,what = "observed"),
                                                                                focal.cor = adjacency_cor(Bootstrap,what = "focal")
                                                 )
                                            )
                                          }
                                   )
                                 }
                          )
                        }
)
stop<- Sys.time()
stop-start

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

# Analysis script backbone ----------------------------------------------

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
n.boot<- 50;

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

Adj<- ADJ[[4]]
total_scan<- TOTAL_SCAN[[4]]

# Parameter choices -------------------------------------------------------


# Listing desired unbiased obs.prob ####
unb.fun_list<- list(identity = identity)
unb.subtype_list<- lapply(seq(0.1,0.9,by = 0.2),
                          function(x){
                            bias.subtype<- function(i,j,Adj){x}
                            bias.subtype
                          }
)
names(unb.subtype_list)<- seq(0.1,0.9,by = 0.2)

# Listing desired trait-based function components ####
trait.fun_list<- list(plus = `+` #,
                      # prod = `*`
)

# Assembling all obs.prob biases into one list of lists ####
obs.bias_list<- list(
  unbias = list(bias.fun_list = unb.fun_list,
                bias.subtype_list = unb.subtype_list,
                type = "unbiased"
  )
)

focal.bias_list<- unlist(c(list(even = "even")),recursive = FALSE)
# focal.bias_list<- unlist(c(list(even = "even"),list(trait = trait.focal.fun_list),list(net = net.focal.fun_list)),recursive = FALSE)

# Generate parameters list for each network once and for all --------------
PARAMETERS.LIST<- initialize_parameter_list(ADJ[4],obs.bias_list,focal.bias_list)
parameters.list<- PARAMETERS.LIST[[1]]
total_scan<- attr(Adj,"total_scan")
use.rare.opti<- attr(Adj,"use.rare.opti")

sort( sapply(ls(),function(x){object.size(get(x))}))
rm("opti.model","standard.model","asnr.Adj")

cl<- snow::makeCluster(5);snow::clusterExport(cl,list = ls());
Bootstrap.list.params<- pbapply::pblapply(seq_along(parameters.list),
                                   function(p){
                                     obs.prob<- parameters.list[[p]]$obs.prob;
                                     mode<- "max"
                                     focal.list<- parameters.list[[p]]$focal.list
                                     Boot_scans(Adj = Adj,n.boot = n.boot,total_scan = total_scan,obs.prob = obs.prob,use.rare.opti = use.rare.opti,
                                                method = "both",focal.list = focal.list,scaled = FALSE,mode = mode,output = "all")
                                   },cl = cl
)
snow::stopCluster(cl)

library(data.table)
method.aij.dt<- function(Bootstrap.list, method){
  rbind_lapply(1:n.boot,
               function(b){
                 mat<- Boot_get.list(Bootstrap.list,method,"adj")[[b]]
                 ij<- data.table(expand.grid(i = 1:nrow(Adj),j = 1:nrow(Adj)))[i<j]
                 n.obs<- attr(mat,"observed_edges")
                 data.table(boot = b,
                            method = method,
                            obs.prob = switch(method,
                                              "theoretical" = 1,
                                              "group" = as.numeric(attr(attr(Bootstrap.list,"attr.list")[["obs.prob"]],"subtype")),
                                              "focal" = 1/nrow(Adj)
                            ),
                            ij,
                            aij = mat[upper.tri(mat)],
                            effort = n.obs[upper.tri(n.obs)])
               }
  )
}

aij.dt<- rbind_lapply(seq_along(Bootstrap.list.params),
                      function(p){
                        Bootstrap.list<- Bootstrap.list.params[[p]]

                        presence.prob<- attr(Bootstrap.list,"attr.list")[["presence.prob"]]
                        obs.prob<- as.numeric(attr(attr(Bootstrap.list,"attr.list")[["obs.prob"]],"subtype"))

                        data.table(param = as.character(p),
                                   presence.prob = rep(presence.prob[upper.tri(presence.prob)],3*n.boot),
                                   rbind(
                                     method.aij.dt(Bootstrap.list,"theoretical"),
                                     method.aij.dt(Bootstrap.list,"group"),
                                     method.aij.dt(Bootstrap.list,"focal")
                                   )
                        )
                      }
)

library(ggplot2)
library(ggstatsplot)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

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
n.boot<- 5;

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

# with.total_scan.inf1000<- sapply(TOTAL_SCAN,function(t) t<1000)
# ADJ<- ADJ[with.total_scan.inf1000]
# TOTAL_SCAN<- TOTAL_SCAN[with.total_scan.inf1000]


# Special treatment of the per week racoon networks -----------------------
# asnr.racoon.path<- list.files("C:/R/Git/asnr/Networks/Mammalia/raccoon_proximity_weighted/",pattern = ".graphml",full.names = TRUE)
# asnr.racoon.Adj<- lapply(asnr.racoon.path,
#                   function(path){
#                     Adj<- import_from_graphml(path,"adjacency")
#                     attr(Adj,"path")<- path
#                     Adj
#                   }
# )
# ADJ.racoon<- asnr.racoon.Adj
# sapply(ADJ.racoon,function(adj) attr(adj,"path"))
#
# TOTAL_SCAN.racoon<- lapply("C:/R/Git/asnr/Networks/Mammalia/raccoon_proximity_weighted/",get.total_scan)
#
# PARAMETERS.LIST.racoon<- lapply(seq_along(ADJ.racoon),
#                                 function(a){
#                                   cat(paste0(a,"/",length(ADJ.racoon)," @ ",Sys.time(),"\n"))
#                                   Adj<- ADJ.racoon[[a]]
#                                   total_scan<- TOTAL_SCAN.racoon[[1]]
#                                   initialize_parameters(Adj,total_scan)
#                                 }
# )

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
trait.subtype_list<- list(identity = function(i,Adj) {i},
                          # double = function(i,Adj){i*2},
                          # normal = function(i,Adj){n<- nrow(Adj);sort(dnorm(1:n,n/2,n/5))[i]},
                          # square = function(i,Adj){i^2},
                          power.4 = function(i,Adj){i^4},
                          # sqrt = function(i,Adj){sqrt(i)},
                          tan = function(i,Adj){tan((i-nrow(Adj)/2)/nrow(Adj))^3},
                          # sin = function(i,Adj){sin((i-nrow(Adj)/2)/nrow(Adj))^3},
                          # tanh = function(i,Adj){tanh(i-nrow(Adj)/2)},
                          sigmoid = function(i,Adj){exp(i-nrow(Adj)/2)/(1+exp(i-nrow(Adj)/2))},
                          log = function(i,Adj){log(i-.99)},
                          exp = function(i,Adj){exp(i)}
)

# Listing desired net-based function components ####
net.fun_list<- list(plus = `+`#,
                    # prod = `*`
)
net.subtype_list<- list(EV = function(Adj){compute.EV(graph = Adj,"max")},
                        # degree = function(Adj){compute.deg(graph = Adj,"max")},
                        strength = function(Adj){compute.strength(graph = Adj,"max")}
)
# focal.list biases

trait.focal.fun_list<- list(identity = function(n,Adj) {1:n},
                            # double = function(n,Adj){1:n*2},
                            # normal = function(n,Adj){sort(dnorm(1:n,n/2,n/5))},
                            # square = function(n,Adj){(1:n)^2},
                            power.4 = function(n,Adj){(1:n)^4},
                            # sqrt = function(n,Adj){sqrt(1:n)},
                            tan = function(n,Adj){tan(((1:n)-n/2)/n)^3-tan((1-n/2)/n)^3},
                            # sin = function(n,Adj){sin(((1:n)-n/2)/n)^3-sin((1-n/2)/n)^3},
                            # tanh = function(n,Adj){tanh((1:n)-n/2)-tanh(1-n/2)},
                            sigmoid = function(n,Adj){exp((1:n)-n/2)/(1+exp((1:n)-n/2))},
                            log = function(n,Adj){log((1:n)-.99)-log(0.01)},
                            exp = function(n,Adj){exp(1:n)}
)

net.focal.fun_list<- list(EV = function(n,Adj){compute.EV(graph = Adj,"max")},
                          degree = function(n,Adj){compute.deg(graph = Adj,"max")},
                          strength = function(n,Adj){compute.strength(graph = Adj,"max")}
)
# Assembling all obs.prob biases into one list of lists ####
obs.bias_list<- list(
  unbias = list(bias.fun_list = unb.fun_list,
                bias.subtype_list = unb.subtype_list,
                type = "unbiased"
  ),
  trait = list(bias.fun_list = trait.fun_list,
               bias.subtype_list = trait.subtype_list,
               type = "trait"
  ),
  net = list(bias.fun_list = net.fun_list,
             bias.subtype_list = net.subtype_list,
             type = "net"
  )
)

focal.bias_list<- unlist(c(list(even = "even"),list(trait = trait.focal.fun_list),list(net = net.focal.fun_list)),recursive = FALSE)

# Generate parameters list for each network once and for all --------------
PARAMETERS.LIST<- initialize_parameter_list(ADJ,obs.bias_list,focal.bias_list)

# Iteration through networks, bootstrap and gather data in data frame -------------------------
cl<- snow::makeCluster(7);snow::clusterExport(cl,list = ls());
start_pblapply<- Sys.time()
data.long<- rbind_lapply(seq_along(ADJ),
                   function(a){
                     cat(paste0("Network ",a,"/",length(ADJ)," @ ",Sys.time(),"\n"))
                     Adj<- ADJ[[a]]
                     total_scan<- attr(Adj,"total_scan")
                     use.rare.opti<- attr(Adj,"use.rare.opti")
                     parameters.list<- PARAMETERS.LIST[[a]]
                     Bootstrap.list<- pbapply::pblapply(seq_along(parameters.list),
                                             function(p){
                                               obs.prob<- parameters.list[[p]]$obs.prob;
                                               mode<- "max"
                                               # mode<- parameters.list[[p]]$mode;
                                               focal.list<- parameters.list[[p]]$focal.list
                                               Boot_scans(Adj = Adj,n.boot = n.boot,total_scan = total_scan,obs.prob = obs.prob,use.rare.opti = use.rare.opti,
                                                          method = "both",focal.list = focal.list,scaled = FALSE,mode = mode,output = "adjacency")
                                             },cl = cl
                     )
                     cat(paste0("\nCompiling data @ ",Sys.time(),"\n"))
                     Get.data(a,Bootstrap.list)
                   }
)
stop_pblapply<- Sys.time()
snow::stopCluster(cl)

stop_pblapply-start_pblapply
# data.long.10.boot<- data.long
# saveRDS(data.long.10.boot,file = ".WIP/data.long.10.boot.rds")
# data.long.100.boot<- data.long
# saveRDS(data.long.100.boot,file = ".WIP/data.long.100.boot.rds")
# data.long<- readRDS(file = ".WIP/data.long.10.boot.rds")
# Draft of data handling, plotting and analysis ---------------------------
library(data.table)
data.long<- data.table(data.long)
library(ggplot2)
mytheme<- theme_bw()+theme(plot.title = element_text(lineheight=.9, face="bold"),axis.line = element_line(colour = "black"),panel.grid.major.x = element_line(colour = "grey95",size = .2),panel.grid.minor.x = element_line(colour = "grey95",linetype = "dotted"),panel.grid.major.y =element_line(colour = "grey95",size = .2),panel.grid.minor.y=element_blank())+theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.75))

data.summary<- data.long[,.(cor=median(cor),sd.cor=sd(cor),
                            # degree=median(degree),sd.degree=sd(degree),
                            strength=median(strength),sd.strength=sd(strength),
                            EV=median(EV),sd.EV=sd(EV),
                            CC=median(CC),sd.CC=sd(CC),
                            Frob=median(Frob),sd.Frob=sd(Frob),
                            Frob.GOF=median(Frob.GOF),sd.Frob.GOF=sd(Frob.GOF),
                            SLap=median(SLap),sd.SLap=sd(SLap),
                            SLap.GOF=median(SLap.GOF),sd.SLap.GOF=sd(SLap.GOF),
                            obs.cor=median(obs.cor),sd.obs.cor=sd(obs.cor)),by = .(Network,obs.prob.type,obs.prob.subtype,focal.list.type,focal.list.subtype,mode,method,scaled)]
summary(data.summary)

data.summary[,para:=paste(obs.prob.type,obs.prob.subtype,focal.list.type,focal.list.subtype,sep = "-"),]

unique(data.summary$para)

# data.summary<- data.summary[scaled==TRUE]

# Exploratory graphs ------------------------------------------------------
data.summary$Network<- factor(data.summary$Network,levels = 1:12)
data.summary$Network.txt<- factor(paste0("Network ",data.summary$Network),levels = paste0("Network ",1:12))

vars<- c("cor","strength","Frob.GOF","SLap.GOF","obs.cor")

ggplot_points.lines.and.boxplot<- function(var){
  plot<- ggplot(data.summary,aes_string("method",var,colour="method",fill="method",group="para"))+
    facet_wrap(.~Network.txt+scaled,ncol=8)+
    geom_hline(yintercept = 0,colour='grey50')+geom_hline(yintercept = 1,colour='grey75',lty='dashed')+
    geom_line(alpha=0.3,colour='grey80')+geom_point(alpha=0.1)+
    geom_violin(aes(group=method),alpha=0.3,colour='grey50')+theme_bw()
  plot
}

ggplot_points.lines.and.boxplot("obs.cor")

for(var in vars){
  (plot<- ggplot_points.lines.and.boxplot(var))
  ggsave(paste0(".WIP/SN plots/Raw ",var," boxplot.png"),plot = plot,width = 10,height = 6,units = "in",dpi = 200)
}

# Draft PCA ---------------------------------------------------------------

data.long.scaled<- data.long;
data.long.scaled[,c("cor","degree","strength","EV","ClustCoef","Frob","Frob.GOF","SLap","SLap.GOF")]<- data.table(scale(data.long.scaled[,c("cor","degree","strength","EV","ClustCoef","Frob","Frob.GOF","SLap","SLap.GOF")]))
PCA<- FactoMineR::PCA(data.long.scaled[Network=="1",c("cor","degree","strength","EV","ClustCoef","Frob","Frob.GOF","SLap","SLap.GOF")])
PCA$eig
PCA$var$contrib
PCA$var$coord

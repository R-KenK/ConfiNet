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
                          double = function(i,Adj){i*2},
                          square = function(i,Adj){i^2},
                          power.4 = function(i,Adj){i^4},
                          sqrt = function(i,Adj){sqrt(i)},
                          tan = function(i,Adj){tan((i-nrow(Adj)/2)/nrow(Adj))^3},
                          sin = function(i,Adj){sin((i-nrow(Adj)/2)/nrow(Adj))^3},
                          tanh = function(i,Adj){tanh(i-nrow(Adj)/2)},
                          sigmoid = function(i,Adj){exp(i-nrow(Adj)/2)/(1+exp(i-nrow(Adj)/2))},
                          log = function(i,Adj){log(i-.99)},
                          exp = function(i,Adj){exp(i)}
)

# Listing desired net-based function components ####
net.fun_list<- list(plus = `+`,
                    prod = `*`
)
net.subtype_list<- list(EV = function(Adj){compute.EV(graph = Adj,"max")},
                        strength = function(Adj){compute.strength(graph = Adj,"max")},
                        degree = function(Adj){compute.deg(graph = Adj,"max")}
)

# Assembling all into one list of lists ####
global_list<- list(
  list(bias.fun_list = unb.fun_list,
       bias.subtype_list = unb.subtype_list,
       type = "unbiased"
  ),
  list(bias.fun_list = trait.fun_list,
       bias.subtype_list = trait.subtype_list,
       type = "trait"
  ),
  list(bias.fun_list = net.fun_list,
       bias.subtype_list = net.subtype_list,
       type = "net"
  )
)

# Assembling all into one list of lists ####
initialize_obs.prob_list(ADJ,global_list)

# Generate parameters list for each network once and for all --------------
PARAMETERS.LIST<- lapply(seq_along(ADJ),
                         function(a){
                           cat(paste0(a,"/",length(ADJ)," @ ",Sys.time(),"\n"))
                           Adj<- ADJ[[a]]
                           total_scan<- attr(Adj,"total_scan")
                           initialize_parameters(Adj,total_scan)
                         }
)

# Iteration through networks, bootstrap and gather data in data frame -------------------------
data.long<- rbind_lapply(seq_along(ADJ),
                   function(a){
                     cat(paste0(a,"/",length(ADJ)," @ ",Sys.time(),"\n"))
                     Adj<- ADJ[[a]]
                     total_scan<- attr(Adj,"total_scan")
                     use.rare.opti<- attr(Adj,"use.rare.opti")
                     parameters.list<- PARAMETERS.LIST[[a]]
                     Bootstrap.list<- lapply(seq_along(parameters.list),
                                             function(p){
                                               obs.prob<- parameters.list[[p]]$obs.prob;
                                               mode<- parameters.list[[p]]$mode;
                                               focal.list<- parameters.list[[p]]$focal.list
                                               boot_progress.param(p,parameters.list = parameters.list)
                                               Boot_scans(Adj = Adj,n.boot = n.boot,total_scan = total_scan,obs.prob = obs.prob,keep = TRUE,
                                                          method = "both",focal.list = focal.list,scaled = TRUE,mode = mode,output = "adjacency",n.cores = 7)
                                             }
                     )
                     cat(paste0("\nCompiling data @ ",Sys.time(),"\n"))
                     Get.data(a,Bootstrap.list,parameters.list)
                   }
)


# Draft of data handling, plotting and analysis ---------------------------
library(data.table)
data.long<- data.table(data.long)
library(ggplot2)
mytheme<- theme_bw()+theme(plot.title = element_text(lineheight=.9, face="bold"),axis.line = element_line(colour = "black"),panel.grid.major.x = element_line(colour = "grey95",size = .2),panel.grid.minor.x = element_line(colour = "grey95",linetype = "dotted"),panel.grid.major.y =element_line(colour = "grey95",size = .2),panel.grid.minor.y=element_blank())+theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.75))

data.long$obs.prob.type<- as.factor(substr(data.long$obs.prob,1,3))
data.long$obs.prob.details<- as.factor(
  substr(data.long$obs.prob,
         nchar(as.character(data.long$obs.prob))-2,
         nchar(as.character(data.long$obs.prob))
  )
)

data.summary<- data.long[,.(cor=median(cor),sd.cor=sd(cor),
                            degree=median(degree),sd.degree=sd(degree),
                            strength=median(strength),sd.strength=sd(strength),
                            EV=median(EV),sd.EV=sd(EV),
                            CC=median(CC),sd.CC=sd(CC),
                            Frob=median(Frob),sd.Frob=sd(Frob),
                            Frob.GOF=median(Frob.GOF),sd.Frob.GOF=sd(Frob.GOF),
                            SLap=median(SLap),sd.SLap=sd(SLap),
                            SLap.GOF=median(SLap.GOF),sd.SLap.GOF=sd(SLap.GOF)),by = .(Network,obs.prob.type,obs.prob.details,focal.list,mode,method)]

# Matrix correlation
ggplot(data.summary[obs.prob.type=="net"],aes(interaction(method,obs.prob.type,obs.prob.details),cor,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = cor-sd.cor,ymax=cor+sd.cor),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type %in% c("net","tra")],aes(interaction(method,obs.prob.details),cor,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = cor-sd.cor,ymax=cor+sd.cor),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type=="unb"],aes(obs.prob.details,cor,colour = method,group=interaction(method,obs.prob.type)))+
  facet_grid(focal.list~Network)+geom_hline(yintercept = 1,lty="dashed",colour="grey50")+
  geom_linerange(aes(ymin = cor-sd.cor,ymax=cor+sd.cor))+geom_line()+geom_point(shape=21,fill="white")+
  scale_y_continuous(limits = c(min(data.summary[obs.prob.type=="unb"]$cor)-max(data.summary[obs.prob.type=="unb"]$sd.cor),1))+mytheme

# degree correlation
ggplot(data.summary[obs.prob.type=="net"],aes(interaction(method,obs.prob.type,obs.prob.details),degree,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = degree-sd.degree,ymax=degree+sd.degree),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type %in% c("net","tra")],aes(interaction(method,obs.prob.details),degree,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = degree-sd.degree,ymax=degree+sd.degree),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type=="unb"],aes(obs.prob.details,degree,colour = method,group=interaction(method,obs.prob.type)))+
  facet_grid(focal.list~Network)+geom_hline(yintercept = 1,lty="dashed",colour="grey50")+
  geom_linerange(aes(ymin = degree-sd.degree,ymax=degree+sd.degree))+geom_line()+geom_point(shape=21,fill="white")+
  scale_y_continuous(limits = c(min(data.summary[obs.prob.type=="unb"]$degree)-max(data.summary[obs.prob.type=="unb"]$sd.degree),1))+mytheme

# strength correlation
ggplot(data.summary[obs.prob.type=="net"],aes(interaction(method,obs.prob.type,obs.prob.details),strength,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = strength-sd.strength,ymax=strength+sd.strength),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type %in% c("net","tra")],aes(interaction(method,obs.prob.details),strength,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = strength-sd.strength,ymax=strength+sd.strength),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type=="unb"],aes(obs.prob.details,strength,colour = method,group=interaction(method,obs.prob.type)))+
  facet_grid(focal.list~Network)+geom_hline(yintercept = 1,lty="dashed",colour="grey50")+
  geom_linerange(aes(ymin = strength-sd.strength,ymax=strength+sd.strength))+geom_line()+geom_point(shape=21,fill="white")+
  scale_y_continuous(limits = c(min(data.summary[obs.prob.type=="unb"]$strength)-max(data.summary[obs.prob.type=="unb"]$sd.strength),1))+mytheme

# eigen-vector correlation
ggplot(data.summary[obs.prob.type=="net"],aes(interaction(method,obs.prob.type,obs.prob.details),EV,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = EV-sd.EV,ymax=EV+sd.EV),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type %in% c("net","tra")],aes(interaction(method,obs.prob.details),EV,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = EV-sd.EV,ymax=EV+sd.EV),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type=="unb"],aes(obs.prob.details,EV,colour = method,group=interaction(method,obs.prob.type)))+
  facet_grid(focal.list~Network)+geom_hline(yintercept = 1,lty="dashed",colour="grey50")+
  geom_linerange(aes(ymin = EV-sd.EV,ymax=EV+sd.EV))+geom_line()+geom_point(shape=21,fill="white")+
  scale_y_continuous(limits = c(min(data.summary[obs.prob.type=="unb"]$EV)-max(data.summary[obs.prob.type=="unb"]$sd.EV),1))+mytheme

# Clusering coefficient distance
ggplot(data.summary,aes(interaction(method,obs.prob.type,obs.prob.details),CC,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = CC-sd.CC,ymax=CC+sd.CC),colour="grey50",width = 0.2)+
  facet_grid(.~Network)+#facet_grid(obs.prob.type+focal.list~Network)+
  geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary,aes(interaction(method,obs.prob.details),CC,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = CC-sd.CC,ymax=CC+sd.CC),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+
  geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary,aes(obs.prob.details,CC,colour = method,group=interaction(method,obs.prob.type)))+
  facet_grid(focal.list~Network)+
  #geom_hline(yintercept = 1,lty="dashed",colour="grey50")+
  geom_linerange(aes(ymin = CC-sd.CC,ymax=CC+sd.CC))+geom_line()+geom_point(shape=21,fill="white")+
  mytheme

# GOF distance
ggplot(data.summary,aes(interaction(method,obs.prob.type,obs.prob.details),GOF,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = GOF-sd.GOF,ymax=GOF+sd.GOF),colour="grey50",width = 0.2)+
  facet_grid(.~Network)+#facet_grid(obs.prob.type+focal.list~Network)+
  geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary,aes(interaction(method,obs.prob.details),GOF,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = GOF-sd.GOF,ymax=GOF+sd.GOF),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+
  geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary,aes(obs.prob.details,GOF,colour = method,group=interaction(method,obs.prob.type)))+
  facet_grid(focal.list~Network)+
  geom_hline(yintercept = 1,lty="dashed",colour="grey50")+
  geom_linerange(aes(ymin = GOF-sd.GOF,ymax=GOF+sd.GOF))+geom_line()+geom_point(shape=21,fill="white")+
  mytheme

# Frob distance
ggplot(data.summary,aes(interaction(method,obs.prob.type,obs.prob.details),Frob,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = Frob-sd.frob,ymax=Frob+sd.frob),colour="grey50",width = 0.2)+
  facet_grid(.~Network)+#facet_grid(obs.prob.type+focal.list~Network)+
  geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary,aes(interaction(method,obs.prob.details),Frob,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = Frob-sd.frob,ymax=Frob+sd.frob),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+
  geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary,aes(obs.prob.details,Frob,colour = method,group=interaction(method,obs.prob.type)))+
  facet_grid(focal.list~Network)+
  geom_hline(yintercept = 1,lty="dashed",colour="grey50")+
  geom_linerange(aes(ymin = Frob-sd.frob,ymax=Frob+sd.frob))+geom_line()+geom_point(shape=21,fill="white")+
  mytheme


# Draft PCA ---------------------------------------------------------------

data.long.scaled<- data.long;
data.long.scaled[,c("cor","degree","strength","EV","ClustCoef","Frob","Frob.GOF","SLap","SLap.GOF")]<- data.table(scale(data.long.scaled[,c("cor","degree","strength","EV","ClustCoef","Frob","Frob.GOF","SLap","SLap.GOF")]))
PCA<- FactoMineR::PCA(data.long.scaled[Network=="1",c("cor","degree","strength","EV","ClustCoef","Frob","Frob.GOF","SLap","SLap.GOF")])
PCA$eig
PCA$var$contrib
PCA$var$coord

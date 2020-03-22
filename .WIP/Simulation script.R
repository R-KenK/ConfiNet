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
source(".WIP/ASNR.tools.R")


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

TOTAL_SCAN<- lapply(asnr.weighted.dir[which(sapply(asnr.list,length)==1)],get.total_scan)

with.total_scan<- !sapply(TOTAL_SCAN,is.null)
ADJ<- ADJ[with.total_scan]
TOTAL_SCAN<- TOTAL_SCAN[with.total_scan]

with.total_scan.inf1000<- sapply(TOTAL_SCAN,function(t) t<1000)
ADJ<- ADJ[with.total_scan.inf1000]
TOTAL_SCAN<- TOTAL_SCAN[with.total_scan.inf1000]
# Parameter choices -------------------------------------------------------

# Generate parameters list for each network once and for all --------------
PARAMETERS.LIST<- lapply(seq_along(ADJ),
                         function(a){
                           cat(paste0(a,"/",length(ADJ)," @ ",Sys.time(),"\n"))
                           Adj<- ADJ[[a]]
                           total_scan<- TOTAL_SCAN[[a]]
                           initialize_parameters(Adj,total_scan)
                         }
)


# Iterated Boot_scans() through parameters.list -------------------------
Boot.list<- lapply(seq_along(ADJ),
                   function(a){
                     cat(paste0(a,"/",length(ADJ)," @ ",Sys.time(),"\n"))
                     Adj<- ADJ[[a]]
                     total_scan<- TOTAL_SCAN[[a]]
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
                     Bootstrap.list
                   }
)

data.long<- rbind_lapply(seq_along(Boot.list),
                         function(B){
                           cat(paste0(B,"/",length(ADJ)," @ ",Sys.time(),"\n"))
                           Adj<- ADJ[[B]]
                           total_scan<- TOTAL_SCAN[[B]]
                           parameters.list<- PARAMETERS.LIST[[B]]
                           Bootstrap.list<- Boot.list[[B]]
                           Get.data(B,Bootstrap.list,parameters.list)
                         }
)

library(ggplot2)
mytheme<- theme_bw()+theme(plot.title = element_text(lineheight=.9, face="bold"),axis.line = element_line(colour = "black"),panel.grid.major.x = element_line(colour = "grey95",size = .2),panel.grid.minor.x = element_line(colour = "grey95",linetype = "dotted"),panel.grid.major.y =element_line(colour = "grey95",size = .2),panel.grid.minor.y=element_blank())+theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.75))
library(data.table)
data.long<- data.table(data.long)

data.long$obs.prob.type<- as.factor(substr(data.long$obs.prob,1,3))
data.long$obs.prob.details<- as.factor(
  substr(data.long$obs.prob,
         nchar(as.character(data.long$obs.prob))-2,
         nchar(as.character(data.long$obs.prob))
  )
)

data.summary<- data.long[,.(cor=median(cor),sd=sd(cor)),by = .(Network,obs.prob.type,obs.prob.details,focal.list,mode,method)]

ggplot(data.summary[obs.prob.type=="net"],aes(interaction(method,obs.prob.type,obs.prob.details),cor,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = cor-sd,ymax=cor+sd),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type %in% c("net","tra")],aes(interaction(method,obs.prob.details),cor,fill = method))+geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin = cor-sd,ymax=cor+sd),colour="grey50",width = 0.2)+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=1)+mytheme
ggplot(data.summary[obs.prob.type=="unb"],aes(obs.prob.details,cor,colour = method,group=interaction(method,obs.prob.type)))+
  facet_grid(focal.list~Network)+geom_hline(yintercept = 1,lty="dashed",colour="grey50")+
  geom_linerange(aes(ymin = cor-sd,ymax=cor+sd))+geom_line()+geom_point(shape=21,fill="white")+
  scale_y_continuous(limits = c(min(data.summary[obs.prob.type=="unb"]$cor)-0.1,1))+mytheme

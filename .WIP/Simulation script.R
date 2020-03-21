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
# Here preferably should be impleemnted as automatic import from ASNR

set.seed(42)
#
# n<- 10;nodes<- as.character(1:n);
# total_scan<- 200;
n.boot<- 10;
#
# Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
# Adj[non.diagonal(Adj)]<- sample((0:round(total_scan*.50)),n*(n-1),replace = TRUE)
# Adj

# devtools::install_github("schochastics/networkdata")
# library(networkdata)
# networkdata.list<- ls("package:networkdata")
# networkdata.list<- networkdata.list[
#   Reduce("|",
#          list(
#            grepl("^animal*",networkdata.list),
#            grepl("^ants*",networkdata.list),
#            grepl("^dolphin*",networkdata.list),
#            grepl("^giraffe*",networkdata.list),
#            grepl("^sheep*",networkdata.list),
#            grepl("^kangaroo*",networkdata.list),
#            grepl("^giraffe*",networkdata.list),
#            grepl("^macaque*",networkdata.list),
#            grepl("^rhesus*",networkdata.list)
#          )
#   )
#   ]
#
# networkdata.list<- paste0("networkdata::",networkdata.list)
# network.list<- lapply(networkdata.list,
#                       function(net){
#                         eval(parse(text = net))
#                       }
# )
#
# test<- sapply(network.list,
#               function(net){
#                 igraph::is.igraph(net)
#               }
# )
#
# sapply(unlist(network.list[!test],recursive = FALSE),
#        function(net){
#          igraph::is.igraph(net)
#        }
# )
# net.list<- c(unlist(network.list[!test],recursive = FALSE),network.list[test])
# length(net.list)
# weighted<- sapply(net.list,
#                   function(net){
#                     igraph::is_weighted(net)
#                   }
# )
# net.weighted<- net.list[weighted]
# sapply(net.weighted,
#        function(net){
#          igraph::edge.attributes(net)
#        }
# )[1:5]

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

# Parameter choices -------------------------------------------------------
# for each variable type, there should be only a non-nested list of parameters. I'll figure out later how to group similar values through factor ifelse() and substring I guess...
# ADJ<- lapply(net.weighted,
#              function(net){
#                Adj<- igraph::get.adjacency(net,attr = "weight",names = TRUE,sparse = FALSE)
#                rownames(Adj)<- as.character(1:nrow(Adj));colnames(Adj)<- as.character(1:ncol(Adj));
#                Adj
#              }
# )
#
# non.binary.sup.to.one<- sapply(ADJ,
#                                function(Adj){
#                                  Reduce("|",as.list(!(Adj %in% c(0,1)))) & Reduce("|",as.list(Adj>=1))
#                                }
# )
#
# ADJ<- ADJ[non.binary.sup.to.one]

initialize_parameters<- function(Adj,total_scan,n.cores=(parallel::detectCores()-1),cl=NULL){
  OBS.PROB<- list(trait.pos = obs.prob_bias(Adj = Adj,obs.prob_fun = prod,bias_fun = NULL,reverse = FALSE,cl = cl),
                  trait.neg = obs.prob_bias(Adj = Adj,obs.prob_fun = function(i,j) 1/prod(i,j),bias_fun = NULL,reverse = FALSE,cl = cl),
                  network.pos = obs.prob_bias(Adj = Adj,obs.prob_fun = prod,
                                              bias_fun = function(node) igraph::strength(igraph::graph.adjacency(Adj,weighted = TRUE))[node],
                                              reverse = FALSE,cl = cl),
                  network.neg = obs.prob_bias(Adj = Adj,obs.prob_fun = prod,
                                              bias_fun = function(node) 1/igraph::strength(igraph::graph.adjacency(Adj,weighted = TRUE))[node],
                                              reverse = FALSE,cl = cl)
  )
  OBS.PROB<- c({unb<- seq(0.1,0.9,by = 0.2);names(unb)<- paste0("unbiased_",unb);as.list(unb)},OBS.PROB) # c() over two lists makes them flat while allowing for shorter calls
  MODE<- if(identical(Adj,t(Adj))) {as.list(c("max"))} else  {as.list(c(directed = "directed",max = "max",min = "min",plus = "plus"))}
  FOCAL.LIST<- list(random = sample(rownames(Adj),total_scan,replace = TRUE),  # consider checking if can be implemented for each boot (leaving it NULL?)
                    even = rep_len(rownames(Adj),length.out = total_scan),
                    biased = "TO IMPLEMENT")

  parameters.comb<- expand.grid(
    list(mode = 1:length(MODE),
         focal.list = 1:length(FOCAL.LIST[1:2]),
         obs.prob = 1:length(OBS.PROB)
    )
  )

  if(is.null(cl)) {cl<- snow::makeCluster(n.cores);doSNOW::registerDoSNOW(cl);on.exit(snow::stopCluster(cl))} # left to avoid error if the function is used alone, but should probably be used internally from Boot_scans() now.
  parameters.list<- foreach::`%dopar%`(
    foreach::foreach(p=1:nrow(parameters.comb)),
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
  )
  parameters.list
}

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
sapply(ADJ,
       function(Adj){
         nrow(Adj)
       }
)

sapply(ADJ,
       function(Adj){
         round(max(Adj))
       }
)

Boot.list<- lapply(seq_along(ADJ),
                   function(a){
                     cat(a,"/",length(ADJ),"\n")
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
                                                          method = "both",focal.list = focal.list,scaled = TRUE,mode = mode,output = "all",n.cores = 7)
                                             }
                     )
                     Bootstrap.list
                   }
)

# Adj<- round(ADJ.test[[1]]/2)
# total_scan<- round(max(Adj)*1.25)
# parameters.list<- initialize_parameters(Adj,total_scan)
# n.boot<-1;
# start<- Sys.time()
# Bootstrap.list<- lapply(1,
#                         function(p){
#                           obs.prob<- parameters.list[[p]]$obs.prob;
#                           mode<- parameters.list[[p]]$mode;
#                           focal.list<- parameters.list[[p]]$focal.list
#                           boot_progress.param(p,parameters.list = parameters.list)
#                           Boot_scans(Adj = Adj,n.boot = n.boot,total_scan = total_scan,obs.prob = obs.prob,keep = TRUE,
#                                      method = "both",focal.list = NULL,scaled = TRUE,mode = mode,output = "all",n.cores = 7)
#                         }
# )
# stop<- Sys.time()
# stop-start
# Bootstrap.list

Bootstrap.list<- Boot.list[[1]]
parameters.list<- initialize_parameters(Adj,total_scan)

Get.data<- function(B,Bootstrap.list,parameters.list){
  rbind_lapply(seq_along(Bootstrap.list),
               function(b){
                 rbind(data.frame(Network = B,boot = b,method = "group",
                                  Boot_get.param(parameters.list[[b]]),
                                  cor = adjacency_cor(Bootstrap = Bootstrap.list[[b]],what = "observed")),
                       data.frame(Network = B,boot = b,method = "focal",
                                  Boot_get.param(parameters.list[[b]]),
                                  cor = adjacency_cor(Bootstrap = Bootstrap.list[[b]],what = "focal")))
               }
  )
}

data.long<- rbind_lapply(seq_along(Boot.list),
                         function(B){
                           cat(paste0(B,"/",length(ADJ.test)," @ ",Sys.time(),"\n"))
                           Adj<- ADJ.test[[B]]
                           total_scan<- round(max(Adj)*1.25)
                           parameters.list<- initialize_parameters(Adj,total_scan,n.cores = 7)
                           Bootstrap.list<- Boot.list[[B]]
                           Get.data(B,Bootstrap.list,parameters.list)
                         }
)

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

ggplot(data.long[obs.prob.type=="net"],aes(interaction(method,obs.prob.type,obs.prob.details),cor,colour = method))+
  facet_grid(Network~focal.list)+geom_jitter(alpha=0.2)+geom_boxplot()+mytheme
ggplot(data.long[obs.prob %in% c("network","trait")],aes(interaction(method,obs.prob),cor,colour = method))+facet_grid(mode~focal.list)+geom_jitter(alpha=0.2)+geom_boxplot()+theme_bw()
ggplot(data.long[!(obs.prob %in% c("network.pos","trait.pos","network.neg","trait.neg"))],aes(interaction(method,obs.prob),cor,colour = method))+facet_grid(mode~focal.list)+geom_jitter(alpha=0.2)+geom_boxplot()+theme_bw()

ggplot(data.summary[obs.prob.type=="net"],aes(interaction(method,obs.prob.type,obs.prob.details),cor,fill = method))+
  facet_grid(Network~focal.list)+geom_bar(stat = "identity",alpha=0.8)+geom_errorbar(aes(ymin = cor-sd,ymax=cor+sd),width = 0.2)+coord_flip()+theme_bw()
ggplot(data.summary[obs.prob.type %in% c("net","tra")],aes(interaction(method,obs.prob.details),cor,fill = method))+
  facet_grid(obs.prob.type+focal.list~Network)+geom_bar(stat = "identity",alpha=0.8)+geom_errorbar(aes(ymin = cor-sd,ymax=cor+sd),width = 0.2)+mytheme
ggplot(data.summary[obs.prob.type=="unb"],aes(obs.prob.details,cor,colour = method,group=interaction(method,obs.prob.type)))+
  facet_grid(focal.list~Network)+geom_linerange(aes(ymin = cor-sd,ymax=cor+sd))+geom_line()+geom_point(shape=21,fill="white")+
  scale_y_continuous(limits = c(min(data.summary[obs.prob.type=="unb"]$cor)-0.01,1))+mytheme

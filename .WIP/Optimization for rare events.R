# Draft of optimization for rare events ----------------------------------------------------------
rm(list = ls(envir = .GlobalEnv))
set.seed(42)

## import needed function, but ultimately should use library(ConfiNet) --------
source("R/matrix.tools.R")
source('R/rare.event.optimization_tools.R')
source("R/Misc_tools.R")
source("R/Binary.prob.R")
source("R/do.scan.R")
source("R/sum_up.scans.R")
source("R/iterate_scans.R")
source("R/Boot_scans.R")
source("R/observable_edges.R")
source("R/Bootstrap_tools.R")
source(".WIP/ASNR.tools.R")

Opti.test.Boot_scan.to.adj_wrapper<- function(Adj,total_scan,max.obs,method,use.rare.opti){
  Bootstrap<- Boot_scans(Adj = Adj,total_scan = total_scan,method = method,n.boot = 1,obs.prob = 0.42,scaled = TRUE,output = "all",use.rare.opti = use.rare.opti,cl = NULL)
  non.diagonal(Boot_get.list(Bootstrap,what = "theoretical",get.format = "adjacency")[[1]],output = "vector")
}

#' Optimization empirical support bootstrap function
#' Perform a bootstrap to estimate the coefficient of correlation and time of calculation of the theoretical networks with or without the optimization for rare event, and for given n, total_scan and max.obs. Parameters of Boot_scans() are set to do only one network estimation, but with the most demanding parameters(scaled = TRUE, method = "both", output = "all", mode, obs.prob and focal.list to be set by scan.default.args)
#'
#' @param n number of nodes in a simulated original network, to be approximated through iterating binary scans.
#' @param total_scan sampling effort, i.e. number of binary scan to perform
#' @param max.obs maximum number of observation possible within the simulated original network
#' @param boot number of bootstrap to perform (for the coefficient of correlation and time measure)
#'
#' @return a data.table with the coefficient of correlation and (median) time of calculation of the theoretical networks, obtained through the Boot_scan() function, with or without the optimization for rare event.
#' @export
#'
#' @examples
#' # internal use, cf hereafter.
cor.bootstrap<- function(n,total_scan,max.obs,boot=100){
  nodes<- as.character(1:n)
  Adj<- matrix(round(runif(n*n,0,max.obs)),nrow = n,ncol = n,dimnames = list(nodes,nodes))
  diag(Adj)<- 0
  rbind_lapply(1:boot,
                 function(b){
                   time<- microbenchmark::microbenchmark(
                     # standard.group = {Adj.std.group<- Opti.test.Boot_scan.to.adj_wrapper(Adj = Adj,total_scan = total_scan,max.obs = max.obs,method = "group",use.rare.opti = FALSE)},
                     # opti.group = {Adj.opti.group<- Opti.test.Boot_scan.to.adj_wrapper(Adj = Adj,total_scan = total_scan,max.obs = max.obs,method = "group",use.rare.opti = TRUE)},
                     # standard.focal = {Adj.std.focal<- Opti.test.Boot_scan.to.adj_wrapper(Adj = Adj,total_scan = total_scan,max.obs = max.obs,method = "focal",use.rare.opti = FALSE)},
                     # opti.focal = {Adj.opti.focal<- Opti.test.Boot_scan.to.adj_wrapper(Adj = Adj,total_scan = total_scan,max.obs = max.obs,method = "focal",use.rare.opti = TRUE)},
                     standard.all = {Adj.std.all<- Opti.test.Boot_scan.to.adj_wrapper(Adj = Adj,total_scan = total_scan,max.obs = max.obs,method = "both",use.rare.opti = FALSE)},
                     opti.all = {Adj.opti.all<- Opti.test.Boot_scan.to.adj_wrapper(Adj = Adj,total_scan = total_scan,max.obs = max.obs,method = "both",use.rare.opti = TRUE)},
                     times = 1,unit = "ms"
                   )
                   time



                   cor<- data.table::data.table(n = n,total_scan = total_scan,max.obs = max.obs,
                                                algorithm = rep_len(c("standard","opti"),length.out = 2),
                                                method = "all",
                                                # method = rep(c("group","focal","both"),each = 2),
                                                cor=c(
                                                  # cor(non.diagonal(Adj,output = "vector"),Adj.std.group),
                                                  # cor(non.diagonal(Adj,output = "vector"),Adj.opti.group),
                                                  # cor(non.diagonal(Adj,output = "vector"),Adj.std.focal),
                                                  # cor(non.diagonal(Adj,output = "vector"),Adj.opti.focal),
                                                  cor(non.diagonal(Adj,output = "vector"),Adj.std.all),
                                                  cor(non.diagonal(Adj,output = "vector"),Adj.opti.all)
                                                )
                   )


                   cbind(cor,time=summary(time)[["median"]],boot=b)
                 }
  )
}



# Parameters list ---------------------------------------------------------
param.comb<- expand.grid(n=seq(4,50,by=2),
                         total_scan=c(500,1000,5000,10000,25000),
                         max.obs=c(5,10,50,100,250))

parameters.list<- lapply(1:nrow(param.comb),
                         function(p){
                           c(n=param.comb[p,]$n,
                             total_scan=param.comb[p,]$total_scan,
                             max.obs=param.comb[p,]$max.obs)
                         }
)

# Validation (and comparison) of the optimizations ------------------------------------------
## Run standard method of doing a scan for each i of 1:total_scan, as well as optimized ones where only scans are performed when the scan is expected to be non full-zero
## then record coefficient of correlation of cumulated sum of the bernouilli trials with original vector of observations
## as well as the time to perform each method

# cl<- snow::makeCluster(7);snow::clusterExport(cl = cl,list = ls(envir = .GlobalEnv));doSNOW::registerDoSNOW(cl)
system.time(
  cor.boot<- rbind_lapply(seq_along(parameters.list),
                            # cor.boot<- rbind_pblapply(seq_along(parameters.list),
                            function(p){
                            n<- parameters.list[[p]]["n"];
                            total_scan<- parameters.list[[p]]["total_scan"];
                            max.obs<- parameters.list[[p]]["max.obs"];
                            cat(paste0("Param: n= ",n," - total_scan= ",total_scan," - max.obs= ",max.obs," (",p,"/",length(parameters.list),") @ ",Sys.time(),"\n"))
                            cor.bootstrap(n = n,total_scan = total_scan,max.obs = max.obs,boot = 30)
                            # },cl = cl
                          }
  )
)
# snow::stopCluster(cl)
cor.boot<- cor.boot[order(algorithm)][,c("algorithm","n","total_scan","max.obs","method","cor","time","boot"),with=FALSE]
cor.boot$algorithm<- factor(cor.boot$algorithm,levels = c("standard","opti"))
cor.boot$method<- factor(cor.boot$method,levels = c("group","focal","both"))

# saveRDS(cor.boot,".WIP/cor.boot.rds")
# cor.boot<- readRDS(".WIP/cor.boot.rds")

library(ggplot2)
library(data.table)
ggplot(cor.boot,aes(interaction(algorithm,total_scan),cor,colour=algorithm))+
  facet_grid(.~max.obs)+
  geom_jitter(alpha=0.2)+geom_boxplot(alpha=0.8)+theme_bw()

cor.summary<- cor.boot[,.(cor=median(cor),perc.5=quantile(cor,.05),perc.95=quantile(cor,.95)),by=.(algorithm,n,total_scan,max.obs)]
cor.summary$algorithm<- relevel(factor(cor.summary$algorithm),"standard")
cor.summary$total_scan.text<- factor(paste0("total_scan = ",cor.summary$total_scan),levels = paste0("total_scan = ",unique(cor.summary$total_scan)))
cor.summary$max.obs.text<- factor(paste0("max.obs = ",cor.summary$max.obs),levels = paste0("max.obs = ",unique(cor.summary$max.obs)))
cor.plot<- ggplot(cor.summary,aes(n,cor,colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,scales = "free")+
  geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+geom_point(alpha=1)+theme_bw()
cor.plot
ggplot(cor.summary,aes(n,cor,colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan~max.obs,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1,position = position_dodge(50/3-1))+
  geom_line(position = position_dodge(50/3-1))+geom_point(alpha=1,position = position_dodge(50/3-1))+theme_bw()
ggsave(filename = ".WIP/optimization_cor.plot.pdf",plot = cor.plot,width = 20,height = 15,units = "in")


time.summary<- cor.boot[,.(time=median(time),perc.5=quantile(time,.05),perc.95=quantile(time,.95)),by=.(algorithm,n,total_scan,max.obs)]
time.summary$algorithm<- relevel(factor(time.summary$algorithm),"standard")
time.summary$total_scan.text<- factor(paste0("total_scan = ",time.summary$total_scan),levels = paste0("total_scan = ",unique(time.summary$total_scan)))
time.summary$max.obs.text<- factor(paste0("max.obs = ",time.summary$max.obs),levels = paste0("max.obs = ",unique(time.summary$max.obs)))
time.plot<- ggplot(time.summary,aes(n,time,colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,scales = "free")+
  geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+geom_point(alpha=1)+theme_bw()
time.plot
ggsave(filename = ".WIP/optimization_time.plot.pdf",plot = time.plot,width = 20,height = 15,units = "in")

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


# Saving experimental computation times (TO CLEAN) -----------------------------------

# saveRDS(cor.boot,".WIP/timemod/cor.boot.rds")
# cor.boot<- readRDS(".WIP/timemod/cor.boot.rds")
# cor.boot.missing<- readRDS(".WIP/timemod/cor.boot,missing.small.n.rds")
# cor.boot.full<- rbind(cor.boot,cor.boot.missing)
# saveRDS(cor.boot.full,".WIP/timemod/cor.boot.full.rds")
# cor.boot<- readRDS(".WIP/timemod/cor.boot.full.rds")
# cor.boot.by.2<- readRDS(".WIP/timemod/cor.boot.by.2.rds")
# cor.boot.full<- readRDS(".WIP/timemod/cor.boot.full.rds")
#
# cor.boot<- rbind(cor.boot.by.2,cor.boot.full)

# Modelling draft #########################################################
library(ggplot2)
library(data.table)

# Import working modelling data -------------------------------------------
cor.boot.by.2<- readRDS(".WIP/timemod/cor.boot.by.2.rds")
cor.boot.full<- readRDS(".WIP/timemod/cor.boot.full.rds")

cor.boot<- rbind(cor.boot.by.2,cor.boot.full)

cor.boot$algorithm<- relevel(factor(cor.boot$algorithm),"standard")
cor.boot$total_scan.text<- factor(paste0("total_scan = ",cor.boot$total_scan),levels = paste0("total_scan = ",unique(cor.boot$total_scan)))
cor.boot$max.obs.text<- factor(paste0("max.obs = ",cor.boot$max.obs),levels = paste0("max.obs = ",unique(cor.boot$max.obs)))

time.sample<- cor.boot[boot<30]
time.summary<- cor.boot[,.(time=mean(time),perc.5=quantile(time,.05),perc.95=quantile(time,.95),var = var(time)),by=.(algorithm,n,total_scan,max.obs,max.obs.text,total_scan.text)]
time.summary[,CI.range:=perc.95-perc.5,]


# exploratory graphs ------------------------------------------------------
ggplot(time.summary,aes(n,time,colour=total_scan,fill=total_scan,group=total_scan))+
  facet_wrap(algorithm~max.obs.text,ncol=length(unique(cor.boot$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_line()+
  geom_point(alpha=1)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()

ggplot(time.summary,aes(total_scan,time,colour=n,fill=n,group=n))+
  facet_wrap(algorithm~max.obs.text,ncol=length(unique(cor.boot$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_line()+
  geom_point(alpha=1)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()

ggplot(time.summary,aes(max.obs,time,colour=n,fill=n,group=n))+
  facet_wrap(algorithm~total_scan.text,ncol=length(unique(cor.boot$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_line()+
  geom_point(alpha=1)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()

ggplot(time.summary,aes(time,fill=algorithm))+facet_grid(algorithm~total_scan+max.obs)+geom_histogram(colour="grey50",bins = 10)+theme_bw()


# Modelling standard method -----------------------------------------------
library(glmmTMB)
time.standard<- time.summary[algorithm=="standard"]
n.poly.std<- poly(time.standard$n,3)
total.poly.std<- poly(time.standard$total_scan,3)
max.poly.std<- poly(time.standard$max.obs,3)
time.standard<- cbind(time.standard,n=n.poly.std,total_scan=total.poly.std,max.obs=max.poly.std)
time.standard[,12:20]<- as.data.frame(scale(time.standard[,12:20,with=FALSE]))

standard.glm<- glm(round(time)~1,data = time.standard,family = "poisson")
standard.glm<- step(standard.glm,scope = list(lower=as.formula("round(time)~1"),
                                              upper=as.formula("round(time)~n.1*n.2*n.3*total_scan.1*total_scan.2*total_scan.3*max.obs.1*max.obs.2*max.obs.3")),
                    direction = "forward")
summary(standard.glm)
plot(standard.glm,which = 1)

standard.wglm<- glm(round(time)~1,data = time.standard,family = "poisson",weights = 1/sqrt(var))
standard.wglm<- step(standard.wglm,
                     scope = list(
                       lower=as.formula("round(time)~1"),
                       upper=as.formula("round(time)~n.1*n.2*n.3*total_scan.1*total_scan.2*total_scan.3*max.obs.1*max.obs.2*max.obs.3")),
                     direction = "both")
summary(standard.wglm)
plot(standard.wglm,which = 1)

lmtest::lrtest(standard.glm,standard.wglm)

# predictions<- predict(standard.glm,se.fit = TRUE,type = "response")
predictions<- predict(standard.wglm,se.fit = TRUE,type = "response")
standard.toplot<- cbind(time.standard,pred = predictions$fit,sd = predictions$se.fit)

# Modelling standard method -----------------------------------------------
library(glmmTMB)
time.opti<- time.summary[algorithm=="opti"]
n.poly.opt<- poly(time.opti$n,3)
total.poly.opt<- poly(time.opti$total_scan,3)
max.poly.opt<- poly(time.opti$max.obs,3)
time.opti<- cbind(time.opti,n=n.poly.opt,total_scan=total.poly.opt,max.obs=max.poly.opt)
time.opti[,12:20]<- as.data.frame(scale(time.opti[,12:20,with=FALSE]))

opti.glm<- glm(round(time)~1,data = time.opti,family = "poisson")
opti.glm<- step(opti.glm,scope = list(lower=as.formula("round(time)~1"),
                                              upper=as.formula("round(time)~n.1*n.2*n.3*total_scan.1*total_scan.2*total_scan.3*max.obs.1*max.obs.2*max.obs.3")),
                    direction = "forward")
summary(opti.glm)
plot(opti.glm,which = 1)

opti.wglm<- glm(round(time)~1,data = time.opti,family = "poisson",weights = 1/sqrt(var))
opti.wglm<- step(opti.wglm,
                     scope = list(
                       lower=as.formula("round(time)~1"),
                       upper=as.formula("round(time)~n.1*n.2*n.3*total_scan.1*total_scan.2*total_scan.3*max.obs.1*max.obs.2*max.obs.3")),
                     direction = "both")
summary(opti.wglm)
plot(opti.wglm,which = 1)

lmtest::lrtest(opti.glm,opti.wglm)

# predictions<- predict(opti.glm,se.fit = TRUE,type = "response")
predictions<- predict(opti.wglm,se.fit = TRUE,type = "response")
opti.toplot<- cbind(time.opti,pred = predictions$fit,sd = predictions$se.fit)


# Gather all result to plot -----------------------------------------------
time.toplot<- rbind(standard.toplot,opti.toplot)

time.plot<- ggplot(time.summary,aes(n,time,colour=total_scan,fill=total_scan,group=total_scan))+
  facet_wrap(algorithm~max.obs.text,ncol=length(unique(cor.boot$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4)+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),colour=NA,alpha=0.2)+
  geom_point(alpha=1)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()
time.plot
ggsave(filename = ".WIP/timemod/expected.time.models.plot.pdf",plot = time.plot,width = 20,height = 10,units = "in",dpi = 300)

time.plot.free<- ggplot(time.summary,aes(n,time,colour=total_scan,fill=total_scan,group=total_scan))+
  facet_wrap(algorithm~max.obs.text,ncol=length(unique(cor.boot$max.obs)),scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4)+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),colour=NA,alpha=0.2)+
  geom_point(alpha=1)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()
time.plot.free
ggsave(filename = ".WIP/timemod/expected.time.models.plot.free.pdf",plot = time.plot.free,width = 20,height = 10,units = "in",dpi = 300)

ggplot(time.summary,aes(total_scan,time,colour=n,fill=n,group=n))+
  facet_wrap(algorithm~max.obs.text,ncol=length(unique(cor.boot$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_point(alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4)+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),colour=NA,alpha=0.2)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()

# export standard and opti models into internal data -----------------------------------------
usethis::use_data(opti.model,n.poly.opt,total.poly.opt,max.poly.opt,standard.model,n.poly.std,total.poly.std,max.poly.std,internal = TRUE,overwrite = TRUE)

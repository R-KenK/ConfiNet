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
param.comb<- expand.grid(n=c(10,50,100,150,200,500,1000,1500),
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

cl<- snow::makeCluster(7);snow::clusterExport(cl = cl,list = ls(envir = .GlobalEnv));doSNOW::registerDoSNOW(cl)
system.time(
  cor.boot<- rbind_pblapply(seq_along(parameters.list),
                          function(p){
                            n<- parameters.list[[p]]["n"];
                            total_scan<- parameters.list[[p]]["total_scan"];
                            max.obs<- parameters.list[[p]]["max.obs"];
                            cat(paste0("Param: n= ",n," - total_scan= ",total_scan," - max.obs= ",max.obs," (",p,"/",length(parameters.list),") @ ",Sys.time(),"\n"))
                            cor.bootstrap(n = n,total_scan = total_scan,max.obs = max.obs,boot = 30)
                          },cl = cl
  )
)
snow::stopCluster(cl)
cor.boot<- cor.boot[order(algorithm)][,c("algorithm","n","total_scan","max.obs","method","cor","time","boot"),with=FALSE]
cor.boot$algorithm<- factor(cor.boot$algorithm,levels = c("standard","opti"))
cor.boot$method<- factor(cor.boot$method,levels = c("group","focal","both"))

# saveRDS(cor.boot,".WIP/cor.boot.rds")
# cor.boot<- readRDS(".WIP/cor.boot.rds")
# cor.boot.missing<- readRDS(".WIP/cor.boot,missing.small.n.rds")
# cor.boot.full<- rbind(cor.boot,cor.boot.missing)
# saveRDS(cor.boot.full,".WIP/cor.boot.full.rds")
# cor.boot<- readRDS(".WIP/cor.boot.full.rds")


# Modelling draft -----------------------------------------------------------
library(ggplot2)
library(data.table)
cor.boot<- readRDS(".WIP/cor.boot.full.rds")
cor.boot$algorithm<- relevel(factor(cor.boot$algorithm),"standard")
cor.boot$total_scan.text<- factor(paste0("total_scan = ",cor.boot$total_scan),levels = paste0("total_scan = ",unique(cor.boot$total_scan)))
cor.boot$max.obs.text<- factor(paste0("max.obs = ",cor.boot$max.obs),levels = paste0("max.obs = ",unique(cor.boot$max.obs)))

time.sample<- cor.boot#[boot<30]
time.summary<- cor.boot[,.(time=mean(time),perc.5=quantile(time,.05),perc.95=quantile(time,.95),var = var(time)),by=.(algorithm,n,total_scan,max.obs,max.obs.text,total_scan.text)]

ggplot(time.summary,aes(n,time,colour=total_scan,fill=total_scan,group=total_scan))+
  facet_wrap(algorithm~max.obs.text,ncol=length(unique(cor.boot$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_point(alpha=1)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()

ggplot(time.summary,aes(total_scan,time,colour=n,fill=n,group=n))+
  facet_wrap(algorithm~max.obs.text,ncol=length(unique(cor.boot$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_point(alpha=1)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()


library(glmmTMB)
time.sample<- time.summary[algorithm=="standard"]
n.poly<- poly(time.sample$n,2)
time.lm<- lm(time~total_scan*n.poly-1,data = time.sample)
summary(time.lm)
plot(time.lm,which = 1)

time.wlm<- lm(time~total_scan*n.poly-1,data = time.sample,weights = 1/var)
summary(time.wlm)
plot(time.wlm,which = 1)

time.glm<- glm(round(time)~total_scan*n.poly-1,data = time.sample,family = poisson(link = "log"))
summary(time.glm)
plot(time.glm,which = 1)

time.wglm<- glm(round(time)~total_scan+total_scan:n.poly-1,data = time.sample,family = poisson(link = "identity"),weights = 1/sqrt(var))
summary(time.wglm)
plot(time.wglm,which = 1)
plot(DHARMa::simulateResiduals(time.wglm))

lmtest::lrtest(time.lm,time.wlm,time.glm,time.wglm)

predict(time.lm,newdata = list(total_scan=1500,n.poly=predict(n.poly,60)))
predict(time.wlm,newdata = list(total_scan=1500,n.poly=predict(n.poly,60)))
predict(time.glm,newdata = list(total_scan=1500,n.poly=predict(n.poly,60)))
predict(time.wglm,newdata = list(total_scan=1500,n.poly=predict(n.poly,60)))

predictions<- predict(time.wglm,se.fit = TRUE,type = "response")
time.toplot<- cbind(time.sample,pred = predictions$fit,sd = predictions$se.fit)

time.sample<- time.summary[algorithm=="opti"]
n.poly<- poly(time.sample$n,2)
time.lm<- lm(time~total_scan*n.poly*max.obs-1,data = time.sample)
summary(time.lm)
plot(time.lm,which = 1)

time.wlm<- lm(time~total_scan*n.poly*max.obs-1,data = time.sample,weights = 1/var)
summary(time.wlm)
plot(time.wlm,which = 1)

time.glm<- glm(round(time)~total_scan*n.poly*max.obs-1,data = time.sample,family = poisson(link = "log"))
summary(time.glm)
plot(time.glm,which = 1)

time.wglm<- glm(round(time)~total_scan+total_scan:n.poly-1,data = time.sample,family = poisson(link = "log"),weights = 1/sqrt(var))
summary(time.wglm)
plot(time.wglm,which = 3)
plot(DHARMa::simulateResiduals(time.wglm))

lmtest::lrtest(time.lm,time.wlm,time.glm,time.wglm)

predict(time.lm,newdata = list(total_scan=1500,n.poly=predict(n.poly,60)))
predict(time.wlm,newdata = list(total_scan=1500,n.poly=predict(n.poly,60)))
predict(time.glm,newdata = list(total_scan=1500,n.poly=predict(n.poly,60)))
predict(time.wglm,newdata = list(total_scan=1500,n.poly=predict(n.poly,60)))

time.wglm<- glm(round(time)~total_scan*n.poly*max.obs-1,data = time.sample,family = poisson(link = "identity"))
summary(time.wglm)
plot(time.wglm,which = 1)
plot(DHARMa::simulateResiduals(time.wglm))


ggplot(time.summary,aes(n,time,colour=total_scan,fill=total_scan,group=total_scan))+
  facet_wrap(algorithm~max.obs.text,ncol=length(unique(cor.boot$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4)+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),colour=NA,alpha=0.2)+
  geom_point(alpha=1)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()


ggplot(time.summary[algorithm=="standard"],aes(total_scan,time,colour=n,fill=n,group=n))+
  facet_wrap(.~max.obs.text,ncol=length(unique(cor.boot[algorithm=="standard"]$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_point(alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4)+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),colour=NA,alpha=0.2)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()

plot(DHARMa::simulateResiduals(time.glm))


ggplot(time.sample,aes(n,time,colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,ncol=length(unique(time.sample$max.obs)),scales = "free")+
  geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+geom_point(alpha=1)+theme_bw()

ggplot(time.sample[algorithm=="standard"],aes(n,time,colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,ncol=length(unique(time.sample$max.obs)),scales = "free")+
  geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+geom_point(alpha=1)+theme_bw()

ggplot(time.summary[algorithm=="standard"],aes(total_scan,time,colour=n,fill=algorithm,group=n))+
  facet_wrap(.~max.obs.text,ncol=length(unique(time.sample$max.obs)),scales = "free")+
  # geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+
  geom_point(alpha=1)+theme_bw()

ggplot(time.summary[algorithm=="standard"],aes(total_scan,time,colour=n,fill=algorithm,group=n))+
  facet_wrap(.~max.obs.text,ncol=length(unique(time.summary$max.obs)),scales = "free")+
  # geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+
  geom_point(alpha=1)+theme_bw()

ggplot(time.summary[algorithm=="standard"],aes(n,time,colour=total_scan,fill=algorithm,group=total_scan))+
  facet_wrap(.~max.obs.text,ncol=length(unique(time.summary$max.obs)),scales = "free")+
  # geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+
  geom_point(alpha=1)+theme_bw()

ggplot(cor.boot[algorithm=="standard"],aes(max.obs,time,colour=n,fill=algorithm,group=n))+
  facet_wrap(.~total_scan,ncol=length(unique(time.summary$total_scan)),scales = "free")+
  # geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+
  geom_jitter(alpha=0.1)+theme_bw()

ggplot(cor.boot[algorithm=="standard"&boot<30],aes(max.obs,time,colour=n,fill=algorithm,group=n))+
  facet_wrap(.~total_scan,ncol=length(unique(time.summary$total_scan)))+
  # geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  # geom_line()+
  geom_point(alpha=0.1)+theme_bw()

ggplot(time.summary[algorithm=="standard"&total_scan==300],aes(n,time,colour=max.obs,fill=algorithm,group=max.obs))+
  # facet_wrap(.~max.obs.text,ncol=length(unique(time.summary$max.obs)),scales = "free")+
  # geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+
  geom_point(alpha=1)+theme_bw()

unique(time.summary$total_scan)


# standard case -----------------------------------------------------------
library(glmmTMB)
time.scaled<- time.summary[algorithm=="standard"]
time.glm.2<- glmmTMB(time~total_scan+total_scan:poly(n,2)-1,data = time.scaled,weights = 1/var)
summary(time.glm.2)

predictions<- predict(time.glm.2,se.fit = TRUE)
time.toplot<- cbind(time.scaled,pred = predictions$fit,sd = predictions$se.fit)

ggplot(cor.boot[algorithm=="standard"&boot<20],aes(n,time,colour=algorithm,fill=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,ncol=length(unique(cor.boot[algorithm=="standard"]$max.obs)),scales = "free")+
  geom_point(alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4,colour="red")+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),alpha=0.2,colour=NA,fill="red")+
  theme_bw()

ggplot(time.summary[algorithm=="standard"],aes(n,time,colour=total_scan,fill=total_scan,group=total_scan))+
  facet_wrap(.~max.obs.text,ncol=length(unique(cor.boot[algorithm=="standard"]$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_point(alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4)+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),colour=NA,alpha=0.2)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()

ggplot(time.summary[algorithm=="standard"],aes(total_scan,time,colour=n,fill=n,group=n))+
  facet_wrap(.~max.obs.text,ncol=length(unique(cor.boot[algorithm=="standard"]$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_point(alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4)+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),colour=NA,alpha=0.2)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()

plot(time.glm)

# opti case -----------------------------------------------------------
library(glmmTMB)
time.scaled<- time.summary[algorithm=="opti"]
Pol.n<- poly(time.scaled$n,2)
time.glm.2<- glm(time~total_scan+total_scan:n+total_scan:I(n^2)+total_scan:max.obs-1,data = time.scaled,weights = 1/(var))
summary(time.glm.2)

predictions<- predict(time.glm.2,se.fit = TRUE)
time.toplot<- cbind(time.scaled,pred = predictions$fit,sd = predictions$se.fit)

ggplot(time.summary[algorithm=="opti"],aes(n,time,colour=total_scan,fill=total_scan,group=total_scan))+
  facet_wrap(.~max.obs.text,ncol=length(unique(cor.boot[algorithm=="opti"]$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_point(alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4)+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),colour=NA,alpha=0.2)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()

ggplot(time.summary[algorithm=="opti"],aes(total_scan,time,colour=n,fill=n,group=n))+
  facet_wrap(.~max.obs.text,ncol=length(unique(cor.boot[algorithm=="opti"]$max.obs)))+#,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1)+
  geom_point(alpha=1)+
  geom_line(data=time.toplot,aes(y = pred),alpha=0.4)+
  geom_ribbon(data=time.toplot,aes(ymin = pred-1.96*sd,ymax = pred+1.96*sd),colour=NA,alpha=0.2)+
  scale_color_gradient(low="lightblue",high = "red")+
  scale_fill_gradient(low="lightblue",high = "red")+
  theme_bw()


predict(time.glm.2,newdata = data.table(algorithm="opti",n=40,total_scan=200000,max.obs=30))


# Draft -------------------------------------------------------------------


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


time.summary<- cor.boot[,.(time=mean(time),perc.5=quantile(time,.05),perc.95=quantile(time,.95),var = var(time)),by=.(algorithm,n,total_scan,max.obs)]
time.summary$algorithm<- relevel(factor(time.summary$algorithm),"standard")
time.summary$total_scan.text<- factor(paste0("total_scan = ",time.summary$total_scan),levels = paste0("total_scan = ",unique(time.summary$total_scan)))
time.summary$max.obs.text<- factor(paste0("max.obs = ",time.summary$max.obs),levels = paste0("max.obs = ",unique(time.summary$max.obs)))
time.plot<- ggplot(time.summary[n<50],aes(n,time,colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,ncol=length(unique(time.summary$max.obs)),scales = "free")+
  geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+geom_point(alpha=1)+theme_bw()
time.plot
ggsave(filename = ".WIP/optimization_time.plot.pdf",plot = time.plot,width = 20,height = 15,units = "in")


# Modelling expected time =f(n,total_scan,max.obs) ------------------------
time.plot.transformed<- ggplot(time.summary,aes(n,sqrt(time),colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,ncol=length(unique(time.summary$max.obs)),scales = "free")+
  geom_ribbon(aes(ymin = sqrt(perc.5),ymax = sqrt(perc.95)),colour=NA,alpha=0.3)+
  geom_line()+geom_point(alpha=1)+theme_bw()
time.plot.transformed
ggsave(filename = ".WIP/optimization_time.transformed.plot.pdf",plot = time.plot.transformed,width = 20,height = 15,units = "in")

time.sqrt.lm<-lm(sqrt(time)~n+total_scan+max.obs,data = time.summary[algorithm=="standard"],weights = 1/var)
summary(time.sqrt.lm)
plot(time.sqrt.lm)

time.sqrt.lm<-lm(sqrt(time)~n+total_scan+max.obs,data = time.summary[algorithm=="opti"],weights = 1/var)
summary(time.sqrt.lm)
plot(time.sqrt.lm)

library(glmmTMB)
time.scaled<- time.summary
time.scaled$n<- scale(time.scaled$n)[seq_along(time.scaled$n)]
time.scaled$total_scan<- scale(time.scaled$total_scan)[seq_along(time.scaled$total_scan)]
time.scaled$max.obs<- scale(time.scaled$max.obs)[seq_along(time.scaled$max.obs)]

fitdistrplus::descdist((time.scaled$time),boot = 100)

time.summary$n2<- time.summary$n*time.summary$n

time.glm<- glm(round(time)~n2+total_scan+max.obs,data = time.summary[algorithm=="opti"],weights = 1/var,family = poisson(link = "sqrt"))
summary(time.glm)
plot(time.glm)

time.plot<- ggplot(time.summary[algorithm=="opti"],aes(n,time,colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,ncol=length(unique(time.summary$max.obs)),scales = "free")+
  geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+geom_point(alpha=1)+theme_bw()
time.plot

time.opti<- cor.boot[algorithm=="opti"][max.obs==5]
time.opti$algorithm<- relevel(factor(time.opti$algorithm),"standard")
time.opti$total_scan.text<- factor(paste0("total_scan = ",time.opti$total_scan),levels = paste0("total_scan = ",unique(time.opti$total_scan)))
time.opti$max.obs.text<- factor(paste0("max.obs = ",time.opti$max.obs),levels = paste0("max.obs = ",unique(time.opti$max.obs)))
time.weight<- time.opti[,.(w = 1/(max(time)-min(time))),by=.(n,total_scan,max.obs)]
time.opti<- merge(time.opti,time.weight,by = c("n","total_scan","max.obs"))

time.lm<- lm(time~n+I(n^3)-1,data = time.opti,weights = w)
time.lm<- lm(time~poly(n,degree = 4,raw = TRUE)-1,data = time.opti,weights = w)

time.opti.summary<- time.opti[,.(logtime=median(log(time)),logn=log(n),w = 1/var(log(time))),by=.(algorithm,n,total_scan,max.obs,total_scan.text,max.obs.text)]
time.loglog<- glmmTMB(logtime~1+(logn:total_scan|logn:total_scan),data = time.opti.summary)
summary(time.loglog)

predictions<- data.table(predict(time.loglog,interval = "prediction"))
confidence<- data.table(predict(time.loglog,interval = "confidence"))
time.toplot<- cbind(time.opti.summary,pred = predictions,conf = confidence)

time.plot<- ggplot(time.toplot,aes(log(n),logtime,colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,ncol=length(unique(time.summary$max.obs)),scales = "free")+
  # geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_point(alpha=1)+
  geom_line(aes(y = pred.fit),alpha=0.4,colour="red")+
  geom_ribbon(aes(ymin = pred.lwr,ymax = pred.upr),alpha=0.2,colour=NA,fill="red")+
  geom_line(aes(y = conf.fit),alpha=0.4,colour="blue")+
  geom_ribbon(aes(ymin = conf.lwr,ymax = conf.upr),alpha=0.2,colour=NA,fill="blue")+
  theme_bw()
time.plot
plot(time.loglog)



plot(time.summary[algorithm=="opti"]$n,time.summary[algorithm=="opti"]$time)
lines(time.summary[algorithm=="opti"]$n,predict(time.nls),lty=2,col="red",lwd=3)

time.plot.transformed<- ggplot(time.scaled[algorithm=="standard"],aes(n,sqrt(time),colour=algorithm,fill=algorithm,group=algorithm))+
  facet_wrap(total_scan.text~max.obs.text,ncol=length(unique(time.scaled$max.obs)),scales = "free")+
  geom_boxplot()+theme_bw()
time.plot.transformed

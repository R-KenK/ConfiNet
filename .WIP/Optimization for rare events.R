# Draft of optimization for rare events ----------------------------------------------------------
set.seed(42)

## import needed function, but ultimately should use library(ConfiNet) --------
source("R/matrix.tools.R")
source("R/Binary.prob.R")
source("R/Bootstrap_tools.R")
library(data.table)
library(microbenchmark)

#' Title
#'
#' @param P
#' @param ...
#' @param n
#'
#' @return
#' @export
#'
#' @examples
P.cond.in.order<- function(P,...,n = length(P)){
  S<- 1-prod(1-P)
  previous.are.zeros<- c(1,cumprod(1-P[1:(n-1)]))
  sapply(1:n,function(i) P[i]/S*previous.are.zeros[i])
}

#' Title
#'
#' @param P
#' @param ...
#' @param n
#'
#' @return
#' @export
#'
#' @examples
do.non.zero.scan<-function(P,...,n = length(P)){
  X<- rep(0,n)
  P.cond<- cumsum(P.cond.in.order(P))
  first<- min(which(runif(1)<P.cond))
  X[first]<- 1
  X[-first]<- rbinom(n-1,1,P[-first])
  X
}

do.non.zero.scan.bis<-function(P,...,n = length(P)){
  X<- rep(0,n)
  rand.order<- sample(1:n,n)
  P.cond<- cumsum(P.cond.in.order(P[rand.order]))
  first<- min(which(runif(1)<P.cond))
  X[rand.order][first]<- 1
  X[rand.order][-first]<- rbinom(n-1,1,P[rand.order][-first])
  X
}

#' Title
#'
#' @param P
#' @param N
#' @param ...
#' @param n
#'
#' @return
#' @export
#'
#' @examples
iterate_rare.scans<- function(P,N,...,n = length(P)){
  scan<- data.table::data.table(matrix(0,N,n));
  non.zero<- rbinom(N,1,1-prod(1-P))==1;
  scan[non.zero,]<- data.table::data.table(rbind_lapply(seq_len(nrow(scan[non.zero,])),function(k) do.non.zero.scan(P)))
  scan
}

iterate_rare.scans.bis<- function(P,N,...,n = length(P)){
  scan<- data.table(matrix(0,N,n));
  non.zero<- rbinom(N,1,1-prod(1-P))==1;
  scan[non.zero,]<- data.table(rbind_lapply(seq_len(nrow(scan[non.zero,])),function(k) do.non.zero.scan.bis(P)))
  scan
}

#' Title
#'
#' @param P
#' @param N
#' @param ...
#' @param n
#'
#' @return
#' @export
#'
#' @examples
iterate_standard.scans<- function(P,N,...,n = length(P)){
  rbind_lapply(1:N,
               function(k){
                 X<- rbinom(n,1,P)
                 if(Reduce("&",X==0)) rep(0,n) else X
               }
  )
}

#' Title
#'
#' @param n
#' @param N
#' @param max.obs
#' @param boot
#'
#' @return
#' @export
#'
#' @examples
cor.bootstrap<- function(n,N,max.obs,boot=100,n.cores=1,.export){
  V<- round(runif(n,0,max.obs))
  P<- Binary.prob(V,N,mode = "vector")$present
  rbind_pblapply(1:boot,
                 function(b){
                   time<- microbenchmark(
                     standard = {X<- iterate_standard.scans(P,N)},
                     opti.ordered = {X.opti<- iterate_rare.scans(P,N)},
                     opti.random = {X.opti.bis<- iterate_rare.scans.bis(P,N)},
                     times = 1,unit = "ms"
                   )
                   time

                   cor<- data.table(n=n,N=N,max.obs=max.obs,method=c("standard","opti.ordered","opti.random"),
                                    cor=c(cor(colSums(X),V),cor(colSums(X.opti),V),cor(colSums(X.opti.bis),V)))

                   cbind(cor,time=summary(time)[["median"]],boot=b)
                 },n.cores = n.cores,.export = .export
  )
}



# Parameters list ---------------------------------------------------------
param.comb<- expand.grid(n=c(10,50,100,150,200,500,1000),
                         N=c(500,1000,10000,50000),
                         max.obs=c(5,10,50,100,250))

parameters.list<- lapply(1:nrow(param.comb),
                         function(p){
                           c(n=param.comb[p,]$n,
                             N=param.comb[p,]$N,
                             max.obs=param.comb[p,]$max.obs)
                         }
)


# Validation (and comparison) of the optimizations ------------------------------------------
## Run standard method of doing a scan for each i of 1:N, as well as optimized ones where only scans are performed when the scan is expected to be non full-zero
## then record coefficient of correlation of cumulated sum of the bernouilli trials with original vector of observations
## as well as the time to perform each method

cor.boot<- rbind_lapply(seq_along(parameters.list),
                        function(p){
                          n<- parameters.list[[p]]["n"];
                          N<- parameters.list[[p]]["N"];
                          max.obs<- parameters.list[[p]]["max.obs"];
                          cat(paste0("Param: n= ",n," - N= ",N," - max.obs= ",max.obs," (",p,"/",length(parameters.list),")\n"))
                          cor.bootstrap(n,N,max.obs,boot = 100,n.cores = 7,
                                        .export = c("parameters.list","cor.bootstrap",
                                                    "iterate_standard.scans","iterate_rare.scans","iterate_rare.scans.bis",
                                                    "do.non.zero.scan","do.non.zero.scan.bis","P.cond.in.order",
                                                    "Binary.prob","rbind_lapply","data.table","microbenchmark"))
                        }

)

cor.boot<- cor.boot[order(method)]
cor.boot$method<- factor(cor.boot$method,levels = c("standard","opti.ordered","opti.random"))

# saveRDS(cor.boot,".WIP/cor.boot.rds")
# cor.boot<- readRDS(".WIP/cor.boot.rds")

library(ggplot2)
library(data.table)
ggplot(cor.boot,aes(interaction(method,N),cor,colour=method))+
  facet_grid(.~max.obs)+
  geom_jitter(alpha=0.2)+geom_boxplot(alpha=0.8)+theme_bw()

cor.summary<- cor.boot[,.(cor=median(cor),perc.5=quantile(cor,.05),perc.95=quantile(cor,.95)),by=.(method,n,N,max.obs)]
cor.plot<- ggplot(cor.summary,aes(n,cor,colour=method,fill=method,group=method))+
  facet_wrap(N~max.obs,scales = "free")+
  geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+geom_point(alpha=1)+theme_bw()
cor.plot
ggplot(cor.summary,aes(n,cor,colour=method,fill=method,group=method))+
  facet_wrap(N~max.obs,scales = "free")+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1,position = position_dodge(50/3-1))+
  geom_line(position = position_dodge(50/3-1))+geom_point(alpha=1,position = position_dodge(50/3-1))+theme_bw()
ggsave(filename = ".WIP/optimization_cor.plot.pdf",plot = cor.plot,width = 12,height = 7,units = "in")

time.summary<- cor.boot[,.(time=median(as.numeric(time)),perc.5=quantile(as.numeric(time),.05),perc.95=quantile(as.numeric(time),.95)),by=.(method,n,N,max.obs)]
time.plot<- ggplot(time.summary,aes(n,time,colour=method,fill=method,group=method))+
  facet_wrap(N~max.obs,scales = "free")+
  geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour=NA,alpha=0.3)+
  geom_line()+geom_point(alpha=1)+theme_bw()
time.plot
ggplot(time.summary,aes(n,time,colour=method,fill=method,group=method))+
  facet_grid(N~max.obs)+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1,position = position_dodge(50/3-1))+
  geom_line(position = position_dodge(50/3-1))+geom_point(alpha=1,position = position_dodge(50/3-1))+theme_bw()
ggsave(filename = ".WIP/optimization_time.plot.pdf",plot = time.plot,width = 12,height = 7,units = "in")

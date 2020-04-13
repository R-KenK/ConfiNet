# Draft of optimization for rare events ----------------------------------------------------------
set.seed(42)
n<- 100;N<- 10000
V<- round(runif(n,0,10))
P<- Binary.prob(V,N)

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
  scan<- data.table(matrix(0,N,n));
  non.zero<- rbinom(N,1,1-prod(1-P))==1;
  scan[non.zero,]<- data.table(rbind_lapply(seq_len(nrow(scan[non.zero,])),function(k) do.non.zero.scan(P)))
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
                 X<-rbinom(n,1,P)
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
cor.bootstrap<- function(n,N,max.obs,boot){
  rbind_pblapply(1:boot,
                 function(b){
                   n<- n;N<- N;max.obs<- max.obs;
                   V<- round(runif(n,0,max.obs)) # because max.obs << N, this should be a scenario where the optimization is relevant
                   P<- Binary.prob(V,N)

                   start.std<- Sys.time()
                   X<- iterate_standard.scans(P,N)
                   stop.std<- Sys.time()

                   start.opt<- Sys.time()
                   X.opti<- iterate_rare.scans(P,N)
                   stop.opt<- Sys.time()

                   start.opt.bis<- Sys.time()
                   X.opti.bis<- iterate_rare.scans.bis(P,N)
                   stop.opt.bis<- Sys.time()

                   rbind(
                     data.table(n=n,N=N,max.obs=max.obs,method="standard",cor=cor(colSums(X),V),time = stop.std-start.std),
                     data.table(n=n,N=N,max.obs=max.obs,method="opti.ordered",cor=cor(colSums(X.opti),V),time = stop.opt-start.opt),
                     data.table(n=n,N=N,max.obs=max.obs,method="opti.random",cor=cor(colSums(X.opti.bis),V),time = stop.opt.bis-start.opt.bis)
                   )
                 },
                 n.cores = n.cores,.export = c("iterate_standard.scans","iterate_rare.scans","iterate_rare.scans.bis",
                                               "do.non.zero.scan","do.non.zero.scan.bis",
                                               "P.cond.in.order","Binary.prob","rbind_lapply","data.table")
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
                          cor.bootstrap(n,N,max.obs,100)
                        }
)

cor.boot<- cor.boot[order(method)]
cor.boot$method<- factor(cor.boot$method,levels = c("standard","opti.ordered","opti.random"))

# saveRDS(cor.boot,".WIP/cor.boot.rds")

library(ggplot2)
ggplot(cor.boot,aes(interaction(method,N),cor,colour=method))+
  facet_grid(.~max.obs)+
  geom_jitter(alpha=0.2)+geom_boxplot(alpha=0.8)+mytheme
ggplot(cor.boot,aes(rareness,cor,colour=method,group=interaction(rareness,method)))+
  geom_jitter(alpha=0.2)+geom_boxplot(alpha=0.8)+mytheme

cor.summary<- cor.boot[,.(cor=mean(cor),perc.5=quantile(cor,.05),perc.95=quantile(cor,.95)),by=.(method,n,N,max.obs)]
ggplot(cor.summary,aes(n,cor,colour=method,fill=method,group=method))+
  facet_grid(N~max.obs)+
  geom_ribbon(aes(ymin = perc.5,ymax = perc.95),colour="white",alpha=0.3)+
  geom_line()+geom_point(alpha=1)+mytheme
ggplot(cor.summary,aes(n,cor,colour=method,fill=method,group=method))+
  facet_grid(N~max.obs)+
  geom_linerange(aes(ymin = perc.5,ymax = perc.95),alpha=1,position = position_dodge(50/3-1))+
  geom_line(position = position_dodge(50/3-1))+geom_point(alpha=1,position = position_dodge(50/3-1))+mytheme


cor.null<- glm(cor~1,data = cor.boot,family = "binomial",weights = rep(n,nrow(cor.boot)))     # Considering the coefficient of correlation is based on n points per calculation of cor
cor.glm<- glm(cor~method,data = cor.boot,family = "binomial",weights = rep(n,nrow(cor.boot)))

summary(cor.null)
summary(cor.glm)

lmtest::lrtest(cor.null,cor.glm)

ggplot(cor.boot,aes(method,time,colour=method))+geom_jitter(alpha=0.2)+geom_boxplot(alpha=0.8)+mytheme
time.null<- lm(as.numeric(time)~1,data = cor.boot)
time.glm<- lm(as.numeric(time)~method,data = cor.boot)

summary(time.null)
summary(time.glm)

lmtest::lrtest(time.null,time.glm)

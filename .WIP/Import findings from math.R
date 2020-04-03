source("R/matrix.tools.R")
source("R/Binary.prob.R")
source("R/do.scan.R")
source("R/sum_up.scans.R")
source("R/iterate_scans.R")
source("R/Boot_scans.R")
source("R/observable_edges.R")
source("R/Bootstrap_tools.R")

is.all.equal <- function(x) {diff(range(x)) < .Machine$double.eps ^ 0.5} # from https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector, efficiently test if vector contain all similar value by comparing the difference between its range and the square root of the internal minimal value

set.seed(42)

n<- 20;nodes<- as.character(1:n);
total_scan<- 10000;
Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
Adj[non.diagonal(Adj)]<- sample(0:5,n*(n-1),replace = TRUE)
Adj<- adjacency_mode(Adj,"min")
focal.list<- sample(nodes,total_scan,replace = TRUE)
table(focal.list)


# Findings from probability analysis --------------------------------------

Bin.prob<- Binary.prob(Adj,total_scan,mode = "max")

(1-prod(Bin.prob$absent))
list.of.non.zero<- sample(c(TRUE,FALSE),total_scan,replace = TRUE,prob = c(1-prod(Bin.prob$absent),prod(Bin.prob$absent)))
table(list.of.non.zero)

P<- matrix(0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
P[non.diagonal(P)]<- Bin.prob$present
prob.test<- non.diagonal(P,"vector")


library(microbenchmark)

test.runif<- function(prob.test,total_scan){
  rbind_lapply(seq_len(total_scan),
               function(k){
                 scan<- as.numeric(runif(length(prob.test)) < prob.test)
                 if(is.all.equal(scan)){
                   rep(0,length(prob.test))
                 }
                 else{
                   scan
                 }
               }
  )
}

test.rbinom<- function(prob.test,total_scan){
  rbind_lapply(seq_len(total_scan),
               function(k){
                 scan<- rbinom(length(prob.test),1,prob.test)
                 if(is.all.equal(scan)){
                   rep(0,length(prob.test))
                 }
                 else{
                   scan
                 }
               }
  )
}

test.sample<- function(prob.test,total_scan){
  rbind_lapply(seq_len(total_scan),
               function(k){
                 scan<- sapply(seq_along(prob.test),function(d) sample(c(1,0),1,replace = TRUE,prob=c(prob.test[d],1-prob.test[d])))
                 if(is.all.equal(scan)){
                   rep(0,length(prob.test))
                 }
                 else{
                   scan
                 }
               }
  )
}

microbenchmark(test.rbinom(prob.test,total_scan),test.runif(prob.test,total_scan),test.sample(prob.test,total_scan),times = 5)

as.numeric(runif(length(prob.test)) < prob.test)

test.runif<- vector(mode="list",length = length(list.of.non.zero))
system.time(
  test.runif<- rbind_lapply(seq_len(total_scan),
                             function(k){
                               scan<- as.numeric(runif(length(prob.test)) < prob.test)
                               if(is.all.equal(scan)){
                                 rep(0,length(prob.test))
                               }
                               else{
                                 scan
                               }
                             }
  )
)


scan.list<- sapply(list.of.non.zero,
                   function(is.non.zero){
                     if(is.non.zero){
                       "I COMPUTED A SCAN"
                     }else{
                       "MAtrIX IS ZERO LOL"
                     }
                   }
)
table(scan.list)
P<- matrix(0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
P[non.diagonal(P)]<- Bin.prob$present



Scan<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))

prob.test<- non.diagonal(P,"vector")
system.time(
  S0<- rbind_lapply(seq_along(scan.list),
                    function(l){
                      table(
                        sapply(1:length(Scan[non.diagonal(Scan)]),
                               function(dyad) {
                                 Scan[non.diagonal(Scan)][dyad]<- sample(c(1,0),1,replace = TRUE,prob=c(prob.test[dyad],1-prob.test[dyad]))
                               }
                        )
                      )
                    }
  )
)
S0<- data.table(S0)
S0[`1`==380]$`1`<- 0
table(S0$`1`>0)


P.cond<- P/(1-prod(1-P))
prob.test<- non.diagonal(P.cond,"vector")
system.time(
  S0.cond<- rbind_lapply(seq_along(list.of.non.zero[list.of.non.zero]),
                  function(l){
                    table(
                      sapply(1:length(Scan[non.diagonal(Scan)]),
                             function(dyad) {
                               Scan[non.diagonal(Scan)][dyad]<- sample(c(1,0),1,replace = TRUE,prob=c(prob.test[dyad],1-prob.test[dyad]))
                             }
                      )
                    )
                  }
  )
)
S0.cond<- data.table(S0.cond)
S0.cond[`1`==380]$`1`<- 0
table(S0.cond$`1`>0)

prob.test

system.time(
  test.rbinom<- rbind_lapply(1:500000,
                       function(k){
                         scan<- rbinom(length(prob.test),1,prob.test)
                         if(is.all.equal(scan)){
                           rep(0,length(prob.test))
                         }
                         else{
                           scan
                         }
                       }
  )
)

system.time(
  test.sample<- sapply(1:500,
                       function(k){
                         sapply(1:length(Scan[non.diagonal(Scan)]),
                                function(dyad) {
                                  Scan[non.diagonal(Scan)][dyad]<- sample(c(1,0),1,replace = TRUE,prob=c(prob.test[dyad],1-prob.test[dyad]))
                                }
                         )
                       }
  )
)
test.rbinom<- t(test.rbinom)
test.sample<- t(test.sample)

test.rbinom<- data.table(test.rbinom)
test.sample<- data.table(test.sample)

test.rbinom$method<- "rbinom"
test.sample$method<- "sample"

test.dt<- rbind(test.rbinom,test.sample)
test.long<- melt.data.table(test.dt,381,1:380,variable.name = "dyad",value.name = "a")
test.long$a.scaled<- as.vector(scale(test.long$a))
test.summary<- test.long[,.(a=mean(a.scaled),sd=sd(a.scaled)),by = .(method,dyad)]
ggplot(test.summary[1:100,],aes(dyad,a,colour=method,group=method))+
  geom_line(position = position_dodge(.75))+
  geom_linerange(aes(ymin = a-sd/20,ymax=a+sd/20),position = position_dodge(.75),alpha=0.8)+
  geom_point(alpha=1,position = position_dodge(.75))+theme_bw()
ggplot(test.summary[1:200,],aes(dyad,a,colour=method,group=method))+
  geom_line()+
  geom_linerange(aes(ymin = a-sd/10,ymax=a+sd/10),alpha=0.8)+
  geom_point(alpha=1)+theme_bw()

# draft toy example from Schieber et al (2007) ----------------------------
#
# N1<- matrix(c(0,0,0,1,0,0,0,0,0,
#               0,0,0,0,1,0,0,0,0,
#               0,0,0,1,0,0,0,0,0,
#               1,0,1,0,1,0,1,0,0,
#               0,1,0,1,0,1,1,0,0,
#               0,0,0,0,1,0,0,0,0,
#               0,0,0,1,1,0,0,1,1,
#               0,0,0,0,0,0,1,0,0,
#               0,0,0,0,0,0,1,0,0),nrow = 9,byrow = TRUE,dimnames = list(as.character(1:9),as.character(1:9)))
# N1
# plot(igraph::graph.adjacency(N1,mode = "max"))
#
#
#
# N2<- matrix(c(0,0,0,1,0,0,0,0,0,
#               0,0,0,0,0,0,0,0,0,
#               0,0,0,1,0,0,0,1,0,
#               1,0,1,0,0,0,0,0,0,
#               0,0,0,0,0,1,1,0,0,
#               0,0,0,0,1,0,0,0,1,
#               0,0,0,0,1,0,0,1,1,
#               0,0,1,0,0,0,1,0,1,
#               0,0,0,0,0,1,1,1,0),nrow = 9,byrow = TRUE,dimnames = list(as.character(1:9),as.character(1:9)))
# N2
# plot(igraph::graph.adjacency(N2,mode = "max"))
#
# N3<- matrix(c(0,0,1,1,0,0,0,0,0,
#               0,0,0,0,1,1,0,0,0,
#               1,0,0,1,0,0,0,0,0,
#               1,0,1,0,0,0,0,0,0,
#               0,1,0,0,0,1,0,0,0,
#               0,1,0,0,1,0,0,0,0,
#               0,0,0,0,0,0,0,1,1,
#               0,0,0,0,0,0,1,0,1,
#               0,0,0,0,0,0,1,1,0),nrow = 9,byrow = TRUE,dimnames = list(as.character(1:9),as.character(1:9)))
# N3
# plot(igraph::graph.adjacency(N3,mode = "max"))
#
# Frobenius_from_adjacency(N1,N2)
# Frobenius_from_adjacency(N1,N3)
# Frobenius_from_adjacency(N2,N3)
#
# Ken_gw(N1,N2,"max")
# Ken_gw(N1,N3,"max")
# Ken_gw(N3,N2,"max")
#
# cor(N1[upper.tri(N1)],N2[upper.tri(N2)])
# cor(N1[upper.tri(N1)],N3[upper.tri(N3)])
# cor(N1[upper.tri(N2)],N3[upper.tri(N3)])
#
#
# g<- igraph::graph.adjacency(Adj,mode = "max",weighted = TRUE)
# g<- igraph::graph.adjacency(N3,mode = "max",weighted = TRUE)
# h<- igraph::graph.adjacency(N1,mode = "max",weighted = TRUE)
# plot(igraph::graph.adjacency(Adj,mode = "max",weighted = TRUE))
# igraph::shortest.paths(g,algorithm = "di")
# m<- igraph::shortest.paths(igraph::graph.adjacency(N3,mode = "max",weighted = TRUE),algorithm = "di")
# m[which(m=="Inf")]<- 9
# setdiff(m,0)
#



# draft on distance -------------------------------------------------------------------

# Bootstrap.list<- Boot_scans(Adj,1,total_scan = total_scan,obs.prob = 0.5,
#                             focal.list = focal.list,scaled = TRUE,keep=TRUE,
#                             method = "both",mode = "directed",output = "adj")
# theo<- Boot_get.list(Bootstrap.list,what = "theo")[[1]]
# group<- Boot_get.list(Bootstrap.list,what = "obs")[[1]]
# focal<- Boot_get.list(Bootstrap.list,what = "foc")[[1]]
#
# g.Adj<- function(Adj,mode="max"){
#   igraph::graph.adjacency(Adj,mode,weighted = TRUE)
# }
#
#
# plot(g.Adj(Adj,"max"))
# plot(g.Adj(theo,"max"))
# plot(g.Adj(group,"max"))
# plot(g.Adj(focal,"max"))
#
# theo<- adjacency_mode(theo,"max")
# group<- adjacency_mode(group,"max")
# focal<- adjacency_mode(focal,"max")
#
#
# plot(non.diagonal(Adj,"vec"),non.diagonal(theo,"vec"))
# plot(non.diagonal(theo,"vec"),non.diagonal(group,"vec"))
# plot(non.diagonal(theo,"vec"),non.diagonal(focal,"vec"))
#
# igraph::shortest.paths(g.Adj(theo),algorithm = "dijkstra")
#
# lapply(seq_len(nrow(theo)),
#        function(i){
#          hist(unique(theo[i,]))
#        }
# )
#
#
# hist(non.diagonal(theo,"vec"))
# hist(non.diagonal(group,"vec"))
#
# # df<- rbind(data.frame(method="theo",weight = non.diagonal(theo,"vec")),
# #            data.frame(method="group",weight = non.diagonal(group,"vec")),
# #            data.frame(method="focal",weight = non.diagonal(focal,"vec")))
# # library(data.table)
# # df<- data.table(df)
# # df$weight.scaled<- scale(df$weight)
# # df[,weight.scaled.relative := scale(weight),by = .(method)]
# # hist(df$weight.scaled)
# # library(ggplot2)
# # ggplot(df,aes(method,weight,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
# # ggplot(df,aes(method,weight.scaled,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
# # ggplot(df,aes(method,weight.scaled.relative,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
#
#
# theo.dist<- non.diagonal(igraph::shortest.paths(g.Adj(theo),algorithm = "dijkstra"),"vec")
# group.dist<- non.diagonal(igraph::shortest.paths(g.Adj(group),algorithm = "dijkstra"),"vec")
# focal.dist<- non.diagonal(igraph::shortest.paths(g.Adj(focal),algorithm = "dijkstra"),"vec")
#
# scaled.inf<- function(X.dist){
#   ifelse(!is.infinite(X.dist),X.dist/sd(X.dist[!is.infinite(X.dist)]),1+max(X.dist[!is.infinite(X.dist)])/sd(X.dist[!is.infinite(X.dist)]))
# }
#
# scaled.inf<- function(X.dist){
#   ifelse(!is.infinite(X.dist),scale(X.dist[!is.infinite(X.dist)]),1+max(scale(X.dist[!is.infinite(X.dist)])))
# }
#
# scaled.inf(theo.dist)
# scaled.inf(group.dist)
# scaled.inf(focal.dist)
#
# sd(focal.dist[!is.infinite(focal.dist)])
#
# test.distances<- rbind(data.frame(method="theo",dist = scaled.inf(theo.dist)),
#                   data.frame(method="group",dist = scaled.inf(theo.dist)),
#                   data.frame(method="focal",dist = scaled.inf(theo.dist)))
# distances<- data.table(test.distances)
# ggplot(distances,aes(method,dist,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
#
#
# across.max<- max(scaled.inf(theo.dist),scaled.inf(group.dist),scaled.inf(focal.dist))
# distances<- rbind(data.frame(method="theo",dist = ifelse(!is.infinite(theo.dist),theo.dist,across.max)),
#                   data.frame(method="group",dist = ifelse(!is.infinite(group.dist),group.dist,across.max)),
#                   data.frame(method="focal",dist = ifelse(!is.infinite(focal.dist),focal.dist,across.max)))
# library(data.table)
# distances<- data.table(distances)
# distances$dist.scaled<- scale(distances$dist,center = FALSE)
# distances[,dist.scaled.relative := scale(dist,center = FALSE),by = .(method)]
# distances[,dist.scaled.max := dist/max(dist),by = .(method)]
# distances[,dist.scaled.min := dist/min(dist),by = .(method)]
# distances[,dist.scaled.range := dist/(max(dist)-min(dist)),by = .(method)]
# hist(distances$dist.scaled)
# library(ggplot2)
# ggplot(distances[dist<5],aes(method,dist,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
# ggplot(distances[method!="focal"&dist!=across.max],aes(method,dist,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
# ggplot(distances[dist!=n],aes(method,dist.scaled,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
# ggplot(distances[dist!=across.max],aes(method,dist.scaled.relative,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
# ggplot(distances[dist!=across.max],aes(method,dist.scaled.max,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
# ggplot(distances[dist!=across.max],aes(method,dist.scaled.min,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
# ggplot(distances[dist<10],aes(method,dist.scaled.range,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
#
# scaled.range<- range(distances$dist.scaled.max)
# 1/min(distances$dist.scaled.max)
# seq(0,1,by = min(distances$dist.scaled.max))
# length(unique(distances[method=="focal"]$dist.scaled.max))
# distances$floored<- round(distances$dist.scaled.max/min(distances[method=="theo"]$dist.scaled.max))
# ggplot(distances,aes(method,floored,colour=method))+geom_jitter()+geom_boxplot(alpha=0.8)+theme_bw()
#
#
#
# max(round(focal.dist[!is.infinite(focal.dist)]/min(focal.dist[!is.infinite(focal.dist)])))+1
#
#
# plot(non.diagonal(theo,"vec"),c(theo.dist))



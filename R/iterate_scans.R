#' Iterate scans
#' Internal use in Boot_scan. Iterate several binary group or focal scans with probabilities derived from a provided adjancecy matrix, to produce a new adjancecy matrix.
#'
#' @param Adj square integers matrix of occurences of dyads. Optional if using presence.prob. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param presence.prob. square probability matrix of presence (as in Bernouilli trial) of dyads. Optional if using Adj and total_scan.
#' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
#' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param n.cores number of threads to use while performingh the bootstrap
#' @param cl Optional cluster object (cf snow package), experimentally set to put the makeCluster and stopCluster out of the bootable function. (WIP, next implementation should rethink this).
#' @param ... additional argument to be used, to use produce a scan in a desired way.
#' @param focal.list Character vector, indicate the list of focals to consider throughout the scans.
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument
#'
#' @return according to output and method: a list of iterated scans, or an adjacency matrix
#' @export
#' @importFrom parallel detectCores
#' @importFrom snow makeCluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom pbapply pblapply
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- as.character(1:n);total_scan<- 42;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
#' Adj
#'
#' presence.prob<- Binary.prob(Adj,50)
#' obs.prob<- matrix(runif(n*n,0,1),n,n);diag(obs.prob)<- 0
#'
#' focal.list<- sample(nodes,total_scan,replace = TRUE)
#' table(focal.list)
#'
#' iterate_scans(Adj,total_scan,scaled = FALSE,method = "both",output = "list",obs.prob = 0.8,n.cores = 1)
#' iterate_scans(Adj,total_scan,scaled = TRUE,method = "group",output = "adjacency",obs.prob = obs.prob,n.cores = 1)
#' iterate_scans(Adj,total_scan,scaled = FALSE,method = "focal",output = "adjacency",n.cores = 1)
#' iterate_scans(Adj,total_scan,focal.list = focal.list,scaled = TRUE,obs.prob = 0.7,method = "both",mode = "directed",output = "all")

iterate_scans<- function(Adj=NULL,total_scan,method=c("theoretical","group","focal","both"),
                         output=c("list","adjacency","all"),scaled=FALSE,...,
                         n.cores=(parallel::detectCores()-1),cl){
  scan.default.args(Adj = Adj,total_scan = total_scan,method = method,...)

  if(missing(cl)) {
    .export<- c("do.scan","non.diagonal","quick.sample","Binary.prob","binary_adjacency_mode","scan.default.args","observable_edges","focal.scan");
    cl<- snow::makeCluster(n.cores);
    snow::clusterExport(cl,list = .export);snow::clusterExport(cl,list = ls(sys.frame(which = 1)),envir = sys.frame(which = 1));
    doSNOW::registerDoSNOW(cl);on.exit(snow::stopCluster(cl))
  }

  switch(method,
         "theoretical" = ,
         "group" = {
           scan_list<- pbapply::pblapply(
             1:total_scan,
             function(i){
               do.scan(presence.prob = presence.prob,method = method,obs.prob = obs.prob,Adj.subfun = Adj.subfun)
             },cl = cl
           )
         },
         "focal" = {
           scan_list<- pbapply::pblapply(
             1:total_scan,
             function(i){
               do.scan(presence.prob = presence.prob,method = "focal",focal = focal.list[i],Adj.subfun = Adj.subfun)
             },cl = cl
           )
         },
         "both" = {
           scan_list<- pbapply::pblapply(
             1:total_scan,
             function(i){
               do.scan(presence.prob = presence.prob,method = "both",obs.prob = obs.prob,focal = focal.list[i],Adj.subfun = Adj.subfun)
             },cl = cl
           )
         }
  )

  switch(output,
         "list" = scan_list,
         "adjacency" = sum_up.scans(scan_list = scan_list,scaled = scaled,method = method,mode = mode),
         "all" = list(
           list = scan_list,
           adjacency = sum_up.scans(scan_list = scan_list,scaled = scaled,method = method,mode = mode)
         )
  )
}


#' Title
#'
#' @param Adj
#'
#' @param total_scan
#' @param method
#' @param output
#' @param scaled
#' @param ...
#' @param n.cores
#' @param cl
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- as.character(1:n);total_scan<- 200;
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:total_scan/100,n*(n-1),replace = TRUE)
#' Adj
#'
#' presence.prob<- Binary.prob(Adj,total_scan)
#' obs.prob<- matrix(runif(n*n,0,1),n,n);diag(obs.prob)<- 0
#'
#' focal.list<- sample(nodes,total_scan,replace = TRUE)
#' table(focal.list)
#'
#' iterate_rare.scans(Adj,total_scan,scaled = FALSE,method = "both",output = "list",obs.prob = 0.8,n.cores = 1)
#' iterate_rare.scans(Adj,total_scan,scaled = TRUE,method = "group",output = "adjacency",obs.prob = obs.prob,n.cores = 1)
#' iterate_rare.scans(Adj,total_scan,scaled = FALSE,method = "focal",output = "adjacency",n.cores = 1)
#' iterate_rare.scans(Adj,total_scan,focal.list = focal.list,scaled = TRUE,obs.prob = 0.7,method = "both",mode = "directed",output = "all")

iterate_rare.scans<- function(Adj=NULL,total_scan,method=c("theoretical","group","focal","both"),
                              output=c("list","adjacency","all"),scaled=FALSE,...,
                              n.cores=(parallel::detectCores()-1),cl){
  scan.default.args(Adj = Adj,total_scan = total_scan,method = method,...)

  if(missing(cl)) {
    .export<- c("do.non.zero.scan","non.diagonal","quick.sample","Binary.prob","binary_adjacency_mode","scan.default.args","observable_edges","adjust.conditional.prob","focal.scan");
    cl<- snow::makeCluster(n.cores);
    snow::clusterExport(cl,list = .export);snow::clusterExport(cl,list = ls(sys.frame(which = 1)),envir = sys.frame(which = 1));
    doSNOW::registerDoSNOW(cl);on.exit(snow::stopCluster(cl))
  }

  scan_list<- simulate.zeroes.non.zeroes(total_scan,presence.prob,method = method)
  non.zeroes<- sapply(scan_list,is.null)

  switch(method,
         "theoretical" = ,
         "group" = {
           scan_list[non.zeroes]<- pbapply::pblapply(
             seq_along(scan_list[non.zeroes]),
             function(i){
               do.non.zero.scan(presence.prob = presence.prob,method = method,obs.prob = obs.prob,Adj.subfun = Adj.subfun)
             },cl = cl
           )
         },
         "focal" = {
           scan_list[non.zeroes]<- pbapply::pblapply(
             seq_along(scan_list[non.zeroes]),
             function(i){
               do.non.zero.scan(presence.prob = presence.prob,method = "focal",focal = focal.list[i],Adj.subfun = Adj.subfun)
             },cl = cl
           )
         },
         "both" = {
           scan_list[non.zeroes]<- pbapply::pblapply(
             seq_along(scan_list[non.zeroes]),
             function(i){
               do.non.zero.scan(presence.prob = presence.prob,method = "both",obs.prob = obs.prob,focal = focal.list[i],Adj.subfun = Adj.subfun)
             },cl = cl
           )
         }
  )

  switch(output,
         "list" = scan_list,
         "adjacency" = sum_up.scans(scan_list = scan_list,scaled = scaled,method = method,mode = mode),
         "all" = list(
           list = scan_list,
           adjacency = sum_up.scans(scan_list = scan_list,scaled = scaled,method = method,mode = mode)
         )
  )
}

iterate_rare.scans(Adj,total_scan,scaled = TRUE,method = "group",output = "all",obs.prob = obs.prob,n.cores = 1)


bench<- microbenchmark::microbenchmark(standard = iterate_scans(Adj,total_scan,scaled = TRUE,method = "both",output = "adjacency",obs.prob = obs.prob,n.cores = 1),
                                       rare = iterate_rare.scans(Adj,total_scan,scaled = TRUE,method = "both",output = "adj",obs.prob = obs.prob,n.cores = 1),
                                       unit = "ms")
summary(bench)
ggplot2::autoplot(bench)

cl<- snow::makeCluster(7);
snow::clusterExport(cl,list = ls(.GlobalEnv))

test.dt<- rbind_lapply(1:100,
       function(b){
         time<- microbenchmark::microbenchmark(
           standard = {standard<- iterate_scans(Adj,2000,scaled = FALSE,method = "both",output = "adjacency",obs.prob = obs.prob,n.cores = 7)},
           rare = {rare<- iterate_rare.scans(Adj,2000,scaled = FALSE,method = "both",output = "adjacency",obs.prob = obs.prob,n.cores = 7)},
           unit = "ms",times = 1
         )
         cor<- data.table::data.table(
           method=c("standard","rare"),
           cor=c(stats::cor(Adj,standard$theoretical),
                 stats::cor(Adj,rare$theoretical))
         )
         cbind(cor,time=summary(time)[["median"]],boot=b)
       }
)

ggplot(test.dt,aes(method,cor,colour=method))+geom_jitter(alpha=0.2)+geom_violin(alpha=0.75)+theme_bw()
ggplot(test.dt,aes(method,time,colour=method))+geom_jitter(alpha=0.2)+geom_violin(alpha=0.75)+theme_bw()

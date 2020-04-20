#' Perform a scan
#' According to edge presence probability, or a reference adjacency matrix `Adj` coupled with a sample effort `total_scan`, perform a scan in the form of a theoretical binary adjacency matrix, and simulate the empirical obtention of a group scan and/or focal scan.
#'
#' @param Adj square integers matrix of occurences of dyads. Optional if using presence.prob. Update: now if presence prob is passed as Adj (thus all(Adj<1) is TRUE), it will be rightfully assigned to presence.prob. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
#' @param ... additional argument to be used, to use produce a scan in a desired way.
#' @param focal Only required for method = "focal" or "both" Character scalar, indicate which focal to consider for the scan.
#' @param presence.prob square probability matrix of presence (as in Bernouilli trial) of dyads. Optional if using Adj and total_scan.
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument
#'
#' @return a list of the theoretical square binary matrix representing the whole group scan, and:
#' \itemize{
#'  \item{"method = theoretical"}{that's all,}
#'  \item{"method = group"}{the empirical group scan square binary matrix}
#'  \item{"method = focal"}{the empirical focal scan binary vector (and the focal identity as a "focal" attribute)}
#'  \item{"method = both"}{or the three above}
#'  }
#'
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
#' Adj
#'
#' presence.prob<- Binary.prob(Adj,50)
#'
#' do.scan(Adj,50,method = "group")
#' do.scan(Adj,50,method = "both")
#' do.scan(presence.prob,method = "focal")

do.scan<-function(Adj=NULL,total_scan=NULL,
                  method = c("theoretical","group","focal","both"),...,check.defaults=TRUE){
  method<- match.arg(method)
  if(check.defaults){scan.default.args(Adj,total_scan,method,...)}

  if(!is.null(Adj)){n<- nrow(Adj);nodes_names<- rownames(Adj)} else {n<- nrow(presence.prob);nodes_names<- rownames(presence.prob)}
  presence.P<- Adj.subfun(presence.prob,"vector");p<- length(presence.P)

  scan<- matrix(0,nrow = n,ncol = n,dimnames = list(nodes_names,nodes_names))
  scan[Adj.subfun(scan)]<- rbinom(p,1,presence.P)
  observable_edges<- function(scan,obs.prob) {scan}
  switch(method,
         "theoretical" = scan,
         "group" = list(theoretical = scan,
                        group = observable_edges(scan,obs.prob)
         ),
         "focal" = list(theoretical = scan,
                        focal = {
                          focal.scan<- scan[focal,];
                          attr(focal.scan,"focal")<- focal;
                          focal.scan
                        }
         ),
         "both" = list(theoretical = scan,
                       group = observable_edges(scan,obs.prob),
                       focal = { # should write a symmetrical to observable_edges() for the focal case
                         focal.scan<- scan[focal,];
                         attr(focal.scan,"focal")<- focal;
                         focal.scan
                       }
         )
  )
}

#' Perform a scan knowing that at least one edge has been observed
#' According to edge presence probability, or a reference adjacency matrix `Adj` coupled with a sample effort `total_scan`, perform a scan in the form of a theoretical binary adjacency matrix, and simulate the empirical obtention of a group scan and/or focal scan. Involves calculation of approximate conditional probabilities. Required for the optimization for rare events, not fully tested in other cases.
#'
#' @param Adj square integers matrix of occurences of dyads. Optional if using presence.prob. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
#' @param ... additional argument to be used, to use produce a scan in a desired way.
#' @param focal Only required for method = "focal" or "both" Character scalar, indicate which focal to consider for the scan.
#' @param presence.prob. square probability matrix of presence (as in Bernouilli trial) of dyads. Optional if using Adj and total_scan.
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument
#'
#' @return a list of the theoretical square binary matrix representing the whole group scan, including at least a one, and:
#' \itemize{
#'  \item{"method = theoretical"}{that's all,}
#'  \item{"method = group"}{the empirical group scan square binary matrix}
#'  \item{"method = focal"}{the empirical focal scan binary vector (and the focal identity as a "focal" attribute)}
#'  \item{"method = both"}{or the three above}
#'  }
#'
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
#' Adj
#' do.non.zero.scan(Adj,500,method = "both",focal=3)
do.non.zero.scan<-function(Adj=NULL,total_scan=NULL,
                           method = c("theoretical","group","focal","both"),...,check.defaults=TRUE){
  method<- match.arg(method)
  if(check.defaults){scan.default.args(Adj,total_scan,method,...)}

  if(!is.null(Adj)){n<- nrow(Adj);nodes_names<- rownames(Adj)} else {n<- nrow(presence.prob);nodes_names<- rownames(presence.prob)}
  presence.P<- Adj.subfun(presence.prob,"vector");p<- length(presence.P)

  scan<- matrix(0,nrow = n,ncol = n,dimnames = list(nodes_names,nodes_names))
  rand.order<- sample(1:p,p)
  P.cond<- cumsum(adjust.conditional.prob(presence.P[rand.order]))
  first.one<- min(which(runif(1)<P.cond))
  scan[Adj.subfun(scan)][rand.order][first.one]<- 1
  scan[Adj.subfun(scan)][rand.order][-first.one]<- rbinom(p-1,1,presence.P[rand.order][-first.one])
  observable_edges<- function(scan) {scan}
  switch(method,
         "theoretical" = scan,
         "group" = list(theoretical = scan,
                        group = observable_edges(scan)
         ),
         "focal" = list(theoretical = scan,
                        focal = {
                          focal.scan<- scan[focal,];
                          attr(focal.scan,"focal")<- focal;
                          focal.scan
                        }
         ),
         "both" = list(theoretical = scan,
                       group = observable_edges(scan),
                       focal = {
                         focal.scan<- scan[focal,];
                         attr(focal.scan,"focal")<- focal;
                         focal.scan
                       }
         )
  )
}

#' Set arguments to default for do.(non.zero.)scan() function when necessary
#'
#' @param Adj square integers matrix of occurences of dyads. Optional if using presence.prob. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
#' @param ... additional argument to be used, to use produce a scan in a desired way.
#' @param focal Only required for method = "focal" or "both" Character scalar, indicate which focal to consider for the scan. Default is a random node.
#' @param presence.prob. square probability matrix of presence (as in Bernouilli trial) of dyads. Optional if using Adj and total_scan.
#' @param obs.prob Default behaviour is obs.prob=1. Otherwise either:
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed.}
#' }
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument. Default is non.diagonal if mode is defaulted to directed.
#'
#' @return nothing, but assign required variable for do.(non.zero.)scan() in their environment.
#'
#' @export
#'
#' @examples
#' # Internal use in do.(non.zero.)scan()
scan.default.args<- function(Adj,total_scan,method,...){
  opt.args<- list(...)

  lapply(names(opt.args),function(name) assign(name,opt.args[[name]],parent.frame(n = 2)))

  if(is.null(opt.args$obs.prob)){assign("obs.prob",1,parent.frame())}else{assign("obs.prob",opt.args$obs.prob,parent.frame(n = 1))}

  if(is.null(opt.args$mode)){assign("mode","directed",parent.frame());opt.args$mode<- "directed"}else{assign("mode",opt.args$mode,parent.frame(n = 1))}

  if(is.null(opt.args$Adj.subfun)){
    assign(
      "Adj.subfun",
      switch(opt.args$mode,
             "directed" = ,
             "undirected" = ,
             "max" = ,
             "min" = ,
             "plus" = non.diagonal,
             "upper" = upper.tri,
             "lower" =  lower.tri,
             "vector" = function(V){rep(TRUE,length(V))}
      ),
      parent.frame(n = 1)
    )
  }else{assign("Adj.subfun",opt.args$Adj.subfun,parent.frame(n = 1))}

  if(is.null(opt.args$presence.prob)){
    if(all(Adj<1&Adj>=0)){
      assign("presence.prob",Adj,parent.frame(n = 1)) # assuming presence prob has been passed as first argument `Adj`
    }
    else{
      if(is.null(total_scan)){
        stop("Please provide either `total_scan` and `Adj`, or `presence.prob` to run a scan.")
      }
      assign("presence.prob",Binary.prob(Adj,total_scan,mode = opt.args$mode),parent.frame(n = 1))
    }
  }else{assign("presence.prob",opt.args$presence.prob,parent.frame(n = 1))}

  if(is.null(opt.args$focal)){
    if(method!="group"){
      if(!is.null(Adj)){
        assign("focal",quick.sample(1:nrow(Adj),1),parent.frame(n = 1))
      }else{
        assign("focal",quick.sample(1:nrow(presence.prob),1),parent.frame(n = 1))
      }
    }else{
      assign("focal",NULL,parent.frame(n = 1))
    }
  }else{assign("focal",opt.args$focal,parent.frame(n = 1))}

  if(is.null(opt.args$focal.list)&!is.null(opt.args$total_scan)){
    if(method!="group"){
      if(!is.null(Adj)){
        assign("focal.list",quick.sample(1:nrow(Adj),total_scan),parent.frame(n = 1))
      }else{
        assign("focal.list",quick.sample(1:nrow(presence.prob),total_scan),parent.frame(n = 1))
      }
    }else{
      assign("focal.list",NULL,parent.frame(n = 1))
    }
  }else{assign("focal.list",opt.args$focal.list,parent.frame(n = 1))}
}

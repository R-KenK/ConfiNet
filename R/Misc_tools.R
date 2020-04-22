#' Row bind list of data frames
#' wrapper to one-function do.call rbind over a lapply list
#'
#' @param X a list. See details \link[base]{lapply}.
#' @param FUN a function to subset data frames (or data tables). See details \link[base]{lapply}.
#'
#' @return a row bound data frame
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' X<- lapply(1:3,function(i) list(int = 42,df = data.frame(x = runif(10,0,1),y = runif(10,0,1))))
#' rbind_lapply(X,function(x) x$df)
rbind_lapply<- function(X,FUN){
  do.call(rbind,lapply(X = X,FUN = FUN))
}

#' Row bind list of data frames
#' wrapper to one-function do.call rbind over a pblapply list
#'
#' @param X a list. See details \link[base]{lapply}.
#' @param FUN a function to subset data frames (or data tables). See details \link[base]{lapply}.
#'
#' @return a row bound data frame
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' X<- lapply(1:3,function(i) list(int = 42,df = data.frame(x = runif(10,0,1),y = runif(10,0,1))))
#' rbind_lapply(X,function(x) x$df)
rbind_pblapply<- function(X,FUN,n.cores,.export=NULL){
  cl<- make_cl(n.cores,.export);on.exit(snow::stopCluster(cl))
  do.call(rbind,pbapply::pblapply(X = X,FUN = FUN,cl = cl))
}

#' Make cluster object wrapper
#'
#' @param n.cores number of threads to use.
#' @param .export vector of variable/function names to use in each clusters
#'
#' @return a cluster object
#' @export
#'
#' @importFrom  snow makeCluster
#' @importFrom  snow clusterExport
#' @importFrom  doSNOW registerDoSNOW
#'
#' @examples
#' # Internal use
make_cl<- function(n.cores,.export){
  cl<- snow::makeCluster(n.cores);snow::clusterExport(cl,list = .export);doSNOW::registerDoSNOW(cl);
  cl
}

#' Quick optimized equivalent to sample(x,size,replace=TRUE)
#'
#' @param x a vector
#' @param size number of elements to sample
#'
#' @return sample of size "size" taken from x
#' @export
#'
#' @examples
#' quick.sample(1:20,5)
#' # microbenchmark::microbenchmark(runif={(1:20)[ceiling(runif(5,0,20))]},quick.sample=quick.sample(1:20,5),sample=sample(1:20,5,replace = TRUE),times = 1000,control = list("warmup"=100))
quick.sample<- function(x,size){
  x[ceiling(runif(size,0,length(x)))]
}

# Rare events optimization related functions ------------------------------

#' Placeholder for calculating expected time using the optimization for rare event
#'
#' @param n number of node of the original network
#' @param total_scan sampling effort
#' @param max.obs maximum value of the weighted adjacency matrix of the original network
#'
#' @return a time value to be compared with the one using the standard approach
#' @export
#'
#' @examples
#' # Internal use in decide_use.rare.opti
opti.expected.time<- function(n,total_scan,max.obs){
  # try and model?
}

#' Placeholder for calculating expected time using the standard scan method
#' i.e. performing all scans
#'
#' @param n number of node of the original network
#' @param total_scan sampling effort
#' @param max.obs maximum value of the weighted adjacency matrix of the original network
#'
#' @return a time value to be compared with the one using the optimization for rare event
#' @export
#'
#' @examples
#' # Internal use in decide_use.rare.opti
standard.expected.time<- function(n,total_scan,max.obs){
  # try and model?
}

#' Decide based on expected times if the otpimization for rare event should be used
#'
#' @param n number of node of the original network
#' @param total_scan sampling effort
#' @param max.obs maximum value of the weighted adjacency matrix of the original network
#'
#' @return a logical value meaning that the optimization for rare events should be use when TRUE is returned
#' @export
#'
#' @examples
#' # Internal use
decide_use.rare.opti<- function(n,total_scan,max.obs=NULL){
  if(is.null(max.obs)&is.matrix(n)){max.obs<- max(n)}
  if(is.matrix(n)){n<- nrow(n)}
  # Figure out how to relate n, N and max.obs through opti.expected.time(n,total_scan,max.obs)
  opti.expected.time(n,total_scan,max.obs)<standard.expected.time(n,total_scan,max.obs)
}

#' Simulate which scan returns an all-zeroes matrix
#'
#' @param total_scan integer, sampling effort
#' @param presence.prob presence probability matrix (or vector)
#'
#' @return a list of zero-matrices (all-zeroes scans) and NULL (non-zeroes scans to be later performed)
#' @export
#'
#' @examples
#' # Internal use
simulate.zeroes.non.zeroes<- function(total_scan,presence.prob){
  scan_list<- vector(mode="list",length = total_scan)
  non.zeroes<- rbinom(total_scan,1,1-prod(1-presence.prob))==1;
  scan_list[!non.zeroes]<- lapply(seq_along(scan_list[!non.zeroes]),function(scan) matrix(0,nrow(presence.prob),ncol(presence.prob)));
  scan_list
}

#' Determine step-by-step conditional probabilities for non-zeroes scans
#' Internal use. Returns the CDF of the probability that the i-th dyad (with probability presence.prob[i]) is the first to yield a 1.
#'
#' @param presence.prob presence probability matrix (or vector)
#'
#' @return a cumulative distribution function of the probability that the i-th dyad (with probability presence.prob[i]) is the first to yield a 1
#' @export
#'
#' @details Workflow is as follows: first simulate.zeroes.non.zeroes() determines which scans are all-zeroes and which are non-zeroes. Then for non zeroes, at a random order each dyad is drawn in order with conditional probability that: (1) there is at least one 1 in the scan, and (2) all the previous coins were zeroes. Once the first one is drawn, the rest are drawn with their regular probabilities. In details, cumulative density probability of each dyad (in a given random order) to be the first one to be a 1 is calculated, and a random draw determine which one is first, set the previous ones to zero, and draw the rest normally. cf. Supplmenentary material X.
#'
#' @examples
#' # Internal use.
adjust.conditional.prob<- function(presence.prob){
  prob.all.zeroes<- 1-prod(1-presence.prob)
  previous.are.zeros<- c(1,cumprod(1-presence.prob[1:(length(presence.prob)-1)]))
  sapply(1:length(presence.prob),function(i) presence.prob[i]/prob.all.zeroes*previous.are.zeros[i])
}

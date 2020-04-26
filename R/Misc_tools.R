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
#' @param n.cores number or cores (threads rather) to use in the cluster.
#' @param .export character vector of the name of the variable and function to export to each worker of the cluster
#' @param cl cluster object
#'
#' @return a row bound data frame
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' X<- lapply(1:3,function(i) list(int = 42,df = data.frame(x = runif(10,0,1),y = runif(10,0,1))))
#' rbind_lapply(X,function(x) x$df)
rbind_pblapply<- function(X,FUN,n.cores=NULL,.export=NULL,cl=NULL){
  if(is.null(cl)){cl<- make_cl(n.cores,.export);on.exit(snow::stopCluster(cl))}
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
  cl<- snow::makeCluster(n.cores);snow::clusterExport(cl,list = .export);
  cl
}

#' Quick optimized equivalent to sample(x,size,replace=TRUE)
#'
#' @param x a vector
#' @param size number of elements to sample
#'
#' @importFrom stats runif
#'
#' @return sample of size "size" taken from x
#' @export
#'
#' @examples
#' quick.sample(1:20,5)
#' # microbenchmark::microbenchmark(runif={(1:20)[ceiling(runif(5,0,20))]},
#' #   quick.sample=quick.sample(1:20,5),sample=sample(1:20,5,replace = TRUE),
#' #   times = 1000,control = list("warmup"=100))
quick.sample<- function(x,size){
  x[ceiling(stats::runif(size,0,length(x)))]
}

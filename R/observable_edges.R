#' Hide unobservable edges
#' Simulate that some dyads might reasonable not be observable during a group scan
#'
#' @param Scan a square matrix
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation obs.probability matrix"}{of same dimension as Scan}
#'  \item{"a dyad observation vector"}{subsetted similarly as Scan (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument
#'
#' @return a similar matrix as Scan where some non diagonal edges have a obs.probability to be NAs.
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 6;nodes<- as.character(1:n);
#' Scan<- matrix(rbinom(n*n,1,0.2),n,n,dimnames = list(nodes,nodes));diag(Scan)<- 0;
#' obs.prob<- matrix(runif(n*n,0,1),n,n);diag(obs.prob)<- 0
#'
#' observable_edges(Scan,obs.prob,non.diagonal)
observable_edges<- function(Scan,obs.prob=NULL,Adj.subfun=NULL){
  observed<- Scan
  if(!is.null(obs.prob)){
    if(length(obs.prob)==1) {
      if(obs.prob<=1 & obs.prob>=0){
        obs.prob<- rep(obs.prob,length(Adj.subfun(Scan,"vector")))
      }else{
        stop("Single observation obs.probability provided should be within [0,1]")
      }
    }else{
      if(is.matrix(obs.prob)) {
        obs.prob<- Adj.subfun(obs.prob,"vector")
      }
      if(length(obs.prob)!=length(Adj.subfun(Scan,"vector"))){
        stop("Matrix or vector obs.prob dimension(s) incompatible with adjacency matrix's")
      }
    }
    missed<- rbinom(length(obs.prob),1,obs.prob)==0
    observed[Adj.subfun(observed)][missed]<- NA
  }
  observed
}

#' Get number of edge observations (for group scans with unobservable individuals)
#' quantify actual edge-wise sampling effort, considering that some weren't observable in all group scans.Internal use.
#'
#' @param scan_list list of binary group scans, with NAs when the dyad was not observable.
#' @param diag integer (mostly), value to replace the diagonal of the output matrix with. Use NULL if you consider self-loop (untested).
#'
#' @return a square matrix with element quantifying how many time a dyad has been sampled
#' @export
#'
#' @examples
#' #internal use.
n.observed_edges<- function(scan_list,diag=0,use.rare.opti=FALSE,obs.prob=NULL,n.zeros = NULL){
  n.observed<- Reduce("+",
                      lapply(scan_list,
                             function(scan){
                               observed<- ifelse(!is.na(scan),1,0) # counting part of the algorithm
                               # if(!is.null(diag)) {diag(observed)<- diag} # doesn't count the diagonal by default. Left the option to count if self loops should be considered
                               observed
                             }
                      )
  )
  if(!use.rare.opti){
    n.observed
  }else{
    n<- nrow(n.observed);
    if(!is.matrix(obs.prob)){obs.prob<- matrix(obs.prob,n,n)}
    rbind_lapply(1:n,function(i) rbinom(n,n.zeros,obs.prob[i,])+n.observed[i,])
  }
}

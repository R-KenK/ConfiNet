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

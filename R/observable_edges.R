#' Hide unobservable edges
#' Simulate that some dyads might reasonable not be observable during a group scan
#'
#' @param Scan a square matrix
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation obs.probability matrix"}{of same dimension as Scan}
#'  \item{"a dyad observation vector"}{subsetted similarly as Scan (through the non.diagonal() function for instance)}
#'  \item{"a general dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param keep logical. Indicate if the original "theoretical" scan should be kept track of.
#'
#' @return a similar matrix as Scan where some non diagonal edges have a obs.probability to be NAs. If keep is TRUE, returns a list with theoretical and observed scan.
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 6;
#' Scan<- matrix(sample(c(1,0),6*6,replace=TRUE),n,n);diag(Scan)<- 0
#' obs.prob<- matrix(runif(6*6,0,1),n,n);diag(obs.prob)<- 0
#' obs.prob.single<- runif(1,0,1)
#'
#' observable_edges(Scan,obs.prob)
#' observable_edges(Scan,obs.prob.single)
#'

observable_edges<- function(Scan,obs.prob=NULL,keep=FALSE){
  if(is.matrix(obs.prob)) {obs.prob<- non.diagonal(obs.prob,"vector")}
  if(length(obs.prob)==1) {ifelse(obs.prob<=1 & obs.prob>=0,obs.prob<- rep(obs.prob,length(non.diagonal(Scan,"vector"))),stop("Single observation obs.probability provided should be within [0,1]"))}

  observable<- sapply(obs.prob,function(p) sample(c(TRUE,FALSE),1,prob = c(p,1-p)))

  observed<- Scan;
  observed[non.diagonal(observed)]<- ifelse(observable,non.diagonal(Scan,"vector"),NA)
  if(keep){
    list("theoretical" = Scan,
         "observed" = observed)
  }else{
    observed
  }
}

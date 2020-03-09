#' Hide unobservable edges
#' Simulate that some dyads might reasonable not be observable during a group scan
#'
#' @param Scan a square matrix
#' @param prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix"}{of same dimension as Scan}
#'  \item{"a dyad observation vector"}{subsetted similarly as Scan (through the non.diagonal() function for instance)}
#'  \item{"a general dyad observation probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param keep logical. Indicate if the original "theoretical" scan should be kept track of.
#'
#' @return a similar matrix as Scan where some non diagonal edges have a probability to be NAs. If keep is TRUE, returns a list with theoretical and observed scan.
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 6;
#' Scan<- matrix(sample(c(1,0),6*6,replace=TRUE),n,n);diag(Scan)<- 0
#' prob<- matrix(runif(6*6,0,1),n,n);diag(prob)<- 0
#' prob.single<- runif(1,0,1)
#'
#' observable_edges(Scan,prob)
#'

observable_edges<- function(Scan,prob=NULL,keep=FALSE){
  if(is.matrix(prob)) {prob<- non.diagonal(prob,"vector")}
  if(length(prob)==1) {ifelse(prob<1 & prob>=0,prob<- rep(prob,length(non.diagonal(Scan,"vector"))),stop("Single observation probability provided should be within [0,1]"))}

  observable<- sapply(prob,function(p) sample(c(TRUE,FALSE),1,prob = c(p,1-p)))

  observed<- Scan;
  observed[non.diagonal(observed)]<- ifelse(observable,non.diagonal(Scan,"vector"),NA)
  if(keep){
    list("theoretical" = Scan,
         "observed" = observed)
  }else{
    observed
  }
}

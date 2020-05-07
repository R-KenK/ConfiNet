# Hiding edges during group scan ------------------------------------------

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
        obs.prob<- rep(obs.prob,length(Scan[Adj.subfun(Scan)]))
      }else{
        stop("Single observation obs.probability provided should be within [0,1]")
      }
    }else{
      if(is.matrix(obs.prob)) {
        obs.prob<- obs.prob[Adj.subfun(obs.prob)]
      }
      if(length(obs.prob)!=length(Scan[Adj.subfun(Scan)])){
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
#' @param use.rare.opti logical: should the optimization for rare event be used?
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation probability matrix (P.obs)"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a systematic dyad observation (P.obs constant for all i,j)"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param n.zeros integer, the attribute outputed by `simulate_zeros.non.zeros`, representing the number of full-zero scans. Used only when use.rare.opti=TRUE
#'
#' @importFrom stats rbinom
#'
#' @return a square matrix with element quantifying how many time a dyad has been sampled
#' @export
#'
#' @examples
#' #internal use.
n.observed_edges<- function(scan_list,diag=NULL,use.rare.opti=FALSE,obs.prob=NULL,n.zeros = NULL){
  n.observed<- Reduce("+",
                      lapply(scan_list,
                             function(scan){
                               observed<- ifelse(!is.na(scan),1,0) # counting part of the algorithm
                               if(!is.null(diag)) {diag(observed)<- diag} # doesn't count the diagonal by default. Left the option to count if self loops should be considered
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

# obs.prob tools ----------------------------------------------------------

#' Produce matrix of probability of observation from user-defined function
#'
#' @param Adj square integers matrix of occurences of dyads.
#' @param obs.prob_fun either a user-defined function of (i,j) that output a probability of presence for the dyad, or a single value to indicate a constant observation probability
#' @param Adj.subfun subsetting function of the adjacency matrix. Default is non.diagonal.
#'
#' @return a matrix of probability of observation for each dyad (obs.prob)
#' @export
#'
#' @examples
#' set.seed(42)
#' n<- 6;nodes<- as.character(1:n);
#' total_scan<- 20;n.boot<- 5;
#' focal.list<- sample(nodes,total_scan,replace = TRUE)
#'
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- round(runif(n*(n-1),0,total_scan*.50))
#'
#' make_obs.prob(Adj)
#' make_obs.prob(Adj,obs.prob_fun = 0.2)
#' make_obs.prob(Adj,obs.prob_fun = function(i,j) i+j)
#' make_obs.prob(Adj,obs.prob_fun = function(i,j){EVs<- compute.EV(Adj,"directed");EVs[i]*EVs[j]})
make_obs.prob<- function(Adj,obs.prob_fun = NULL,
                         Adj.subfun = non.diagonal){
  if(is.numeric(obs.prob_fun)){
    if(length(obs.prob_fun)==1 & obs.prob_fun>0 & obs.prob_fun<1){
      return(obs.prob_fun)
    }else{
      stop("incompatible numeric obs.prob_fun.")
    }
  }

  n<- nrow(Adj);

  if(is.null(obs.prob_fun)){
    obs.prob<- matrix(runif(n,0,1),n,n,dimnames = list(rownames(Adj),colnames(Adj)))
    diag(obs.prob)<- 0
    return(obs.prob)
  }

  dyads<- expand.grid(row = 1:n,col = 1:n)
  obs.prob<- matrix(nrow = n,ncol = n,dimnames = list(rownames(Adj),colnames(Adj)),
                    data =  sapply(1:nrow(dyads),
                                   function(ij) {
                                     i<- dyads[["row"]];j<- dyads[["col"]];
                                     obs.prob_fun(i,j)
                                   }
                    )
  )
  diag(obs.prob)<- 0;P<- obs.prob[Adj.subfun(obs.prob)]
  if(any(P<=0)|any(P>=1)){
    obs.prob[Adj.subfun(obs.prob)]<- proportional.prob(P)
  }
  obs.prob
}

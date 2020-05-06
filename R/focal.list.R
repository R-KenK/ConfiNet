#' Focal scan from theoretical
#'
#' @param scan.theoretical binary matrix (the theoretical scan within do.scan())
#' @param focal either the row index or the name of the focal
#'
#' @return a focal scan as a binary matrix with matching values for the row of the focal, and NAs otherwise, and the name of the focal as a "focal" attribute.
#' @export
#'
#' @examples
#' # Internal use in do.(non.zero.)scan()
focal.scan<- function(scan.theoretical,focal){
  focal.scan<- scan.theoretical;nodes<- rownames(scan.theoretical)
  if(is.character(focal)){
    focal.index<- match(focal,nodes)
    if(is.na(focal.index)){stop("focal name not recognized.")}
    focal.name<- focal
  }else if(is.numeric(focal)){
    focal.index<- focal
    focal.name<- nodes[focal.index]
  }else{stop("focal format unrecognized")}
  focal.scan[-focal.index,]<- NA
  attr(focal.scan,"focal")<- focal.name;
  focal.scan
}

#' Produce focal.list
#'
#' @param Adj square integers matrix of occurences of dyads.
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param focal.prob_fun a user-defined function of (n) that output a probability of being focal for each node
#'
#' @return a vector of focals (as integers)
#' @export
#'
#' @examples
#' set.seed(42)
#' n<- 6;nodes<- as.character(1:n);
#' total_scan<- 20;n.boot<- 5;
#'
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- round(runif(n*(n-1),0,total_scan*.50))
#'
#' make_focal.list(Adj,total_scan)
#' make_focal.list(Adj,total_scan,focal.prob_fun = function(n) 1:n)
#' make_focal.list(Adj,total_scan,focal.prob_fun = function(n) compute.strength(Adj,"directed"))
make_focal.list<- function(Adj,total_scan,
                           focal.prob_fun = NULL){
  n<- nrow(Adj);

  if(is.null(focal.prob_fun)){
    return(ceiling(runif(total_scan,0,n)))
  }
  sample(1:n,total_scan,replace = TRUE,prob = focal.prob_fun(n))
}

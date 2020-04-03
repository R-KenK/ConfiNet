#' Single group scan
#' Generate a random whole group or focal scan based on provided reference adjacency matrix and sampling effort.
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.
#' @param focal Character scalar, indicate which focal to consider for the scan.
#' @param obs.prob either :
#' \itemize{
#'  \item{"a dyad observation obs.probability matrix"}{of same dimension as Adj}
#'  \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
#'  \item{"a general dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
#' }
#' @param keep logical. Relevant if group scans are performed. Indicate if the original "theoretical" group scan should be kept track of.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param method Character scalar, specifies if the function should return a whole-group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument
#' @param prob output of Binary.prob function. Internal use when wrapped with Boot_scans, to speed up bootstrap process
#'
#' @return a square binary matrix representing the whle group scan
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:30,n*(n-1),replace = TRUE)
#' Adj
#'
#' do.scan(Adj,42)
#' do.scan(Adj,42,"c",method = "focal")
#' Reduce("+",lapply(1:42,function(s) do.scan(Adj,42)))
#'
#' focal.list<- sample(nodes,42,replace = TRUE)
#' table(focal.list)
#' L<- lapply(1:42,function(s) do.scan(Adj,42,focal.list[s],method = "both"))
#'
#' list(group = round(Reduce("+",lapply(L,function(l) l$group))/42,3),
#'      focal = {foc<- do.call(rbind,
#'                             lapply(nodes,
#'                                    function(node) {
#'                                      L.focal<- L[sapply(L,function(l) attr(l$focal,"focal"))==node]
#'                                      round(
#'                                            Reduce("+",
#'                                                   lapply(L.focal,
#'                                                          function(l) {
#'                                                            l$focal
#'                                                          }
#'                                                   )
#'                                            )/length(L.focal),3
#'                                      )
#'                                    }
#'                              )
#'                      );
#'               row.names(foc)<-nodes
#'                              foc}
#' )
do.scan<- function(Adj,total_scan,focal=NULL,obs.prob=NULL,keep=FALSE,
                   mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),
                   method = c("group","focal","both"),Adj.subfun=NULL,prob=NULL){
  if(nrow(Adj)==ncol(Adj)) {n<- nrow(Adj);nodes_names<- row.names(Adj)} else {stop("Adj is not a square matrix")}
  mode<- match.arg(mode);
  method<- match.arg(method);

  if(method=="group" & !is.null(focal)){warning("Focal input but group scan chosen as expected method: is it desired behaviour?")}

  if(method %in% c("focal","both") & is.null(focal)){
    warning("No focal input: random one selected.")
    focal<- sample(row.names(Adj),1)
  }

  Scan<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes_names,nodes_names))
  Adj.subfun<- switch(mode,
                      "directed" = ,
                      "undirected" = ,
                      "max" = ,
                      "min" = ,
                      "plus" = non.diagonal,
                      "upper" = upper.tri,
                      "lower" =  lower.tri
  )
  prob<- Binary.prob(Adj=Adj,total_scan=total_scan,mode = mode)
  Scan[Adj.subfun(Scan)]<-  rbinom(nrow(prob),1,prob$present)
  Scan<- binary_adjacency_mode(Scan,mode)

  if(method == "group") {
    if(is.null(obs.prob)){
      return(Scan)
    }else{
      return(observable_edges(Scan = Scan,obs.prob = obs.prob,keep = keep))
    }
  }


  Focal.scan<- Scan[c(focal),];
  attr(Focal.scan,"focal")<- focal;
  if(method == "focal") {return(Focal.scan)}

  if(method == "both") {
    return(
      list(
        "group" = {
          if(is.null(obs.prob)){
            Scan
          }else{
            observable_edges(Scan = Scan,obs.prob = obs.prob,keep = keep)
          }
        },
        "focal" = Focal.scan
      )
    )
  }

  if(!(method %in% c("group","focal","both"))) {stop("How did you reach here?")}
}

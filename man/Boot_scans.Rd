% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Boot_scans.R
\name{Boot_scans}
\alias{Boot_scans}
\title{Bootstrap scans
Bootstrap iterated binary group or focal scans with probabilities derived from a provided adjancecy matrix.}
\usage{
Boot_scans(
  Adj,
  n.boot,
  total_scan,
  method = c("group", "focal", "both"),
  focal.list = NULL,
  scaled = FALSE,
  obs.prob = 1,
  keep = FALSE,
  mode = c("directed", "undirected", "max", "min", "upper", "lower", "plus"),
  output = c("list", "adjacency", "all"),
  n.cores = (parallel::detectCores() - 1),
  cl = NULL
)
}
\arguments{
\item{Adj}{square integers matrix of occurences of dyads. WIP: implement method for association matrices...}

\item{n.boot}{integer, number of bootstrap to perform.}

\item{total_scan}{integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.}

\item{method}{Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).}

\item{focal.list}{Character vector, indicate the list of focals to consider throughout the scans.}

\item{scaled}{logical, specifies if adjacency data should be scaled by sampling effort.}

\item{obs.prob}{either :
\itemize{
 \item{"a dyad observation obs.probability matrix"}{of same dimension as Adj}
 \item{"a dyad observation vector"}{subsetted similarly as Adj (through the non.diagonal() function for instance)}
 \item{"a general dyad observation obs.probability"}{should be in [0,1], assumed to be the case when only one value is inputed)}
}}

\item{keep}{logical. Relevant if group scans are performed. Indicate if the original "theoretical" group scan should be kept track of.}

\item{mode}{Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.}

\item{output}{Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix}

\item{n.cores}{number of threads to use while performingh the bootstrap}

\item{cl}{Optional cluster object (cf snow package), otherwise created according to n.cores.}
}
\value{
according to output and method: a list of iterated scans, or of adjacency matrix
}
\description{
Bootstrap scans
Bootstrap iterated binary group or focal scans with probabilities derived from a provided adjancecy matrix.
}
\examples{
set.seed(42)

n<- 5;nodes<- letters[1:n];
Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
Adj

focal.list<- sample(nodes,42,replace = TRUE)
table(focal.list)

Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,scaled = TRUE,
           method = "group",mode = "directed",output = "list")
Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,scaled = TRUE,
           method = "focal",mode = "directed",output = "adj")
# system.time( #single threaded
#   Boot_scans(Adj,10,total_scan = 42,focal.list = focal.list,scaled = TRUE,
#              method = "both",mode = "directed",output = "adj",n.cores = 1)
# )
# system.time( #multi threaded (with number of thread existing - 1)
#   Boot_scans(Adj,10,total_scan = 42,focal.list = focal.list,
#             scaled = TRUE,obs.prob=0.7,keep=TRUE,
#             method = "both",mode = "directed",output = "adj")
# )
}

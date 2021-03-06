% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix.tools.R
\name{obs.prob_bias}
\alias{obs.prob_bias}
\title{Create probability of edge observation matrix
Create probability of edge observation matrix with user-defined functions. Designed for trait-based or network-based biases.}
\usage{
obs.prob_bias(Adj, obs.prob_fun, bias_fun = NULL, reverse = FALSE)
}
\arguments{
\item{Adj}{a reference adjacency (square) matrix}

\item{obs.prob_fun}{a function to be applied to calculate each cell of the output matrix. The function is applied to the row and column coordinates i and j (e.g. sum(i,j)), or is composed internally with bias_fun() like this: original_fun(bias_fun(i),bias_fun(j))}

\item{bias_fun}{function to implement a bias into the edge observability probability. Can be "trait-based" if related to a node (and in the result to each of its edges), to an edge, or "network-based" if related to a network metric (e.g. a node's strength)}

\item{reverse}{logical. Apply function reverse_i.n() to the output before returning it.}
}
\value{
a matrix of "probability" of edge obervability. The function can return a matrix of any real number (with zeros at the diagonal), and the scaling to [0,1] is currently handled in observable_edges().
}
\description{
Create probability of edge observation matrix
Create probability of edge observation matrix with user-defined functions. Designed for trait-based or network-based biases.
}
\examples{
set.seed(42)
n<- 6;nodes<- as.character(1:n);
total_scan<- 20;n.boot<- 5;
focal.list<- sample(nodes,total_scan,replace = TRUE)

Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
Adj[non.diagonal(Adj)]<- sample((0:round(total_scan*.50)),n*(n-1),replace = TRUE)

traits<- rnorm(nrow(Adj),0,1)
trait.bias_fun<- function(x) {traits[x]}

obs.prob_bias(Adj,sum,bias_fun = trait.bias_fun)
}

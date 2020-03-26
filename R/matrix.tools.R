# matrix subset tools -----------------------------------------------------

#' Non Diagonal Part of a Matrix
#'
#' Similarly to upper.tri and lower.tri, returns a matrix of logicals to identify the non-diagonal of a square matrix
#'
#' @param M square matrix.
#' @param output choose to either return a matrix of logical values (TRUEs, FALSEs on the diagonal, i.e. comparable behaviour than upper.tri for instance) or a vector of the subsetted values of the matrix.
#'
#' @return either return a matrix of logical values (TRUEs, FALSEs on the diagonal) or a vector of the subsetted values of the matrix.
#' @export
#'
#' @examples
#' M<- matrix(sample(1:10,16,replace = TRUE),4,4)
#' non.diagonal(M)
#'
non.diagonal<- function(M,output=c("matrix.logical","vector.values")) {
  output<- match.arg(output)
  if(dim(M)[1]==dim(M)[2]) logicals<- upper.tri(M,diag = FALSE)|lower.tri(M,diag = FALSE) else stop("Matrix provided is not a square matrix.")
  switch(output,
         "matrix.logical" = logicals,
         "vector.values" = M[logicals])
}

#' Diagonal Part of a Matrix
#'
#' Similarly to upper.tri and lower.tri, returns a matrix of logicals to identify the diagonal of a square matrix
#'
#' @param M matrix or other R object with length(dim(x)) == 2. For back compatibility reasons, when the above is not fulfilled, as.matrix(x) is called first.
#' @param output choose to either return a matrix of logical values (TRUEs on the diagonal, i.e. comparable behaviour than diag() for instance) or a vector of the subsetted values of the matrix.
#'
#' @return square logical matrix with diagonal of TRUEs
#' @export
#'
#' @examples
#' M<- matrix(sample(1:10,16,replace = TRUE),4,4)
#' diagonal<- function(M) {upper.tri(M,diag = TRUE)&!upper.tri(M,diag = FALSE)}
#' diagonal(M)
#'
diagonal<- function(M,output=c("matrix.logical","vector.values")) {
  output<- match.arg(output)
  if(dim(M)[1]==dim(M)[2]) logicals<- upper.tri(M,diag = TRUE)&!upper.tri(M,diag = FALSE) else stop("Matrix provided is not a square matrix.")
  switch(output,
         "matrix.logical" = logicals,
         "vector.values" = M[logicals])
}

#' Identify coordinates of non null non diagonal elements of a matrix
#'
#' @param M a square matrix
#'
#' @return an array of row and col indices of the positive and non diagonal elements
#' @export
#'
#' @examples
#' set.seed(42)
#' M<- matrix(sample(0:10,16,replace = TRUE),4,4)
#'
#' non.zero.non.diag(M)
#'
non.zero.non.diag<- function(M) {which(M>0&!diagonal(M),arr.ind = TRUE,useNames = TRUE)}


# matrix sum tool ---------------------------------------------------------

#' Matrix sum removing NAs
#' Equivalent to the elmement-wise matrix addition, but replacing NAs by zeros. Internal use
#'
#' @param X first matrix
#' @param Y second matrix, similarly dimensioned
#'
#' @return the sum of X and Y, replacing NAs by zeros
#' @export
#'
#' @examples
#' #Internal use
matrix_sum_na.rm<-function(X,Y) {ifelse(!is.na(X),X,0)+ifelse(!is.na(Y),Y,0)}

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
n.observed_edges<- function(scan_list,diag=0){
  Reduce("+",
         lapply(scan_list,
                function(scan){
                  observed<- ifelse(!is.na(scan),1,0) # counting part of the algorithm
                  if(!is.null(diag)) {diag(observed)<- diag} # doesn't count the diagonal by default. Left the option to count if self loops should be considered
                  observed
                }
         )
  )
}

# Adjacency mode tools ----------------------------------------------------

#' Make Adjacency fit the selected mode
#' From a directed adjacency matrix, make it fit the selected mode.
#'
#' @param Adj an adjacency matrix
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return an adjacency matrix fitting the selected mode
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:30,n*(n-1),replace = TRUE)
#' Adj<- iterate_scans(Adj,42,method="group",mode="directed",output = "adjacency",n.cores = 1)
#' adjacency_mode(Adj,"max")
adjacency_mode<- function(Adj,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus")){
  mode<- match.arg(mode)
  switch(mode,
         "undirected" = ,
         "max" = ifelse(Adj>=t(Adj),Adj,t(Adj)),
         "min" = ifelse(Adj<=t(Adj),Adj,t(Adj)),
         "plus" = Adj+t(Adj),
         "directed" = ,
         "upper" = ,
         "lower" =  Adj
  )
}

#' Make Binary Adjacency fit the selected mode
#' From a directed binary adjacency matrix, make it fit the selected mode.
#'
#' @param Adj a binary adjacency matrix
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return a binary adjacency matrix fitting the selected mode
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
#' Adj<- do.scan(Adj,42)
#' binary_adjacency_mode(Adj,"plus")
binary_adjacency_mode<- function(Adj,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus")){
  mode<- match.arg(mode)
  switch(mode,
         "undirected" = ,
         "max" = ifelse(Adj+t(Adj)>=1,1,0), #conserve a connection between nodes if there's one in either directions (either adjacency triangle)
         "min" = ifelse(Adj+t(Adj)==2,1,0), #only conserve a connection between nodes who have one in both directions (each adjacency triangle)
         "plus" = Adj+t(Adj),
         "directed" = ,
         "upper" = ,
         "lower" =  Adj
  )
}


# obs.prob tools ----------------------------------------------------------

#' "Reverse" matrix rows and column order
#' For a matrix M with dim(M)=(n,m), swap element M[i,j] swapped by M[n-i+1,m-j+1]
#'
#' @param M a matrix
#'
#' @return a similarly dimensioned matrix with element M[i,j] swapped by M[n-i+1,m-j+1]
#' @export
#'
#' @examples
#' reverse_i.n(matrix(runif(15),3,5))
reverse_i.n<- function(M){
  Reversed<- matrix(0,nrow = nrow(M),ncol = ncol(M),dimnames = list(rownames(M),colnames(M)))
  for(i in 1:nrow(M)){
    for(j in 1:ncol(M)){
      Reversed[i,j]<- M[nrow(M)-i+1,ncol(M)-j+1]
    }
  }
  Reversed
}

#' Create probability of edge observation matrix
#' Create probability of edge observation matrix with user-defined functions. Designed for trait-based or network-based biases.
#'
#' @param Adj a reference adjacency (square) matrix
#' @param obs.prob_fun a function to be applied to calculate each cell of the output matrix. The function is applied to the row and column coordinates i and j (e.g. sum(i,j)), or is composed internally with bias_fun() like this: original_fun(bias_fun(i),bias_fun(j))
#' @param bias_fun function to implement a bias into the edge observability probability. Can be "trait-based" if related to a node (and in the result to each of its edges), to an edge, or "network-based" if related to a network metric (e.g. a node's strength)
#' @param reverse logical. Apply function reverse_i.n() to the output before returning it.
#'
#' @return a matrix of "probability" of edge obervability. The function can return a matrix of any real number (with zeros at the diagonal), and the scaling to [0,1] is currently handled in observable_edges().
#' @export
#'
#' @examples
#' set.seed(42)
#' n<- 6;nodes<- as.character(1:n);
#' total_scan<- 20;n.boot<- 5;
#' focal.list<- sample(nodes,total_scan,replace = TRUE)
#'
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample((0:round(total_scan*.50)),n*(n-1),replace = TRUE)
#'
#' traits<- rnorm(nrow(Adj),0,1)
#' trait.bias_fun<- function(x) {traits[x]}
#'
#' obs.prob_bias(Adj,sum,bias_fun = trait.bias_fun)
obs.prob_bias<- function(Adj,obs.prob_fun,bias_fun = NULL,reverse = FALSE){
  n<- nrow(Adj);
  obs.prob<- matrix(0,n,n,dimnames = list(rownames(Adj),colnames(Adj)));
  if(!is.null(bias_fun)) {
    original_fun<- obs.prob_fun;
    obs.prob_fun<- function(i,j) {original_fun(bias_fun(i),bias_fun(j))};
  }

  for(i in seq_len(n)) {
    for(j in seq_len(n)) {
      if(i!=j){
        obs.prob[i,j]<- obs.prob_fun(i,j);
      }
    }
  }
  if(!reverse) obs.prob else reverse_i.n(obs.prob)
}


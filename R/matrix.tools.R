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

#' Get number of edge observations
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
                  observed<- ifelse(!is.na(scan),1,0)
                  if(!is.null(diag)) {diag(observed)<- diag}
                  observed
                }
         )
  )
}

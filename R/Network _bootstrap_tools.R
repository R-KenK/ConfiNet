

#' TO WRITE
#'
#' @param scan.list TO WRITE
#' @param n.boot TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
bootstrap_scan.list <- function(scan.list,n.boot = 100) {
  replicate(n.boot,resample_scan.list(scan.list),simplify = FALSE)
}

#' TO WRITE
#'
#' @param scan.list TO WRITE
#' @param n.boot TO WRITE
#'
#' @return TO WRITE
#' @export
#'
#' @examples
#' # TO WRITE
bootstrap_metric <- function(scan.list,metric.fun,n.boot = 100) {
  replicate(n.boot,resample_scan.list(scan.list),simplify = FALSE)
}

#' Perform one resampling with replacement
#'
#' @param scan.list a list containing matrices (snPackMat compatible)
#'
#' @return a new list of the same size, but which elements have been resampled with replacement
#' @noRd
resample_scan.list <- function(scan.list) {
  l <- length(scan.list)
  resampled.indices <- quick_sample(1:l,l)
  scan.list[resampled.indices]
}



# Draft of optimization for rare events ----------------------------------------------------------
set.seed(42)
n<- 1000;N<- 10000
V<- round(runif(n,0,5))
P<- Binary.prob(V,N)

#' Title
#'
#' @param P
#' @param ...
#' @param n
#'
#' @return
#' @export
#'
#' @examples
P.cond.in.order<- function(P,...,n = length(P)){
  S<- 1-prod(1-P)
  previous.are.zeros<- c(1,cumprod(1-P[1:(n-1)]))
  sapply(1:n,function(i) P[i]/S*previous.are.zeros[i])
}

#' Title
#'
#' @param P
#' @param ...
#' @param n
#'
#' @return
#' @export
#'
#' @examples
do.non.zero.scan<-function(P,...,n = length(P)){
  X<- rep(0,n)
  P.cond<- cumsum(P.cond.in.order(P))
  first<- min(which(runif(1)<P.cond))
  X[first]<- 1
  X[-first]<- rbinom(n-1,1,P[-first])
  X
}

#' Title
#'
#' @param P
#' @param N
#' @param ...
#' @param n
#'
#' @return
#' @export
#'
#' @examples
iterate_rare.scans<- function(P,N,...,n = length(P)){
  scan<- data.table(matrix(0,N,n));
  non.zero<- rbinom(N,1,1-prod(1-P))==1;
  scan[non.zero,]<- data.table(rbind_lapply(seq_len(nrow(scan[non.zero,])),function(i) do.non.zero.scan(P)))
  scan
}

system.time(X<- iterate_rare.scans(P,N))
plot(colSums(X),V)
abline(0,1,col="red")
cor(colSums(X),V)

system.time(
  X_bis<- rbind_lapply(1:N,
                       function(k){
                         X<-rbinom(n,1,P)
                         if(Reduce("&",X==0)){rep(0,n)} else {X}
                       }
  )
)
plot(colSums(X_bis),V)
abline(0,1,col="red")
cor(colSums(X_bis),V)

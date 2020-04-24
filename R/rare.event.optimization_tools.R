# Rare events optimization related functions ------------------------------

#' Placeholder for calculating expected time using the optimization for rare event
#'
#' @param n number of node of the original network
#' @param total_scan sampling effort
#' @param max.obs maximum value of the weighted adjacency matrix of the original network
#'
#' @return a time value to be compared with the one using the standard approach
#' @export
#'
#' @examples
#' # Internal use in decide_use.rare.opti
opti.expected.time<- function(n,total_scan,max.obs){
  # try and model?
}

#' Placeholder for calculating expected time using the standard scan method
#' i.e. performing all scans
#'
#' @param n number of node of the original network
#' @param total_scan sampling effort
#' @param max.obs maximum value of the weighted adjacency matrix of the original network
#'
#' @return a time value to be compared with the one using the optimization for rare event
#' @export
#'
#' @examples
#' # Internal use in decide_use.rare.opti
standard.expected.time<- function(n,total_scan,max.obs){
  # try and model?
}

#' Decide based on expected times if the otpimization for rare event should be used
#'
#' @param n number of node of the original network
#' @param total_scan sampling effort
#' @param max.obs maximum value of the weighted adjacency matrix of the original network
#'
#' @return a logical value meaning that the optimization for rare events should be use when TRUE is returned
#' @export
#'
#' @examples
#' # Internal use
decide_use.rare.opti<- function(n,total_scan,max.obs=NULL){
  if(is.null(max.obs)&is.matrix(n)){max.obs<- max(n)}
  if(is.matrix(n)){n<- nrow(n)}
  # Figure out how to relate n, N and max.obs through opti.expected.time(n,total_scan,max.obs)
  opti.expected.time(n,total_scan,max.obs)<standard.expected.time(n,total_scan,max.obs)
}

# #' Simulate which scan returns an all-zeros matrix
# #' (old: will probably be deprecated because inefficient)
# #'
# #' @param total_scan integer, sampling effort
# #' @param presence.prob presence probability matrix (or vector)
# #' @param method Character scalar, specifies if the function should return a theoretical perfect group scan, an  empirical group scan (a similarly dimensioned matrix as Adj), or a focal scan (a vector representing the given focal's row in the group scan matrix).
# #'
# #' @return a list of zero-matrices (all-zeros scans) and NULL (non-zeros scans to be later performed)
# #' @export
# #'
# #' @examples
# #' # Internal use
# simulate_zeros.non.zeros.old<- function(total_scan,presence.prob,method){
#   nodes<- rownames(presence.prob);
#   zero.mat<- matrix(0,nrow(presence.prob),ncol(presence.prob),dimnames = list(nodes,nodes))
#   zero.list.element<- switch(method,
#                              "theoretical" = list(theoretical = zero.mat),
#                              "group" = list(theoretical = zero.mat,group = zero.mat),
#                              "focal" = list(theoretical = zero.mat,focal = zero.mat),
#                              "both" = list(theoretical = zero.mat,group = zero.mat,focal = zero.mat)
#   )
#   scan_list<- vector(mode="list",length = total_scan)
#   non.zeros<- rbinom(total_scan,1,1-prod(1-presence.prob))==1;
#   scan_list[!non.zeros]<- lapply(scan_list[!non.zeros],function(scan) zero.list.element)
#   scan_list
# }

#' Simulate which scan returns an all-zeros matrix
#'
#' @param total_scan integer, sampling effort
#' @param presence.prob presence probability matrix (or vector)
#'
#' @return a list of NULL representing the non-zero scan to run with an attribute `n.zero` being the number of full-zero scans
#' @export
#' @importFrom stats rbinom
#'
#' @examples
#' # Internal use
simulate_zeros.non.zeros<- function(total_scan,presence.prob){
  zero.non.zero<- table(ifelse(stats::rbinom(total_scan,1,1-prod(1-presence.prob))==1,"n.non.zeros","n.zeros"));
  scan_list<- vector(mode="list",length = zero.non.zero["n.non.zeros"])
  attr(scan_list,"n.zeros")<- zero.non.zero["n.zeros"]
  scan_list
}

#' Determine step-by-step conditional probabilities for non-zeros scans
#' Internal use. Returns the CDF of the probability that the i-th dyad (with probability presence.prob[i]) is the first to yield a 1.
#'
#' @param presence.prob presence probability matrix (or vector)
#'
#' @return a cumulative distribution function of the probability that the i-th dyad (with probability presence.prob[i]) is the first to yield a 1
#' @export
#'
#' @details Workflow is as follows: first simulate.zeros.non.zeros() determines which scans are all-zeros and which are non-zeros. Then for non zeros, at a random order each dyad is drawn in order with conditional probability that: (1) there is at least one 1 in the scan, and (2) all the previous coins were zeros. Once the first one is drawn, the rest are drawn with their regular probabilities. In details, cumulative density probability of each dyad (in a given random order) to be the first one to be a 1 is calculated, and a random draw determine which one is first, set the previous ones to zero, and draw the rest normally. cf. Supplmenentary material X.
#'
#' @examples
#' # Internal use.
adjust.conditional.prob<- function(presence.prob){
  prob.all.zeros<- 1-prod(1-presence.prob)
  previous.are.zeros<- c(1,cumprod(1-presence.prob[1:(length(presence.prob)-1)]))
  sapply(1:length(presence.prob),function(i) presence.prob[i]/prob.all.zeros*previous.are.zeros[i])
}

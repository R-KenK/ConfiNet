# Simulation workflow tools -----------------------------------------------

#' Make a list of obs.prob matrices given list of functions to be composed
#'
#' @param bias.fun_list function of (i,j,Adj) to conjugate two individual values (e.g. function(i,j,Adj) {bias.fun(bias.subtype(i,Adj),bias.subtype(j,Adj))})
#' @param bias.subtype_list function of (i,Adj) to get an individual value (can be trait- or net-based, or unbiased)(e.g. function(i,j,Adj) {bias.fun(bias.subtype(i,Adj),bias.subtype(j,Adj))})
#' @param type character, either "unbiased" (`bias.subtype_list` passed as a [0,1] numeric), "trait" for trait-based bias (function(i,j,Adj) {bias.fun(bias.subtype(i,Adj),bias.subtype(j,Adj))}), or "net" for net-based bias (centrality<- bias.subtype(Adj);bias.fun(centrality[i],centrality[j])).
#'
#' @return a list of elements containing: `obs.prob_fun` as the function to use to build the matrix obs.prob, and attributes to keep track of (sub)type and functions involved in its making
#' @export
#'
#' @examples
#' # Internal use
make_obs.prob.bias_list<- function(bias.fun_list,bias.subtype_list,type = c("unbiased","trait","net")){
  bias.comb<- expand.grid(fun=seq_along(bias.fun_list),subtype=seq_along(bias.subtype_list))
  lapply(1:nrow(bias.comb),
         function(par){
           bias.fun<- bias.fun_list[bias.comb[par,"fun"]]
           bias.subtype<- bias.subtype_list[bias.comb[par,"subtype"]]
           list(
             obs.prob_fun = switch(type,
                                   "unbiased" = {
                                     function(i,j,Adj) {
                                       n<- nrow(Adj)
                                       obs.prob<- matrix(bias.subtype[[1]](i,j,Adj),n,n)
                                       diag(obs.prob)<-0
                                       obs.prob
                                     }
                                   },
                                   "trait" = {
                                     function(i,j,Adj) {bias.fun[[1]](bias.subtype[[1]](i,Adj),bias.subtype[[1]](j,Adj))}
                                   },
                                   "net" = {
                                     function(i,j,Adj) {
                                       centrality<- bias.subtype[[1]](Adj);
                                       bias.fun[[1]](centrality[i],centrality[j])
                                     }
                                   }
             ),
             type = type,
             subtype = names(bias.subtype),
             fun = names(bias.fun),
             centrality.fun = if(type=="net") {function(Adj) bias.subtype[[1]](Adj)} else {function(Adj) compute.EV(Adj,"max")}
           )
         }
  )
}

#' Gather all obs.prob list from a list of lists
#'
#' @param global_list list of lists containing each a bias.fun_list, bias.subtype_list, and type elements, to be processed internally by `make_obs.prob.bias_list`

#'
#' @return a list of elements containing: `obs.prob_fun` as the function to use to build the matrix obs.prob, and attributes to keep track of (sub)type and functions involved in its making
#' @export
#'
#' @examples
#' # Internal use
make_global_obs.prob<- function(global_list){
  unlist(
    lapply(global_list,
           function(bias_list){
             make_obs.prob.bias_list(bias_list[["bias.fun_list"]],bias_list[["bias.subtype_list"]],type = bias_list[["type"]])
           }

    ),recursive = FALSE
  )
}

#' Add attribute to each obs.prob
#'
#' @param obs.prob the probability of observation matrix
#' @param type its type ("unbiased","trait", or "net")
#' @param subtype its subtype function used to make it
#' @param fun its conjugation function used to make it
#'
#' @return an obs.prob matrix with attributes to keep track of (sub)type and functions involved in its making
#' @export
#'
#' @examples
#' # Internal use
add_obs.prob_attr<- function(obs.prob,type,subtype,fun){
  attr(obs.prob,"type")<- type
  attr(obs.prob,"subtype")<- subtype
  attr(obs.prob,"fun")<- fun
  obs.prob
}

#' Create a list of parameter combinations
#' to be used in simulations of network sampling
#'
#' @param ADJ list of networks' adjacency matrices
#' @param obs.bias_list list of lists containing each a bias.fun_list, bias.subtype_list, and type elements, to be processed internally by `make_obs.prob.bias_list`
#' @param focal.bias_list a (flat) list of functions of (n,Adj) to be passed as focal.list argument of Bootscan (can be also "even" as the default special case).
#'
#' @return a list of lists containing: (1) a matrix `obs.prob` created for each matrix of `ADJ` from functions in `obs.bias_list`, and (2) a `focal.list` being either a function or "even". Both elements have differents attributes to keep track of (sub)type of bias.
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.

initialize_parameter_list<- function(ADJ,obs.bias_list,focal.bias_list){
  obs.bias_list<- make_global_obs.prob(obs.bias_list)
  lapply(seq_along(ADJ),
         function(a){
           cat(paste0("\n",a,"/",length(ADJ)))
           Adj<- ADJ[[a]]
           param.comb<- expand.grid(obs.bias = seq_along(obs.bias_list),focal.bias = seq_along(focal.bias_list))
           lapply(1:nrow(param.comb),
                  function(par){
                    obs.bias<- obs.bias_list[[ param.comb[par,"obs.bias"] ]]
                    obs.prob_fun<- obs.bias[["obs.prob_fun"]]
                    type<- obs.bias[["type"]]
                    subtype<- obs.bias[["subtype"]]
                    fun<- obs.bias[["fun"]]
                    obs.prob<- make_obs.prob(Adj,obs.prob_fun = obs.prob_fun)
                    obs.prob<- add_obs.prob_attr(obs.prob = obs.prob,type = type,subtype = subtype,fun = fun)

                    focal.bias<- focal.bias_list[[ param.comb[par,"focal.bias"] ]]
                    focal.type.name<- names(focal.bias_list[ param.comb[par,"focal.bias"] ])
                    attr(focal.bias,"type")<- gsub("[.].*$","",focal.type.name)
                    attr(focal.bias,"subtype")<- gsub("^.*[.]","",focal.type.name)
                    list(obs.prob = obs.prob,focal.list = focal.bias)
                  }
           )
         }
  )
}

# Visual check of the obs.prob generation ---------------------------------
# n<- 20
# Adj<- matrix((1:(n*n))^2,n,n,byrow = TRUE);diag(Adj)<- 0
# Adj<- matrix(round(runif((n*n),0,5)),n,n,byrow = TRUE);diag(Adj)<- 0
#
# lapply(make_global_obs.prob(obs.bias_list),
#        function(biais){
#          obs.prob_fun<- biais[["obs.prob_fun"]]
#          type<- biais[["type"]]
#          subtype<- biais[["subtype"]]
#          fun<- biais[["fun"]]
#          obs.prob<- make_obs.prob(Adj,obs.prob_fun = obs.prob_fun)
#          obs.prob<- add_obs.prob_attr(obs.prob = obs.prob,type = type,subtype = subtype,fun = fun)
#          par(mfrow=c(2,2))
#          plot(1:n,rowSums(obs.prob),main = paste(attributes(obs.prob)[c("type","subtype","fun")],sep = " "),sub = "P(i)=f(i)")
#          plot(non.diagonal(Adj,"vect"),non.diagonal(obs.prob,"vect"),main = paste(attributes(obs.prob)[c("type","subtype","fun")],sep = " "),sub = "P(i,j)=f(i,j)")
#          plot(1:n,biais[["centrality.fun"]](Adj),main = paste(attributes(obs.prob)[c("type","subtype","fun")],sep = " "),sub = "Centrality(i)=f(i)")
#          plot(compute.EV(Adj,"max"),rowSums(obs.prob),main = paste(attributes(obs.prob)[c("type","subtype","fun")],sep = " "),sub = "P(i)=f(centrality(i))")
#        }
# )
# par(mfrow=c(1,1))


# Bootstrap "object" specific tools ---------------------------------------

#' Keep bootstrap parameters as attributes
#' Internal use in Boot_scans(). Add "method", "keep", "mode", "output" attributes to be more easily retrieved by the get function
#'
#' @param Bootstrap Boot_scans() intermediate output
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#' @param keep logical. Indicate if the original "theoretical" scan should be kept track of.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
#'
#' @return list which structure depends on chosen parameters, with parameters stored as attributes.
#' @export
#'
#' @examples
#' #Internal
Bootstrap_add.attributes<- function(Bootstrap,method,keep,mode,output){
  attr(Bootstrap,"method")<- method;
  attr(Bootstrap,"keep")<- keep;
  attr(Bootstrap,"mode")<- mode;
  attr(Bootstrap,"output")<- output;
  Bootstrap
}

#' Retrieve specific data from Boot_scans() output
#' Subset rich Bootstrap output choosing what's needed
#'
#' @param Bootstrap Bootstrap output object
#' @param what character scalar, data type requested
#' @param format character scalar, output format type requested ("both","list" or "adjacency").
#'
#' @return list which structure depends on chosen data type and Bootstrap attributes
#' @export
#'
#' @examples
#'
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:42,n*(n-1),replace = TRUE)
#' Adj
#'
#' focal.list<- sample(nodes,42,replace = TRUE)
#' Bootstrap<- Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,
#'                        scaled = FALSE,obs.prob=0.7,keep=TRUE,
#'                        method = "group",mode = "directed",output = "list")
#' Boot_get.list(Bootstrap,"theoretical")
#' Boot_get.list(Bootstrap,"observed")
#'
#' Bootstrap<- Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,
#'                        scaled = TRUE,obs.prob=0.7,keep=TRUE,
#'                        method = "both",mode = "directed",output = "all")
#' Boot_get.list(Bootstrap,"observed","both")
#' Boot_get.list(Bootstrap,"observed","adjacency")
#' Boot_get.list(Bootstrap,"observed","list")
Boot_get.list<- function(Bootstrap,what=c("group","focal","theoretical","observed"),format = c("both","list","adjacency")){
  what<- match.arg(what)
  format<- match.arg(format)

  method<- attr(Bootstrap,"method");
  keep<- attr(Bootstrap,"keep");
  output<- attr(Bootstrap,"output");

  if((what %in% c("theoretical","observed")) & (keep==FALSE)) {
    warning(paste0("keep was set to FALSE in Bootstrap. Do you mean to retrieve group?"))
    what<- "group";
  }
  if((what %in% c("group","focal")) & (method!="both" & method!=what)) {stop("Requested List doesn't exist in provided Bootstrap.")}
  if(output=="list"&format=="adjacency") {stop("WIP.")}

  switch(output,
         "list" =  switch(method,
                          "group" = if(keep){
                            lapply(Bootstrap,function(B) lapply(B,function(l) l[[what]]))
                          }else{
                            lapply(Bootstrap,function(B) lapply(B,function(l) l))
                          },
                          "focal" = lapply(Bootstrap,function(B) lapply(B,function(l) l)),
                          "both" = {
                            if(keep){
                              switch(what,
                                     "group" = ,
                                     "focal" = lapply(Bootstrap,function(B) lapply(B,function(l) l[[what]])),
                                     "theoretical" = ,
                                     "observed" = lapply(Bootstrap,function(B) lapply(B,function(l) l[["group"]][[what]])),
                              )
                            }else{
                              lapply(Bootstrap,function(B) lapply(B,function(l) l[[what]]))
                            }
                          }
         ),
         "adjacency" = switch(method,
                              "group" = if(keep){
                                lapply(Bootstrap,function(B) B[[what]])
                              }else{
                                lapply(Bootstrap,function(B) B)
                              },
                              "focal" = lapply(Bootstrap,function(B) B),
                              "both" = {
                                if(keep){
                                  switch(what,
                                         "group" = ,
                                         "focal" = lapply(Bootstrap,function(B) B[[what]]),
                                         "theoretical" = ,
                                         "observed" = lapply(Bootstrap,function(B) B[["group"]][[what]]),
                                  )
                                }else{
                                  lapply(Bootstrap,function(B) B[[what]])
                                }
                              }
         ),
         "all" = {
           all<- list(
             list = switch(method,
                           "group" = if(keep){
                             lapply(Bootstrap,function(B) lapply(B$list,function(l) l[[what]]))
                           }else{
                             lapply(Bootstrap,function(B) lapply(B$list,function(l) l))
                           },
                           "focal" = lapply(Bootstrap,function(B) lapply(B$list,function(l) l)),
                           "both" = {
                             if(keep){
                               switch(what,
                                      "group" = ,
                                      "focal" = lapply(Bootstrap,function(B) lapply(B$list,function(l) l[[what]])),
                                      "theoretical" = ,
                                      "observed" = lapply(Bootstrap,function(B) lapply(B$list,function(l) l[["group"]][[what]])),
                               )
                             }else{
                               lapply(Bootstrap,function(B) lapply(B$list,function(l) l[[what]]))
                             }
                           }
             ),
             adjacency = switch(method,
                                "group" = if(keep){
                                  lapply(Bootstrap,function(B) B$adjacency[[what]])
                                }else{
                                  lapply(Bootstrap,function(B) B$adjacency)
                                },
                                "focal" = lapply(Bootstrap,function(B) B$adjacency),
                                "both" = {
                                  if(keep){
                                    switch(what,
                                           "group" = ,
                                           "focal" = lapply(Bootstrap,function(B) B$adjacency[[what]]),
                                           "theoretical" = ,
                                           "observed" = lapply(Bootstrap,function(B) B$adjacency[["group"]][[what]]),
                                    )
                                  }else{
                                    lapply(Bootstrap,function(B) B$adjacency[[what]])
                                  }
                                }
             )
           )
           switch(format,
                  "both" = all,
                  "list" = all$list,
                  "adjacency" = all$adjacency
           )
         }
  )
}

#' Row bind list of data frames
#' wrapper to one-function do.call rbind over a lapply list
#'
#' @param X a list. See details \link[base]{lapply}.
#' @param FUN a function to subset data frames (or data tables). See details \link[base]{lapply}.
#'
#' @return a row bound data frame
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' X<- lapply(1:3,function(i) list(int = 42,df = data.frame(x = runif(10,0,1),y = runif(10,0,1))))
#' rbind_lapply(X,function(x) x$df)
rbind_lapply<- function(X,FUN){
  do.call(rbind,lapply(X = X,FUN = FUN))
}

#' Bootstrap specific progress bar
#' Provide feedbacks on the simulation testing situation (which parameters, which combination of parameters over the whole list). Internal use.
#'
#' @param p integer, index of the current combination of parameters within parameters.list.
#' @param parameters.list data frame of parameters for all the simulations
#'
#' @return display some feedbacks on the simulation testing situation (which parameters, which combination of parameters over the whole list).
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
boot_progress.param<- function(p,parameters.list = parameters.list){
  cat(paste0("obs.prob = ",attr(parameters.list[[p]]$obs.prob,"name")
             ," - focal.list = ",attr(parameters.list[[p]]$focal.list,"name"),
             " - mode = ",attr(parameters.list[[p]]$mode,"name")," (",p,"/",length(parameters.list),")","\n"))
}

#' Retrieve specific simulation parameters of given Bootstrap
#' Internal use. To ease the recollection of a given bootstrap performed through Boot_scans() iterations alongside a parameters.list
#'
#' @param parameters a list containing a combination of obs.prob, focal.list, and mode parameters. Element of parameters.list
#'
#' @return a data frame referencing the simulation parameters
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
Boot_get.param<- function(parameters){
  data.frame(obs.prob = as.factor(attr(parameters$obs.prob,"name")),
             focal.list = as.factor(attr(parameters$focal.list,"name")),
             mode = as.factor(ifelse(is.null(attr(parameters$mode,"name")),parameters$mode,attr(parameters$mode,"name")))
  )
}

#' Calculate desired data from simulation
#' Internal use. To ease the recollection of a given bootstrap performed through Boot_scans() iterations alongside a parameters.list
#'
#' @param Bootstrap.list Bootstrap output object
#' @param method sampling method to retrieve data for.
#'
#' @return a data frame of desired data
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
Boot_calc.data<- function(Bootstrap.list,method = c("group","focal")){
  data.frame(cor = adjacency_cor(Bootstrap.list = Bootstrap.list,method = method),
             degree = centrality_cor(Bootstrap.list = Bootstrap.list,method = method,centrality.fun = compute.deg),
             strength = centrality_cor(Bootstrap.list = Bootstrap.list,method = method,centrality.fun = compute.strength),
             EV = centrality_cor(Bootstrap.list = Bootstrap.list,method = method,centrality.fun = compute.EV),
             gw = adj_distance(Bootstrap.list = Bootstrap.list,method = method,dist.fun = Ken_gw),
             frob = adj_distance(Bootstrap.list = Bootstrap.list,method = method,dist.fun = Frobenius_from_adjacency),
             L_p.inv = adj_distance(Bootstrap.list = Bootstrap.list,method = method,dist.fun = Pseudo.inv_Laplacian.dist)
             # HERE IMPLEMENT OTHER STATISTICAL APPROACHES: i.e. NETWORK DISTANCES, METRICS CORRELATION
  )
}


# analyses tools ----------------------------------------------------------

#' Calculate the correlation coefficient between "flattened" adjacency matrices
#' Provide the correlation coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#'
#' @param Bootstrap.list an output of Boot_scan()
#' @param method character scalar, indicate if the function should output the coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#'
#' @return a vector of correlation coefficients
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
#' #Internal use in Simulation_script.R.
adjacency_cor<- function(Bootstrap.list,method = c("group","focal")){
  method<- match.arg(method)
  if(method=="group" & attr(Bootstrap.list,"keep")==TRUE) {method<- "observed"}
  n.boot = length(Bootstrap.list)
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           stats::cor(
             c(Boot_get.list(Bootstrap.list,"theoretical","adjacency")[[b]]),   # c() flattens the matrix to consider it like a vector
             c(Boot_get.list(Bootstrap.list,method,"adjacency")[[b]])
           )
         }
  )
}

#' Wrapper to calculate coefficient of correlation between centrality indices
#'
#' @param Bootstrap.list an output of Boot_scan()
#' @param method character scalar, indicate if the function should output the coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#' @param centrality.fun a centrality function taking an igraph object (or an adjacency matrix but transforms it to an igraph object first) and returns an ordered vector of vertex centralities
#'
#' @return a vector of coefficient of correlations between theoretical and empirical methods derived network node centralities.
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
#' #Internal use in Simulation_script.R.
centrality_cor<- function(Bootstrap.list,method = c("group","focal"),centrality.fun){
  method<- match.arg(method)
  if(method=="group" & attr(Bootstrap.list,"keep")==TRUE) {method<- "observed"}
  n.boot = length(Bootstrap.list)
  mode<- attr(Bootstrap.list,"mode")
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           stat::cor(
             centrality.fun(Boot_get.list(Bootstrap.list,"theoretical","adjacency")[[b]],mode=mode),
             centrality.fun(Boot_get.list(Bootstrap.list,method,"adjacency")[[b]],mode=mode)
           )
         }
  )
}

compute.EV<- function(graph,mode){
  if(is.matrix(graph)){graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  EV<- igraph::eigen_centrality(graph, weights = igraph::E(graph)$weight,scale = FALSE)$vector
  names(EV)<- igraph::vertex_attr(graph)[[1]] # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  EV
}
compute.deg<- function(graph,mode){
  if(is.matrix(graph)){graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  deg<- igraph::degree(graph)
  names(deg)<- igraph::vertex_attr(graph)[[1]] # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  deg
}
compute.strength<- function(graph,mode){
  if(is.matrix(graph)){graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  igraph::strength(graph)
}
compute.flowbet<- function(graph,mode){
  if(is.matrix(graph)){graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  graph.network<- intergraph::asNetwork(graph)
  FB<- sna::flowbet(graph.network)
  names(FB)<- igraph::vertex_attr(graph)[[1]] # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  FB
}

#' Wrapper to calculate adjacency/network distances
#'
#' @param Bootstrap.list an output of Boot_scan()
#' @param method character scalar, indicate if the function should output the coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#' @param centrality.fun a centrality function taking an igraph object (or an adjacency matrix but transforms it to an igraph object first) and returns an ordered vector of vertex centralities
#'
#' @return a vector of coefficient of correlations between theoretical and empirical methods derived network node centralities.
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
#' #Internal use in Simulation_script.R.
adj_distance<- function(Bootstrap.list,method = c("group","focal"),dist.fun,mode=NULL){
  method<- match.arg(method)
  if(method=="group" & attr(Bootstrap.list,"keep")==TRUE) {method<- "observed"}
  n.boot = length(Bootstrap.list)
  mode<- attr(Bootstrap.list,"mode")
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           dist.fun(Boot_get.list(Bootstrap.list,"theoretical","adjacency")[[b]],Boot_get.list(Bootstrap.list,method,"adjacency")[[b]],mode = mode)
         }
  )
}
# arg
# adj_distance(tost,"group",Frobenius_from_adjacency)
# adj_distance(tost,"group",Pseudo.inv_Laplacian.dist)
# adj_distance(tost,"group",Ken_gw)

#' Prototype Gorov-Wasserstein network/adjacency matrix distance
#' USE WITH CAUTION. THEORETICAL BACKGROUND NOT WELL DEFINED
#'
#' @param X an adjacency matrix
#' @param Y another adjacency matrix to calculate the distance with X
#'
#' @return a GW distance ?
#' @export
#' @importFrom gwDist solve_FLB_Rglpk
#' @importFrom gwDist gwDist
#'
#' @examples
#' #Internal prototype
Ken_gw<- function(X,Y,mode){
  Adj.subfun<- switch(mode,
                      "directed" = non.diagonal,
                      "plus" = ,
                      "undirected" = ,
                      "max" = ,
                      "min" = ,
                      "upper" = upper.tri,
                      "lower" =  lower.tri
  )
  points_X<- X[Adj.subfun(X)];points_Y<- Y[Adj.subfun(Y)]
  d_X<- sapply(seq_along(points_X), function(i) sapply(seq_along(points_X), function(j) sqrt((points_X[i]-points_X[j])^2)));
  d_Y<- sapply(seq_along(points_Y), function(i) sapply(seq_along(points_Y), function(j) sqrt((points_Y[i]-points_Y[j])^2)));
  prob_X<- prob_Y<- rep(1/length(points_X),length(points_X));
  sol <- gwDist::solve_FLB_Rglpk(X = points_X,Y = points_Y,d_X = d_X,d_Y = d_Y,mu_X = prob_X,mu_Y = prob_Y,p = 1)$solution
  gwDist::gwDist(sol, d_X, d_Y, prob_X, prob_Y)$optimum
}

Frobenius_from_adjacency<- function(X,Y=0,...){
  sqrt(sum((X-Y)^2))
}

# Frobenius_from_adjacency(theo)

L_p.inv<- function(Adj){
  G<-igraph::graph.adjacency(Adj,mode = "max",weighted = TRUE)
  L<- igraph::laplacian_matrix(G,sparse = FALSE)
  corpcor::pseudoinverse(L)
}

Pseudo.inv_Laplacian.dist<- function(X,Y,...){
  L_p.inv.X<- L_p.inv(X)
  L_p.inv.Y<- L_p.inv(Y)
  Frobenius_from_adjacency(L_p.inv.X,L_p.inv.Y)^2/(Frobenius_from_adjacency(L_p.inv.X)*Frobenius_from_adjacency(L_p.inv.Y))
}

# Simulation workflow tools -----------------------------------------------

#' Create a list of parameter combinations
#' to be used in simulations of network sampling
#'
#' @param Adj reference adjacency matrix (especially for network-based observation bias)
#' @param total_scan integer, total number of scan to be performed in each simulation
#' @param n.cores number of cores to use to generate (in parallel) the parameter combinations list
#' @param cl optional cluster object to use, otherwise automatically created when n.cores is inputted
#'
#' @return a list of lists of the different parameter combinations.
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
initialize_parameters<- function(Adj,total_scan,n.cores=(parallel::detectCores()-1),cl=NULL){
  OBS.PROB<- list(trait.pos = obs.prob_bias(Adj = Adj,obs.prob_fun = prod,bias_fun = NULL,reverse = FALSE),
                  trait.neg = obs.prob_bias(Adj = Adj,obs.prob_fun = function(i,j) 1/prod(i,j),bias_fun = NULL,reverse = FALSE),
                  network.pos = obs.prob_bias(Adj = Adj,obs.prob_fun = prod,
                                              bias_fun = function(node) igraph::strength(igraph::graph.adjacency(Adj,weighted = TRUE))[node],
                                              reverse = FALSE),
                  network.neg = obs.prob_bias(Adj = Adj,obs.prob_fun = prod,
                                              bias_fun = function(node) 1/igraph::strength(igraph::graph.adjacency(Adj,weighted = TRUE))[node],
                                              reverse = FALSE)
  )
  OBS.PROB<- c({unb<- seq(0.1,0.9,by = 0.2);names(unb)<- paste0("unbiased_",unb);as.list(unb)},OBS.PROB) # c() over two lists makes them flat while allowing for shorter calls
  MODE<- if(identical(Adj,t(Adj))) {as.list(c("max"))} else  {as.list(c(directed = "directed",max = "max",min = "min",plus = "plus"))}
  FOCAL.LIST<- list(random = sample(rownames(Adj),total_scan,replace = TRUE),  # consider checking if can be implemented for each boot (leaving it NULL?)
                    even = rep_len(rownames(Adj),length.out = total_scan),
                    biased = "TO IMPLEMENT")

  parameters.comb<- expand.grid(
    list(mode = 1:length(MODE),
         focal.list = 1:length(FOCAL.LIST[1:2]),
         obs.prob = 1:length(OBS.PROB)
    )
  )

  if(is.null(cl)) {cl<- snow::makeCluster(n.cores);doSNOW::registerDoSNOW(cl);on.exit(snow::stopCluster(cl))} # left to avoid error if the function is used alone, but should probably be used internally from Boot_scans() now.
  p<-NULL; #irrelevant bit of code, only to remove annoying note in R CMD Check...
  parameters.list<- foreach::`%dopar%`(
    foreach::foreach(p=1:nrow(parameters.comb)),
    list(
      obs.prob = {
        obs.prob<- OBS.PROB[[parameters.comb[p,"obs.prob"]]];
        attr(obs.prob,"name")<- names(OBS.PROB)[parameters.comb[p,"obs.prob"]];
        obs.prob
      },
      focal.list = {
        focal.list<- FOCAL.LIST[[parameters.comb[p,"focal.list"]]];
        attr(focal.list,"name")<- names(FOCAL.LIST)[parameters.comb[p,"focal.list"]];
        focal.list
      },
      mode = {
        mode<- MODE[[parameters.comb[p,"mode"]]];
        attr(mode,"name")<- names(MODE)[parameters.comb[p,"mode"]];
        mode
      }
    )
  )
  parameters.list
}

#' Retrieve data from list of network bootstrap results
#' To be used in a loop/lapply. format them in a ready to be analyzed data frame
#'
#' @param B index of the network from which data are collected.
#' @param Bootstrap.list list of bootstrap of a given network, through different parameters
#' @param parameters.list list of parameters used in bootstrap. Should be loopable similarly to Bootstrap.list
#'
#' @return a data frame of collected data per network and per parameter combination
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
Get.data<- function(B,Bootstrap.list,parameters.list){
  rbind_lapply(seq_along(Bootstrap.list),
               function(b){
                 parameters.list.df<- Boot_get.param(parameters.list[[b]])
                 rbind(data.frame(Network = B,boot = b,method = "group",
                                  parameters.list.df,
                                  Boot_calc.data(Bootstrap.list[[b]],method = "group")),
                       data.frame(Network = B,boot = b,method = "focal",
                                  parameters.list.df,
                                  Boot_calc.data(Bootstrap.list[[b]],method = "focal")))
               }
  )
}
#
# test <- Boot_calc.data(tost,method = "group")
# plot(test$gw,test$frob)
# plot(test$frob,test$gw)
# test$frob.inv<- 1/test$frob
# test$gw.inv<- 1/test$gw
# FactoMineR::PCA(test[,c("cor","strength","EV","gw.inv","frob.inv")])

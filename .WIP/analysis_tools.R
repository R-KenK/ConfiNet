# Analyses tools ----------------------------------------------------------
#' Retrieve data from list of network bootstrap results
#' To be used in a loop/lapply. format them in a ready to be analyzed data frame.
#'
#' @param B index of the network from which data are collected.
#' @param Bootstrap.list list of bootstrap of a given network, through different parameters. Now should receive an unscaled bootstraplist that will also be post scaled.
#'
#' @return a data frame of collected data per network and per parameter combination
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
Get.data<- function(B,Bootstrap.list){
  rbind_lapply(seq_along(Bootstrap.list),
               function(par){
                 Bootstrap<- Bootstrap.list[[par]]
                 parameters.list.df<- Boot_get.param(Bootstrap)
                 Bootstrap.post.scaled<- Bootstrap.list_post.scale(Bootstrap)
                 parameters.list.df.post.scaled<- Boot_get.param(Bootstrap.post.scaled)
                 rbind(
                   data.frame(Network = B,boot = seq_along(Bootstrap),method = "group",
                              parameters.list.df,
                              Boot_calc.data(Bootstrap,method = "group")),
                   data.frame(Network = B,boot = seq_along(Bootstrap),method = "group",
                              parameters.list.df.post.scaled,
                              Boot_calc.data(Bootstrap.post.scaled,method = "group")),
                   data.frame(Network = B,boot = seq_along(Bootstrap),method = "focal",
                              parameters.list.df,
                              Boot_calc.data(Bootstrap,method = "focal")),
                   data.frame(Network = B,boot = seq_along(Bootstrap),method = "focal",
                              parameters.list.df.post.scaled,
                              Boot_calc.data(Bootstrap.post.scaled,method = "focal"))
                 )
               }
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
             # degree = centrality_cor(Bootstrap.list = Bootstrap.list,method = method,centrality.fun = compute.deg),
             strength = centrality_cor(Bootstrap.list = Bootstrap.list,method = method,centrality.fun = compute.strength),
             EV = centrality_cor(Bootstrap.list = Bootstrap.list,method = method,centrality.fun = compute.EV),
             CC = net.metric.diff(Bootstrap.list = Bootstrap.list,method = method,network.fun = weighted.clustering.coeff),
             Frob = adj_distance(Bootstrap.list = Bootstrap.list,method = method,dist.fun = Frobenius_from_adjacency),
             Frob.GOF = adj_gof(Bootstrap.list = Bootstrap.list,method = method,dist.fun = Frobenius_from_adjacency),
             SLap = adj_distance(Bootstrap.list = Bootstrap.list,method = method,dist.fun = Laplacian_spectral.dist),
             SLap.GOF = adj_gof(Bootstrap.list = Bootstrap.list,method = method,dist.fun = Laplacian_spectral.dist),
             obs.cor = adjacency_cor(Bootstrap.list = Bootstrap.list,method = method,get.format = "observed_edges")#,
             # HERE IMPLEMENT OTHER STATISTICAL APPROACHES: i.e. NETWORK DISTANCES, METRICS CORRELATION
  )
}

#' Calculate the correlation coefficient between "flattened" adjacency matrices
#' Provide the correlation coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#'
#' @param Bootstrap.list an output of Boot_scan()
#' @param method character scalar, indicate if the function should output the coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#' @param get.format character scalar, compare requested data with the one from theoretical (default is "adjacency", but now can be "observed_edges").
#'
#' @return a vector of correlation coefficients
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
#' #Internal use in Simulation_script.R.
adjacency_cor<- function(Bootstrap.list,method = c("group","focal"),get.format = c("adjacency","observed_edges")){
  method<- match.arg(method)
  get.format<- match.arg(get.format)
  n.boot = length(Bootstrap.list)
  mode<- Bootstrap_get.attr(Bootstrap.list,"mode")
  if(mode=="directed"){Adj.subfun<- non.diagonal}else{Adj.subfun<- upper.tri}
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           theoretical<- Boot_get.list(Bootstrap.list,"theoretical",get.format)[[b]]
           empirical<- Boot_get.list(Bootstrap.list,method,get.format)[[b]]
           stats::cor(
             theoretical[Adj.subfun(theoretical)],
             empirical[Adj.subfun(empirical)]
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
#' @examples
#' #Internal use in Simulation_script.R.
centrality_cor<- function(Bootstrap.list,method = c("group","focal"),centrality.fun){
  method<- match.arg(method)
  n.boot = length(Bootstrap.list)
  mode<- Bootstrap_get.attr(Bootstrap.list,"mode")
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           cor(
             centrality.fun(Boot_get.list(Bootstrap.list,"theoretical","adjacency")[[b]],mode=mode),
             centrality.fun(Boot_get.list(Bootstrap.list,method,"adjacency")[[b]],mode=mode)
           )
         }
  )
}

#' Compute eigenvector centrality values from graph or adjacency matrix
#'
#' @param graph an igraph object (or an adjacency matrix)
#' @param mode optional, only if `graph` is an adjacency matrix. Othewise character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return a vector of eigenvector centrality values for each node
#' @export
#'
#' @examples
#' # Internal use
compute.EV<- function(graph,mode=NULL){
  if(is.matrix(graph)){graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  EV<- igraph::eigen_centrality(graph, weights = igraph::E(graph)$weight,scale = FALSE)$vector
  if(!is.null(names(EV))){names(EV)<- igraph::vertex_attr(graph)[[1]]} # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  EV
}

#' Compute node degree from graph or adjacency matrix
#'
#' @param graph an igraph object (or an adjacency matrix)
#' @param mode optional, only if `graph` is an adjacency matrix. Othewise character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return a vector of degree values for each node
#' @export
#'
#' @examples
#' # Internal use
compute.deg<- function(graph,mode=NULL){
  if(is.matrix(graph)){graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  deg<- igraph::degree(graph)
  if(!is.null(names(deg))) {names(deg)<- igraph::vertex_attr(graph)[[1]]} # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  deg
}

#' Compute node strength from graph or adjacency matrix
#'
#' @param graph an igraph object (or an adjacency matrix)
#' @param mode optional, only if `graph` is an adjacency matrix. Othewise character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return a vector of strength values for each node
#' @export
#'
#' @examples
#' # Internal use
compute.strength<- function(graph,mode=NULL){
  if(is.matrix(graph)){graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  stren<-igraph::strength(graph)
  if(!is.null(names(stren))) {names(stren)<- igraph::vertex_attr(graph)[[1]]} # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  stren
}

#' Compute node flow betweenness from graph or adjacency matrix
#'
#' @param graph an igraph object (or an adjacency matrix)
#' @param mode optional, only if `graph` is an adjacency matrix. Othewise character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return a vector of flow betweenness values for each node
#' @export
#'
#' @examples
#' # Internal use
compute.flowbet<- function(graph,mode=NULL){
  if(is.matrix(graph)){graph<- igraph::graph.adjacency(graph,mode = mode,weighted = TRUE,add.colnames = TRUE)}
  graph.network<- intergraph::asNetwork(graph)
  FB<- sna::flowbet(graph.network)
  names(FB)<- igraph::vertex_attr(graph)[[1]] # dirty: does not actually test if the order of the vertex centrality is the same as the name, but I suspect igraph does that by default...
  FB
}

#' Wrapper to calculate difference between whole-group network metric
#'
#' @param Bootstrap.list an output of Boot_scan()
#' @param method character scalar, indicate if the function should output the coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#' @param network.fun a whole-group network metric function taking an igraph object (or an adjacency matrix but transforms it to an igraph object first) and returns a metric
#'
#' @return the difference between metrics of networks obtained via theoretical and empirical methods.
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
net.metric.diff<- function(Bootstrap.list,method = c("group","focal"),network.fun){
  method<- match.arg(method)
  n.boot = length(Bootstrap.list)
  mode<- Bootstrap_get.attr(Bootstrap.list,"mode")
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           "-"(
             network.fun(Boot_get.list(Bootstrap.list,method,"adjacency")[[b]],mode=mode),
             network.fun(Boot_get.list(Bootstrap.list,"theoretical","adjacency")[[b]],mode=mode)
           )
         }
  )
}

#' Weighted clustering coefficient
#'
#' @param graph a graph
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return the clustering coefficient
#' @export
#'
#' @importFrom DirectedClustering ClustF
#'
#' @examples
#' # Internal use.
weighted.clustering.coeff<- function(graph,mode){
  if(!is.matrix(graph)){graph<- igraph::get.adjacency(graph,type = "both",sparse = FALSE,attr = "weight")}
  if(isSymmetric.matrix(graph)){
    DirectedClustering::ClustF(mat = graph,type = "undirected")$GlobalCC
  }else{
    DirectedClustering::ClustF(mat = graph,type = "directed")$GlobaltotalCC
  }
}

#' Coefficient of variation of association indices
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return Coefficient of variation of association indices
#' @export
#'
#' @examples
#' # Internal use.
association.index_CV<- function(Adj,mode){
  if(mode=="directed"){Adj<- Adj[non.diagonal(Adj)]}else{Adj<- Adj[upper.tri(Adj)]}
  sd(Adj)/mean(Adj)
}

#' Wrapper to calculate adjacency/network distances
#'
#' @param Bootstrap.list an output of Boot_scan()
#' @param method character scalar, indicate if the function should output the coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#' @param dist.fun a matrix/network distance function taking an igraph object (or an adjacency matrix but transforms it to an igraph object first) and returns an ordered vector of vertex centralities
#'
#' @return a vector of coefficient of correlations between theoretical and empirical methods derived network node centralities.
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
#' #Internal use in Simulation_script.R.
adj_distance<- function(Bootstrap.list,method = c("group","focal"),dist.fun){
  method<- match.arg(method)
  n.boot = length(Bootstrap.list)
  mode<- Bootstrap_get.attr(Bootstrap.list,"mode")
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           dist.fun(Boot_get.list(Bootstrap.list,"theoretical","adjacency")[[b]],Boot_get.list(Bootstrap.list,method,"adjacency")[[b]],mode = mode)
         }
  )
}

#' Wrapper to calculate goodness of fit (GOF) from adjacency/network distances
#'
#' @param Bootstrap.list an output of Boot_scan()
#' @param method character scalar, indicate if the function should output the GOF between the theoretical adjacency matrix and either the empirical one from group or focal method
#' @param dist.fun a matrix/network distance function taking an igraph object (or an adjacency matrix but transforms it to an igraph object first) and returns an ordered vector of vertex centralities
#'
#' @return a vector of GOF values between theoretical and empirical methods.
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
adj_gof<- function(Bootstrap.list,method = c("group","focal"),dist.fun){
  method<- match.arg(method)
  n.boot = length(Bootstrap.list)
  scaled<- Bootstrap_get.attr(Bootstrap.list,"scaled")
  mode<- Bootstrap_get.attr(Bootstrap.list,"mode")
  total_scan<- Bootstrap_get.attr(Bootstrap.list,"total_scan")
  Adj.subfun<- switch(mode,"directed" = ,"undirected" = ,"max" = ,"min" = ,"plus" = non.diagonal,
                      "upper" = upper.tri,"lower" =  lower.tri)
  sapply(1:n.boot,
         function(b) {
           Adj.theo<- Boot_get.list(Bootstrap.list,"theoretical","adjacency")[[b]]
           Adj.method<- Boot_get.list(Bootstrap.list,method,"adjacency")[[b]]
           dist.method<- dist.fun(Adj.theo,Adj.method,mode = mode)
           dist.null<- dist.fun(Adj.theo,generate.null.adj(Adj = Adj.theo,total_scan = total_scan,scaled = scaled,Adj.subfun = Adj.subfun,mode = mode),mode = mode)
           1-dist.method/dist.null # cf. Shore and Lubin (2015) p. 19
         }
  )
}

#' Generate a random matrix as a null following similarly designed matrix generation
#' used in adj_gof
#'
#' @param Adj square integers matrix of occurences of dyads. Optional if using presence.prob. Update: now if presence prob is passed as Adj (thus all(Adj<1) is TRUE), it will be rightfully assigned to presence.prob. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision. Optional if using presence.prob.
#' @param Adj.subfun subsetting function of the adjacency matrix. Driven by igraph "mode" argument
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @importFrom stats rbinom
#'
#' @return a weighted adjacency matrix randomly drawn from binomial trials, comparably to how adjacency matrices are simulated with iterate_scans()
#' @export
#'
#' @examples
#' # internal use in adj_gof
generate.null.adj<- function(Adj,total_scan,
                             Adj.subfun = non.diagonal,scaled = TRUE,
                             mode = c("directed", "undirected", "max","min", "upper", "lower", "plus")){
  n<- nrow(Adj);if(is.null(row.names(Adj))) {nodes<- as.character(1:n)} else {nodes<- row.names(Adj)}
  Null<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
  Null[Adj.subfun(Null)]<- stats::rbinom(length(Null[Adj.subfun(Null)]),total_scan,prob = 0.5)
  adjacency_mode(Null,mode)/ifelse(scaled,total_scan,1)
}

#' Simple Frobenius distance calculation
#'
#' @param X an adjacency matrix
#' @param Y another adjacency matrix, otherwise the distance measuredd is ||X||
#' @param ... for compatibility when mode has to be inputted
#'
#' @return the Frobenius distance
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
# Frobenius_from_adjacency_old<- function(X,Y=0,...){
#   sqrt(sum((X-Y)^2))
# }

#' Frobenius distance calculation
#' considers if the adjacency matrix is symetric or not (and avoid redundant information)
#'
#' @param X an adjacency matrix
#' @param Y another adjacency matrix, otherwise the distance measuredd is ||X||
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. Default here is directed. Possible values are: directed, undirected, upper, lower, max, min, plus. Added vector too. See details \link[igraph]{graph_from_adjacency_matrix}.
#'
#' @return the Frobenius distance
#' @export
#'
#' #Internal use in Simulation_script.R.
Frobenius_from_adjacency<- function(X,Y=0,mode){
  if(mode=="directed"){Adj.subfun<- non.diagonal}else{Adj.subfun<- upper.tri}
  X<- X[Adj.subfun(X)];Y<- Y[Adj.subfun(Y)];
  sqrt(sum((X-Y)^2))
}

Laplacian_spectrum<- function(Adj,mode){
  G<- igraph::graph.adjacency(Adj,weighted = TRUE,mode = mode)
  zapsmall(eigen(igraph::laplacian_matrix(G))$value)
}

#' Graph Laplacian spectral distance
#'
#' @param X an adjacency matrix
#' @param Y another adjacency matrix, otherwise the distance measuredd is ||X||
#' @param ... for compatibility when mode has to be inputted
#'
#' @return the Laplacian spectral distance
#' @export
#'
#' @examples
#' #Internal use for now only.
Laplacian_spectral.dist<- function(X,Y=0,...){
  sqrt(sum(Laplacian_spectrum(X,...)-Laplacian_spectrum(Y,...))^2)
}

# WIP ON HOLD ----
#' #' WIP ON HOLD: revamping Boot_calc.data to be flexible in what it calculates
#' #'
#' #' @param Bootstrap.list
#' #' @param fun.list
#' #' @param args.list
#' #' @param cname.list
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' Boot_calc.data2<- function(Bootstrap.list,fun.list,args.list = NULL,cname.list = NULL){
#'   if(is.null(args.list)){args.list<- attr(fun.list,"args.list")}
#'
#'   df<- cbind_lapply(seq_along(fun.list),
#'                     function(i){
#'                       fun<- fun.list[[i]];args<- args.list[i]
#'                       data.frame(fun(args))
#'                     }
#'   )
#'   if(is.null(cname.list)){cnames<- attr(fun.list,"colnames")}
#'   colnames(df)<- cname.list
#'   df
#' }
#'
#' adjacency_cor_old<- function(Bootstrap.list,method = c("group","focal"),get.format = c("adjacency","observed_edges")){
#'   method<- match.arg(method)
#'   get.format<- match.arg(get.format)
#'   n.boot = length(Bootstrap.list)
#'   sapply(1:n.boot,  # needs function to gather and structure in a data frame
#'          function(b) {
#'            stats::cor(
#'              c(Boot_get.list(Bootstrap.list,"theoretical",get.format)[[b]]),   # c() flattens the matrix to consider it like a vector
#'              c(Boot_get.list(Bootstrap.list,method,get.format)[[b]])
#'            )
#'          }
#'   )
#' }
#'
#' #' Title
#' #'
#' #' @param Bootsrap.df
#' #' @param ADJ
#' #' @param net.fun.list
#' #' @param mode
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' Boot_merge.from_Adj<- function(Bootsrap.df,ADJ,net.fun.list,mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),cnames = NULL){
#'   mode<- match.arg(mode)
#'
#'   net.metrics.df<- rbind_lapply(unique(Bootsrap.df[["Network"]]),
#'                function(i){
#'                  Adj<- ADJ[[i]]
#'                  graph<- igraph::graph_from_adjacency_matrix(Adj,mode = mode)
#'                  net.metrics<- do.call(cbind,
#'                                           lapply(net.fun.list,
#'                                                  function(net.fun){
#'                                                    data.frame(Network = i,net.fun(graph))
#'                                                  }
#'                                           )
#'                  )
#'                  if(is.null(cnames)){cnames<- attr(net.fun.list,"colnames")}
#'                  colnames(net.metrics)[-c(1)]<- cnames
#'                  net.metrics
#'                }
#'   )
#'   merge(Bootsrap.df,net.metrics.df,by = "Network")
#' }
#'
#' make_net.fun.list<- function(...){
#'
#' }

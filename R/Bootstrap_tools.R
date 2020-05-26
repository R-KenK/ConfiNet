
# Bootstrap object data management ----------------------------------------

#' Keep bootstrap parameters (and more) as attributes
#' Internal use in Boot_scans(). Add "method", "keep", "mode", "output" attributes to be more easily retrieved by the get function
#'
#' @param Bootstrap Boot_scans() intermediate output
#' @param ... Any named attribute to be included to the Bootstrap object. Mostly determined by the code of Boot_scan(). Will produce error later if no name is attributed to any
#'
#' @return Bootstrap object, with stored attribute for later retrieval through Bootstrap_get.attributes(). A special attribute, attr.list stores those who have been properly inputted, to account for those with special structure.
#' @export
#'
#' @examples
#' #Internal
Bootstrap_add.attributes<- function(Bootstrap,...){
  attr(Bootstrap,"attr.list")<- list(...)
  Bootstrap
}

#' Retrieve bootstrap object's attributes and reassign them in caller frame
#'
#' @param Bootstrap Boot_scans() output
#' @param a character, parameter stored as attribute to retrieve from `Bootstrap`.
#'
#' @return the value of the parameter `a`.
#' @export
#'
#' @examples
#' #Internal
Bootstrap_get.attr<- function(Bootstrap,a){
  attr(Bootstrap,"attr.list")[[a]]
}


#' Retrieve data from specific method from Boot_scans() output
#' Subset rich Bootstrap output choosing what's needed
#'
#' @param Bootstrap Bootstrap output object
#' @param what character scalar, data type requested ("theoretical","group" or "focal")
#' @param get.format character scalar, output format type requested ("list","adjacency", "observed_edges" or "all").
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
#'                        scaled = FALSE,obs.prob=0.7,
#'                        method = "group",mode = "directed",output = "list")
#' Boot_get.list(Bootstrap,"theoretical")
#' Boot_get.list(Bootstrap,"group")
#' Boot_get.list(Bootstrap,"group","adj")
#' Boot_get.list(Bootstrap,"group","obs")
#'
#' Bootstrap<- Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,
#'                        scaled = FALSE,obs.prob=0.7,
#'                        method = "group",mode = "directed",output = "adjacency")
#' Boot_get.list(Bootstrap,"theoretical","adj")
#' # Boot_get.list(Bootstrap,"group","list")
#' Boot_get.list(Bootstrap,"group","adj")
#' Boot_get.list(Bootstrap,"group","obs")
#'
#' Bootstrap<- Boot_scans(Adj,3,total_scan = 42,focal.list = focal.list,
#'                        scaled = TRUE,obs.prob=0.7,keep=TRUE,
#'                        method = "both",mode = "directed",output = "all")
#' Boot_get.list(Bootstrap,"focal","all")
#' Boot_get.list(Bootstrap,"group","adjacency")
#' Boot_get.list(Bootstrap,"focal","obs")
#' Boot_get.list(Bootstrap,"group","list")
Boot_get.list<- function(Bootstrap,what=c("theoretical","group","focal"),
                         get.format = c("list","adjacency","observed_edges","all")){
  what<- match.arg(what)
  get.format<- match.arg(get.format)

  method<- Bootstrap_get.attr(Bootstrap,"method");
  output<- Bootstrap_get.attr(Bootstrap,"output");
  total_scan<- Bootstrap_get.attr(Bootstrap,"total_scan");
  n.boot<- Bootstrap_get.attr(Bootstrap,"n.boot");
  scaled<- Bootstrap_get.attr(Bootstrap,"scaled");
  mode<- Bootstrap_get.attr(Bootstrap,"mode");
  use.rare.opti<- Bootstrap_get.attr(Bootstrap,"use.rare.opti");
  focal.list<- Bootstrap_get.attr(Bootstrap,"focal.list");

  if(what!="theoretical" & method!="both" & what!=method){stop("Element requested unavailable in `",substitute(Bootstrap),"`.")}

  switch(output,
         "list" = switch(get.format,
                         "list" = lapply(Bootstrap,function(boot) lapply(boot, function(scan) scan[[what]])),
                         "adjacency" = lapply(Bootstrap,
                                              function(boot){
                                                sum_up.scans(scan_list = boot,focal.list = focal.list,
                                                             scaled = scaled,method = what,mode = mode,use.rare.opti = use.rare.opti)[[what]]
                                              }
                         ),
                         "observed_edges" = lapply(Bootstrap,
                                                   function(boot){
                                                     attr(
                                                       sum_up.scans(scan_list = boot,focal.list = focal.list,
                                                                  scaled = scaled,method = what,mode = mode,use.rare.opti = use.rare.opti)[[what]],
                                                       "observed_edges"
                                                     )
                                                   }
                         ),
                         "all" = {
                           list(
                             list = lapply(Bootstrap,function(boot) lapply(boot, function(scan) scan[[what]])),
                             adjacency = lapply(Bootstrap,
                                                function(boot){
                                                  sum_up.scans(scan_list = lapply(boot, function(scan) scan[[what]]),
                                                               scaled = scaled,method = what,mode = mode,use.rare.opti = use.rare.opti)
                                                }
                             )
                           )
                         }
         ),
         "adjacency" = switch(get.format,
                              "list" = stop("Only summed-up adjacency matrices have been stored in Bootstrap object."),
                              "adjacency" = lapply(Bootstrap,function(boot) boot[[what]]),
                              "observed_edges" = lapply(Bootstrap,function(boot) attr(boot[[what]],"observed_edges")),
                              "all" = {
                                warning("Only summed-up adjacency matrices have been stored in Bootstrap object.")
                                lapply(Bootstrap,function(boot) boot[[what]])
                              }
         ),
         "all" = switch(get.format,
                        "list" = lapply(Bootstrap,function(boot) lapply(boot[["list"]], function(scan) scan[[what]])),
                        "adjacency" = lapply(Bootstrap,function(boot) boot[["adjacency"]][[what]]),
                        "observed_edges" = lapply(Bootstrap,function(boot) attr(boot[["adjacency"]][[what]],"observed_edges")),
                        "all" = {
                          list(
                            list = lapply(Bootstrap,function(boot) lapply(boot[["list"]], function(scan) scan[[what]])),
                            adjancency = lapply(Bootstrap,function(boot) boot[["adjacency"]][[what]])
                          )
                        }
         )
  )
}

#' Scale unscaled bootstrap list
#'
#' @param Bootstrap Bootstrap output object (/!\ WITH output = "adjacency", other output to be implemented if relevant)
#'
#' @return a Bootstrap output object with scaled adjacencies (comparable )
#' @export
#'
#' @examples
#' # Internal use.
Bootstrap.list_post.scale<- function(Bootstrap){
  attr.list<- attr(Bootstrap,"attr.list");
  attr.list[["output"]]<- "adjacency"
  attr.list[["scaled"]]<- TRUE

  Bootstrap.scaled<- lapply(Bootstrap,
                            function(boot){
                              Adj.scaled<- lapply(names(boot),
                                                  function(method){
                                                    observed_edges<- attr(boot[[method]],"observed_edges")
                                                    observed_edges<- ifelse(observed_edges!=0,observed_edges,1) # supposedly boot[[method]]==0 when observed_edges==0, but avoid dividing by 0
                                                    scaled<- boot[[method]]/observed_edges
                                                    diag(scaled)<-0
                                                    scaled
                                                  }
                              )
                              names(Adj.scaled)<- names(boot)
                              Adj.scaled
                            }
  )
  attr(Bootstrap.scaled,"attr.list")<- attr.list
  Bootstrap.scaled
}


#' Retrieve specific simulation parameters of given Bootstrap
#' Internal use. To ease the recollection of a given bootstrap performed through Boot_scans() iterations alongside a parameters.list
#'
#' @param Bootstrap a Bootstrap object (with correct parameter attributes)
#'
#' @return a data frame referencing the simulation parameters
#' @export
#'
#' @examples
#' #Internal use in Simulation_script.R.
Boot_get.param<- function(Bootstrap){
  method<- Bootstrap_get.attr(Bootstrap,"method");
  output<- Bootstrap_get.attr(Bootstrap,"output");
  total_scan<- Bootstrap_get.attr(Bootstrap,"total_scan");
  n.boot<- Bootstrap_get.attr(Bootstrap,"n.boot");
  scaled<- Bootstrap_get.attr(Bootstrap,"scaled");
  mode<- Bootstrap_get.attr(Bootstrap,"mode");
  use.rare.opti<- Bootstrap_get.attr(Bootstrap,"use.rare.opti");

  obs.prob<- Bootstrap_get.attr(Bootstrap,"obs.prob");
  focal.list<- Bootstrap_get.attr(Bootstrap,"focal.list");

  data.frame(obs.prob.type = attr(obs.prob,"type"),
             obs.prob.subtype = attr(obs.prob,"subtype"),
             obs.prob.fun = attr(obs.prob,"fun"),
             focal.list.type = attr(focal.list,"type"),
             focal.list.subtype = attr(focal.list,"subtype"),
             total_scan = total_scan,
             mode = mode,
             scaled = scaled,
             use.rare.opti = use.rare.opti
  )
}

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

# analyses tools ----------------------------------------------------------

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
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           stats::cor(
             c(Boot_get.list(Bootstrap.list,"theoretical",get.format)[[b]]),   # c() flattens the matrix to consider it like a vector
             c(Boot_get.list(Bootstrap.list,method,get.format)[[b]])
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
#
# generate.null.adj<- function(Adj,total_scan,
#                              Adj.subfun = non.diagonal,scaled = TRUE,
#                              mode = c("directed", "undirected", "max","min", "upper", "lower", "plus")){
#   n<- nrow(Adj);if(is.null(row.names(Adj))) {nodes<- as.character(1:n)} else {nodes<- row.names(Adj)}
#   Null<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#   Null[Adj.subfun(Null)]<- rbinom(length(Null[Adj.subfun(Null)]),total_scan,prob = 0.5)
#   adjacency_mode(Null,mode)/ifelse(scaled,total_scan,1)
# }

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
Frobenius_from_adjacency<- function(X,Y=0,...){
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

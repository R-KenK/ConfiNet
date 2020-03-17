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
#' @param Bootstrap Boot_scans() final output
#' @param what character scalar, data type requested
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
#' Boot_get.list(Bootstrap,"focal")$adjacency
#' Boot_get.list(Bootstrap,"observed")$list
Boot_get.list<- function(Bootstrap,what=c("group","focal","theoretical","observed")){
  what<- match.arg(what)

  method<- attr(Bootstrap,"method");
  keep<- attr(Bootstrap,"keep");
  output<- attr(Bootstrap,"output");

  if((what %in% c("theoretical","observed")) & (keep==FALSE)) {
    warning(paste0("keep was set to FALSE in Bootstrap. Do you mean to retrieve group?"))
    what<- "group";
  }
  if((what %in% c("group","focal")) & (method!="both" & method!=what)) {stop("Requested List doesn't exist in provided Bootstrap.")}

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
         "all" = list(
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
             mode = as.factor(attr(parameters$mode,"name"))
  )
}

#' Calculate the correlation coefficient between "flattened" adjacency matrices
#' Provide the correlation coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#'
#' @param Bootstrap.list an output of Boot_scan()
#' @param n.boot number of Bootstrap to go through
#' @param what character scalar, indicate if the function should output the coefficient between the theoretical adjacency matrix and either the empirical one from group or focal method
#'
#' @return a vector of correlation coefficients
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
#' #Internal use in Simulation_script.R.
adjacency_cor<- function(Bootstrap.list,what = c("observed","focal"),n.boot = length(Bootstrap.list)){
  what<- match.arg(what)
  sapply(1:n.boot,  # needs function to gather and structure in a data frame
         function(b) {
           stats::cor(c(Boot_get.list(Bootstrap.list,"theoretical")$adjacency[[b]]),   # c() flattens the matrix to consider it like a vector
                      c(Boot_get.list(Bootstrap.list,what)$adjacency[[b]]))
         }
  )
}
# HERE IMPLEMENT OTHER STATISTICAL APPROACHES: i.e. NETWORK DISTANCES, METRICS CORRELATION

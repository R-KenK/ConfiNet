#' Bootstrap scans
#' Perform a bootstrap of group scans with probabilities derived from a provided adjancecy matrix, to produce a new adjancecy matrix.
#'
#' @param Adj square integers matrix of occurences of dyads. WIP: implement method for association matrices...
#' @param total_scan integer, sampling effort. Note that 1/total_scan should be relatively small, increasingly small with increasing precision.
#' @param method Character scalar, specify if the function should use a whole group or a focal scan sampling method (or both).
#' @param focal.list Character vector, indicate the list of focals to consider throughout the scans.
#' @param scaled logical, specifies if adjacency data should be scaled by sampling effort.
#' @param mode Character scalar, specifies how igraph should interpret the supplied matrix. See also the weighted argument, the interpretation depends on that too. Possible values are: directed, undirected, upper, lower, max, min, plus. See details \link[igraph]{graph_from_adjacency_matrix}.
#' @param output Character scalar, specify if the function should return the list of scans, or reduce them into the bootstrapped adjacency matrix
#' @param n.cores number of threads to use while performingh the bootstrap
#' @param cl Optional cluster object (cf snow package), experimentally set to put the makeCluster and stopCluster out of the bootable function. (WIP, next implementation should rethink this).
#'
#' @return according to output a list of scans, or the bootstrapped adjacency matrix
#' @export
#' @importFrom parallel detectCores
#' @importFrom snow makeCluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach `%dopar%`
#' @importFrom foreach foreach
#'
#' @examples
#' set.seed(42)
#'
#' n<- 5;nodes<- letters[1:n];
#' Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
#' Adj[non.diagonal(Adj)]<- sample(0:80,n*(n-1),replace = TRUE)
#' Adj
#'
#' focal.list<- sample(nodes,80,replace = TRUE)
#' table(focal.list)
#' Boot_scan(Adj,80,scaled = FALSE,method = "group",output = "adj",n.cores = 1)
#' Boot_scan(Adj,80,focal.list,scaled = TRUE,
#'                 method = "focal",mode = "directed",output = "list",n.cores = 1)
#'
#' Bootstrap<- lapply(
#'   1:5,
#'     function(b){
#'         Boot_scan(Adj = Adj,total_scan = 80,focal.list = focal.list,scaled = TRUE,
#'                         method = "both",mode = "directed",output = "adj",n.cores = 1)
#'     }
#' )

Boot_scan<- function(Adj,total_scan,method=c("group","focal","both"),focal.list=NULL,
                     scaled=FALSE,
                     mode = c("directed", "undirected", "max","min", "upper", "lower", "plus"),
                     output=c("list","adjacency"),n.cores=(parallel::detectCores()-1),cl=NULL){
  b<-NULL; #irrelevant bit of code, only to remove annoying note in R CMD Check...
  method<- match.arg(method)
  output<- match.arg(output)
  mode<- match.arg(mode)

  if(is.null(focal.list)){
    focal.list<- sample(rownames(Adj),total_scan,replace=TRUE)
  }else{
    if(length(focal.list)!=total_scan) {stop("Provided focal.list's length doesn't match total number of scans to perform.")}
  }

  if(is.null(cl)){
    cl<- snow::makeCluster(n.cores)
    doSNOW::registerDoSNOW(cl);on.exit(snow::stopCluster(cl))
  }



  switch(method,
         "group" = {
           boot.list<- foreach::`%dopar%`(
             foreach::foreach(b=1:total_scan,
                              .export = c("do.scan","non.diagonal","Binary.prob")),
             do.scan(Adj = Adj,total_scan = total_scan,
                     focal = NULL,
                     mode = mode, output = "group")
           )
         },
         "focal" = {
           boot.list<- foreach::`%dopar%`(
             foreach::foreach(b=1:total_scan,
                              .export = c("do.scan","non.diagonal","Binary.prob")),
             do.scan(Adj = Adj,total_scan = total_scan,
                     focal = focal.list[b],
                     mode = mode,output = "focal")
           )
         },
         "both" = {
           boot.list<- foreach::`%dopar%`(
             foreach::foreach(b=1:total_scan,
                              .export = c("do.scan","non.diagonal","Binary.prob")),
             do.scan(Adj = Adj,total_scan = total_scan,
                     focal = focal.list[b],
                     mode = mode,output = "both")
           )
         }
  )

  switch(output,
         "list" = return(boot.list),
         "adjacency" = return(
           switch(method,
                  "group" = {Reduce("+",boot.list)/ifelse(scaled,length(boot.list),1)},
                  "focal" = {
                    boot.adj<- do.call(rbind,
                                       lapply(rownames(Adj),
                                              function(node) {
                                                boot.list.focal<- boot.list[sapply(boot.list,function(b) attr(b,"focal"))==node]
                                                Reduce("+",boot.list.focal)/ifelse(scaled,length(boot.list.focal),1)
                                              }

                                       )
                    );
                    row.names(boot.adj)<-rownames(Adj)
                    boot.adj
                  },
                  "both" = {
                    list(
                      group = {
                        boot.list.group<- lapply(boot.list,function(l) l$group);
                        Reduce("+",boot.list.group)/ifelse(scaled,length(boot.list.group),1)
                      },
                      focal = {
                        boot.adj<- do.call(rbind,
                                           lapply(rownames(Adj),
                                                  function(node) {
                                                    boot.list.focal<- boot.list[sapply(boot.list,function(l) attr(l$focal,"focal"))==node]
                                                    boot.list.focal<- lapply(boot.list.focal,function(l) l$focal)
                                                    Reduce("+",boot.list.focal)/ifelse(scaled,length(boot.list.focal),1)
                                                  }
                                           )
                        );
                        row.names(boot.adj)<- rownames(Adj)
                        boot.adj
                      }
                    )
                  }
           )
         )
  )
}

# n<- 15;nodes<- letters[1:n];
# Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
# Adj[non.diagonal(Adj)]<- sample(0:80,n*(n-1),replace = TRUE)
# Adj
#
# focal.list<- sample(nodes,80,replace = TRUE)
# table(focal.list)
#
# Boot_scan(Adj = Adj,total_scan = 80,focal.list = focal.list,scaled = TRUE,
#                 method = "both",mode = "directed",output = "adj",n.cores = 4)
#
# cl<- snow::makeCluster(7)
# doSNOW::registerDoSNOW(cl);
# system.time(
#   Boots<- lapply(
#     1:10,
#     function(b){
#       Boot_scan(Adj = Adj,total_scan = 80,focal.list = focal.list,scaled = TRUE,
#                       method = "both",mode = "directed",output = "adj",n.cores = 7,cl=cl)
#     }
#   )
# )
# snow::stopCluster(cl)
#
# list(group = round(Reduce("+",lapply(boot.list,function(l) l$group))/42,3),
#      focal = {test<- do.call(rbind,
#                              lapply(nodes,
#                                     function(node) {
#                                       boot.list.focal<- boot.list[sapply(boot.list,function(l) attr(l$focal,"focal"))==node]
#                                       round(
#                                         Reduce("+",
#                                                lapply(boot.list.focal,
#                                                       function(l) {
#                                                         l$focal
#                                                       }
#                                                )
#                                         )/length(boot.list.focal),3
#                                       )
#                                     }
#
#                              )
#      );
#      row.names(test)<-nodes
#      test}
#
# )

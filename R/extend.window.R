# Extend the time-window
#
# Extend the provided time-window to ensure that the extended one contains at least the minimum number of observation.
#
# @param Obs data.table of dyadic observations.
# @param di inferior date boundary of the desired time-window. Format should be "YYYY-MM-DD"
# @param ds superior date boundary of the desired time-window. Format should be "YYYY-MM-DD"
# @param min.obs minimum number (integer) of observations (rows) desired in the returned time-window, taking into consideration that the output might be randomly subsampled keeping only pctdraw \% of Obs.
# @param rdraw logical. Indicates if the output will be randomly subsampled or not.
# @param pctdraw if rdraw is TRUE, the proportion in (0,1) of Obs that will be kept in the later subsample
# @param mode character string, either "dates" or "data.table", to choose the desirted output format
#
# @return Returns either the time-window as a subset data table of dyadic observation, or a vector of the two (closest) dates so that a time-window bounded by these will contain at least min.obs rows.
# @export
# @importFrom data.table data.table
#
# @examples
# set.seed(42);ID<- letters[1:10];
#
# D.Obs<- round(runif(n = 40,min = 0,max = 15))
# Obs<- data.table::data.table(date=as.Date(D.Obs-15,origin = "1991-03-30"),
#                              id=sample(ID,40,replace = TRUE),
#                              tar=sample(ID,40,replace = TRUE))[order(date)]
#
# extend.window(Obs,"1991-03-20","1991-03-25",15)
# extend.window(Obs,"1991-03-20","1991-03-25",20,rdraw = TRUE,pctdraw = .90,mode="data.table")
# extend.window<- function(Obs,di,ds,min.obs,rdraw=FALSE,pctdraw=1,mode="dates"){
#   min.obs<- min.obs/ifelse(rdraw,pctdraw,1)
#
#
#   if(!exists("OD")) OD<- unique(Obs[order(date)]$date);if(is.null(Obs$date)) stop("No date column in Obs.")
#   if(nrow(Obs)<min.obs) stop("Minimum set too low, even when window is extended to the max.")
#   if(!(di %in% OD)|!(ds %in% OD)) {di<- closest.date(OD,di);ds<- closest.date(OD,ds);}
#
#   while (nrow(Obs[date>=di&date<=ds])<min.obs) {
#     if(sample(0:1,1)==0){
#       ifelse(rk.closest.date(OD,di)>1,di<- OD[rk.closest.date(OD,di)-1],ds<- OD[rk.closest.date(OD,ds)+1])
#     }else{
#       ifelse(rk.closest.date(OD,ds)<length(OD),ds<- OD[rk.closest.date(OD,ds)+1],di<- OD[rk.closest.date(OD,di)-1])
#     }
#   }
#   switch (mode,
#           "d" =,
#           "D" =,
#           "date" =,
#           "Date" =,
#           "Dates" =,
#           "dates" = return(c(di, ds)),
#           "dt" =,
#           "DT" =,
#           "data.frame" =,
#           "data frame" =,
#           "data table" =,
#           "data.table" = return(Obs[date>=di&date<=ds]),
#           stop("Incorrect mode. Accepted ones are dates and data.table")
#   )
# }

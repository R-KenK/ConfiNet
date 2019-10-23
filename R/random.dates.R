#' Random dates
#'
#' Randomly draw dates to bound random time-windows.
#'
#' @param Obs data.table of dyadic observations. Must contain a date column.
#' @param duration desired duration of the time-window, in unit.
#' @param extension bidirectional extension of the time-window.
#' @param unit unit in duration and extension are provided, with the smallest unit considered being days to perform Dates calculations. (EXTEND TO OTHER TIME RESOLUTION LATER?)
#' @param mode type of random draw. So far only "first-half" is implemented: draws an inferior date in the first half of the existing dates, then draw the closest date that is duration+extension (in days) away.
#'
#' @return Returns a vector of two random (closest) dates from Obs, separated apart from ~duration+extension days.
#' @export
#' @importFrom data.table data.table
#'
#' @examples
#' set.seed(42);ID<- letters[1:10];
#' D.Obs<- round(runif(n = 100,min = 0,max = 60))
#' Obs<- data.table::data.table(date=as.Date(D.Obs,origin = "1991-03-30"),
#'                              id=sample(ID,100,replace = TRUE),
#'                              tar=sample(ID,100,replace = TRUE))[order(date)]
#' random.dates(Obs,10,5)

random.dates<- function(Obs, duration=4,extension=2,unit="day",mode="first-half"){
  data.table::data.table(Obs)
  switch (unit,
          "day" = unit<- 1,
          "week" = unit<- 7,
          "month" = unit<- 30.5,
          stop("Incorrect unit. Accepted units are: day, week, and month")
  )
  if(is.null(Obs$date)) stop("No date column in Obs.")

  OD<- unique(Obs[order(date)]$date)
  switch (mode,    # TODO: implement, if makes sense, other ways of randomly draw dates for RTWS
    "first-half" = {di<-OD[sample(1:round(length(OD)/2,0),1)];ds<- OD[rk.closest.date(OD,di+floor((duration+extension)*unit))];}
  )
  c(di,ds)
}

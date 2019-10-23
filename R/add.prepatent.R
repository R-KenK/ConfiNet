#' Add prepatent period to data table
#'
#' @param DT data table containing at least a date column.
#' @param mean.lag time lag substracted from date, in unit.
#' @param duration desired duration of the time-window, in unit.
#' @param extension bidirectional extension of the time-window, especially for a subsequent random time-window shift.
#' @param unit unit in which mean.lag and duration are provided, with the smallest unit considered being days to perform Dates calculations. (EXTEND TO OTHER TIME RESOLUTION LATER?)
#'
#' @return returns the original data table TD after appending to it columns dinf and dsup, representing the boundaries of the small or large time-window (according to extension being zero or not) for each datapoint
#'
#' @export
#' @importFrom data.table data.table
#' @examples
#' set.seed(42)
#' D<- round(runif(n = 10,min = 1,max = 30))
#' DT<- data.table::data.table(date=as.Date(D,origin = "1991-03-30"))
#' add.prepatent(DT,mean.lag = 10, duration = 12,extension = 5,unit = "day")
add.prepatent<- function(DT,mean.lag,duration=4,extension=2,unit=c("day","week","month")){
  switch (unit,
          "day" = unit<- 1,
          "week" = unit<- 7,
          "month" = unit<- 30.5,
          stop("Incorrect unit. Accepted units are: day, week, and month")
  )
  if(is.null(DT$date)) stop("No date column in DT.")

  DT$dinf<- DT$date-floor((mean.lag+duration/2+extension)*unit)
  DT$dsup<- DT$date-ceiling((mean.lag-duration/2-extension)*unit)
  DT
}

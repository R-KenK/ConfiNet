% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.prepatent.R
\name{add.prepatent}
\alias{add.prepatent}
\title{Add prepatent period to data table}
\usage{
add.prepatent(
  DT,
  mean.lag,
  duration = 4,
  extension = 2,
  unit = c("day", "week", "month")
)
}
\arguments{
\item{DT}{data table containing at least a date column.}

\item{mean.lag}{time lag substracted from date, in unit.}

\item{duration}{desired duration of the time-window, in unit.}

\item{extension}{bidirectional extension of the time-window, especially for a subsequent random time-window shift.}

\item{unit}{unit in which mean.lag and duration are provided, with the smallest unit considered being days to perform Dates calculations. (EXTEND TO OTHER TIME RESOLUTION LATER?)}
}
\value{
returns the original data table TD after appending to it columns dinf and dsup, representing the boundaries of the small or large time-window (according to extension being zero or not) for each datapoint
}
\description{
Add prepatent period to data table
}
\examples{
set.seed(42)
D<- round(runif(n = 10,min = 1,max = 30))
DT<- data.table::data.table(date=as.Date(D,origin = "1991-03-30"))
add.prepatent(DT,mean.lag = 10, duration = 12,extension = 5,unit = "day")
}

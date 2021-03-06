% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random.dates.R
\name{random.dates}
\alias{random.dates}
\title{Random dates}
\usage{
random.dates(
  Obs,
  duration = 4,
  extension = 2,
  unit = "day",
  mode = "first-half"
)
}
\arguments{
\item{Obs}{data.table of dyadic observations. Must contain a date column.}

\item{duration}{desired duration of the time-window, in unit.}

\item{extension}{bidirectional extension of the time-window.}

\item{unit}{unit in duration and extension are provided, with the smallest unit considered being days to perform Dates calculations. (EXTEND TO OTHER TIME RESOLUTION LATER?)}

\item{mode}{type of random draw. So far only "first-half" is implemented: draws an inferior date in the first half of the existing dates, then draw the closest date that is duration+extension (in days) away.}
}
\value{
Returns a vector of two random (closest) dates from Obs, separated apart from ~duration+extension days.
}
\description{
Randomly draw dates to bound random time-windows.
}
\examples{
set.seed(42);ID<- letters[1:10];
D.Obs<- round(runif(n = 100,min = 0,max = 60))
Obs<- data.table::data.table(date=as.Date(D.Obs,origin = "1991-03-30"),
                             id=sample(ID,100,replace = TRUE),
                             tar=sample(ID,100,replace = TRUE))[order(date)]
random.dates(Obs,10,5)
}

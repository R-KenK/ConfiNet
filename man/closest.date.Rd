% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/closest.date.R
\name{closest.date}
\alias{closest.date}
\title{Closest date}
\usage{
closest.date(DV, d)
}
\arguments{
\item{DV}{vector of dates. Format should be "YYYY-MM-DD"}

\item{d}{date to match. Format should be "YYYY-MM-DD"}
}
\value{
date d within DV
}
\description{
Provides which date is the closest in a dataset
}
\examples{
set.seed(42)
D<- round(runif(n = 10,min = 1,max = 30))
DT<- as.Date(D,origin = "1991-03-30")
closest.date(DT,"1991-03-30")
}

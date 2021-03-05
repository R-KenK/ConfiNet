# Closest date
#
# Provides which date is the closest in a dataset
#
# @param DV vector of dates. Format should be "YYYY-MM-DD"
# @param d date to match. Format should be "YYYY-MM-DD"
#
# @return date d within DV
# @export
#
# @examples
# set.seed(42)
# D<- round(runif(n = 10,min = 1,max = 30))
# DT<- as.Date(D,origin = "1991-03-30")
# closest.date(DT,"1991-03-30")
# closest.date<- function(DV,d){
#   d<- as.Date(d, origin = "1970-01-01");
#   return(as.Date(DV[which.min(abs(DV-d))]))
# }

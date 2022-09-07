library(nlmeVPC)
library(testit)


data(origdata)
data(simdata)
obj=aqrVPC(origdata,simdata)

assert("aqrVPC plot succeeded", "ggplot" %in% class(obj))

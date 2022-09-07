library(nlmeVPC)
library(testit)


data(origdata)
data(simdata)
obj=quantVPC(origdata,simdata,prob=0.5)

assert("quantVPC plot succeeded", "ggplot" %in% class(obj))

library(nlmeVPC)
library(testit)


data(origdata)
data(simdata)
obj=bootVPC(origdata,simdata)

assert("bootVPC plot succeeded", "ggplot" %in% class(obj))

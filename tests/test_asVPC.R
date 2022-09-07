library(nlmeVPC)
library(testit)


data(origdata)
data(simdata)
obj=asVPC(origdata,simdata,type="CI",N_hist=3,weight_method="distance")

assert("vpc plot succeeded", "ggplot" %in% class(obj))

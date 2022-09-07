library(nlmeVPC)
library(testit)


data(origdata)
data(simdata)
obj=VPCgraph(origdata,simdata,type="CI",X_name="TIME",Y_name="DV")

assert("VPCgraph plot succeeded", "ggplot" %in% class(obj))

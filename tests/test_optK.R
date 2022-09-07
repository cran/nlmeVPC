library(nlmeVPC)
library(testit)


data(origdata)
result = optK(origdata$TIME)

assert("Find correct optK", result$K == 8)

\name{origdata}
\alias{origdata}
\docType{data}
\title{
Pharmacokinetics of Theophylline with a different schedule of time. 
}
\description{
The simulated Theoph data frame has 132 rows and 3 columns of data from an 
experiment on the pharmacokinetics of theophylline.
}

\arguments{
    \item{ID}{An ordered factor with levels 1, ..., 12 identifying the subject
    on whom the observation was made. The ordering is by increasing the maximum
    concentration of theophylline observed.}
    \item{TIME}{Time since drug administration when the sample was drawn (hr).}
    \item{DV}{Theophylline concentration in the sample (mg/L).}
}

\examples{
data(origdata)
dim(origdata)
}
\keyword{datasets}

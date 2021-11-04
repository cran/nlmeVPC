#' @title  Find the best cutoff values of binning for the visual predictive checks
#' @description Find the best cutoff values for a given number of bins by various
#' rules.
#' @usage FindBestCut(X,
#'             K,
#'             beta=0.2,...)
#' @param X A numeric vector to divide into K bins.
#' @param K Number of bins.
#' @param beta Additional parameter in the penalty. 
#' For more detailed explanation, see reference. 
#' @param ... Arguments to be passed to methods.
#' @return The best cutoff values to make K bins using X and
#'         the minimum within sums of square values for the binning
#' @references Lavielle, M. and Bleakley, K.
#'             (2011). Automatic data binning for improved visual 
#'                     diagnosis of pharmacometric models. 
#'                     Journal of pharmacokinetics and pharmacodynamics, 38(6), 861-871.
#' @references VPC automatic binning algorithm in PsN 5.0.0 manual. 
#' @export
#' @examples
#' data(origdata)
#' FindBestCut(origdata$TIME,K=10)

FindBestCut = function(X,
                       K,
                       beta=0.2,...){
  origX = X
  X_n = table(X)
  X = sort(as.numeric(names(X_n)))
  
  N = length(X)
  if(N<=K){
      nest = 2:K
      CUT = (X[nest]+X[nest-1])/2
      return(list(cutoffs=CUT,WSS = 0))
  } 
  II = matrix(0,ncol=N,nrow=K)
  ntrans = matrix(0,ncol=K,nrow=(N-1))
  for(L in 1:(N-K+1)){
    tX= X[1:L]
    tX_n = X_n[1:L]
    tX= rep(tX,tX_n)
    II[1,L] = ifelse(length(tX)<=1,0,length(tX)*stats::var(tX)^beta)
  }

  for(k in 1:(K-1)){
    maxW = 0
    if(k < (K-1)){
      for(L in k:(N-K+k)){
        J = rep(10000,L)
        for(n in k:L){
          tX= X[(n+1):(L+1)]
          tX_n = X_n[(n+1):(L+1)]
          tX= rep(tX,tX_n)
          Del = ifelse(length(tX)<=1,0,length(tX)*stats::var(tX)^beta)
          J[n+1] = II[k,n]+Del
        }
        II[k+1,L+1] = min(J[1:(L+1)],na.rm=TRUE)
        ntrans[L,k] = which.min(J[1:(L+1)])

      }
    } else{
      L = N-1
      J = rep(10000,L)
      for(n in k:(N-1)){
        tX= X[(n+1):N]
        tX_n = X_n[(n+1):N]
        tX= rep(tX,tX_n)
        Del = ifelse(length(tX)<=1,0,length(tX)*stats::var(tX)^beta)
        J[n+1] = II[k,n]+Del
      }
      Imin = min(J[1:N])
      ntrans[N-1,k] = which.min(J[1:N])

    }
  }
  temp = ntrans[N-1,K-1]
  nest=temp
  for(i in (K-2):1){
    temp = ntrans[temp-2,i]
    nest = c(nest,temp)
  }
  nest = sort(nest)
  CUT = (X[nest]+X[nest-1])/2
  return(list(cutoffs=CUT,WSS = Imin))
}
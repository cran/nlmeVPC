#' @title Find the optimal number of bins
#' @description This function automatically finds the optimal number of bins using dynamic programming. 
#' @usage optK(X,
#'      Kmethod = "cluster",
#'      maxK = 10,
#'      beta = 0.2,
#'      lambda = 0.3,
#'      R = 4,
#'      C1 = 2.5,
#'      C2 = 7.8, ...)
#' @param Kmethod The way to calculate the penalty in automatic binning."cluster" or "kernel". 
#' @param maxK The maximum number of bins.
#' @param beta Additional parameter for automatic binning. For more detailed explanation, see reference.
#' @param lambda Additional parameter for automatic binning. For more detailed explanation, see reference.
#' @param R Additional parameter for automatic binning. For more detailed explanation, see reference. 
#' @param C1 Additional parameter for automatic binning. For more detailed explanation, see reference.
#' @param C2 Additional parameter for automatic binning. For more detailed explanation, see reference.
#' @param X Numeric vector corresponding to Y.
#' @param ... Arguments to be passed to methods.
#' @return The optimal number of bins, the result of binning, and the summary of
#'         binning including the penalty values up to the maximum number of bins
#'         are returned.
#' @references Lavielle, M. and Bleakley, K.
#'             (2011). Automatic data binning for improved visual 
#'                     diagnosis of pharmacometric models. 
#'                     Journal of pharmacokinetics and pharmacodynamics, 38(6), 861-871.
#' @export
#' @examples
#' \donttest{
#' data(origdata)
#' optK(origdata$TIME)
#' }

optK <- function(X,
                 Kmethod = "cluster",
                 maxK = 10,
                 beta = 0.2,
                 lambda = 0.3,
                 R = 4,
                 C1 = 2.5,
                 C2 = 7.8, ...){

   X <- X[!is.na(X)]
   N <- length(X)
   maxK <- min(c(maxK,round(length(X)/10),length(unique(X))))
   keep.tot <- matrix(0,ncol=3,nrow=(maxK-2))
   for(k in 3:maxK){
      bestCUT <- FindBestCut(X,k,...)
      if(bestCUT$WSS!=0){
        time_bin <- makeCOVbin(X,k,cutoffs=bestCUT$cutoffs)
        W.k <- (tapply(X,time_bin$COV_bin,function(x)
                   ifelse(length(x)<=1,0,stats::var(x))))*
                               (table(time_bin$COV_bin)-1)
        W <- sum(W.k)
        B.k <- ((tapply(X,time_bin$COV_bin,mean)-mean(X))^2)*
                               (table(time_bin$COV_bin))
        B <- sum(B.k)      
        Phi <- sum(unlist(lapply(bestCUT$cutoffs,function(xx)
                 sum(tapply(X,time_bin$COV_bin,function(x) {
                 nk <-length(x)
                 hk <- ifelse(nk<=1 | stats::sd(x)==0,1,stats::sd(x)/(nk^0.2))
                 hk <- ifelse(nk>1 & stats::sd(x)!=0 & 
                                timeDate::kurtosis(x)< C1,hk/R,hk)
                 z <- (xx-x)/hk
                 sum(exp(-z*z/2)/sqrt(2*pi))/(nk*hk)})))))
       alpha <- C2 * max(W.k)
       Kopt <- W+alpha*Phi
       BWratio <- ifelse(W==0,0,(B/(k-1))/(W/(N-k)))
       KBratio <- ifelse(Kopt==0,0,BWratio/Kopt)
       Jopt <- ifelse(bestCUT$WSS==0,0,log(bestCUT$WSS))+lambda*beta*k
     } else{
        Jopt <- KBratio <- 0
     }
     keep.tot[k-2,] <- c(k,Jopt,KBratio)
   }
   colnames(keep.tot) <- c("k","Jopt","KBratio")
   keep.tot <- data.frame(keep.tot)
   if(Kmethod=="cluster"){
      K <- min(keep.tot$k[which(keep.tot$Jopt==min(keep.tot$Jopt,na.rm=TRUE))])
   } else if(Kmethod=="kernel"){
      K <- min(keep.tot$k[which(keep.tot$KBratio==
                                  min(keep.tot$KBratio,na.rm=TRUE))])      
   }
   CUT <- FindBestCut(X,K,...)$cutoffs
   time_bin <- makeCOVbin(X,K,cutoffs=CUT)
   return(list(K=K,time_bin=time_bin,cutoffs = CUT,penaltyC = keep.tot))
}


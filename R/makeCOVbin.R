#' @title Discretise numeric data into categorical variable
#' @description Discretise numeric value into a categorical variable using  
#'  the user-defined breaks. If cutoffs and the number of bins (K) are NULL, 
#'   find the best number of bins using optK function and find the best 
#'   cutoff values using FindBestCut function. 
#' @usage makeCOVbin(X,
#'            K,
#'            cutoffs,
#'            adjust0bin = TRUE, ...)
#' @param X A numeric vector corresponding to Y.
#' @param K Number of bins.
#' @param cutoffs A numeric vector of two or more unique cut points.
#' @param adjust0bin Adjust bin with 0 observation if TRUE.
#' @param ... Arguments to be passed to methods.
#' @return The result of binning and the summary of the binning results
#' @references Lavielle, M. and Bleakley, K.
#'             (2011). Automatic data binning for improved visual 
#'                     diagnosis of pharmacometric models. 
#'                     Journal of pharmacokinetics and pharmacodynamics, 38(6), 861-871.
#' @export
#' @examples
#' data(origdata)
#' CUT = FindBestCut(origdata$TIME,8)$cutoffs
#' makeCOVbin(origdata$TIME,K=8,cutoffs=CUT)

makeCOVbin <- function(X,
                       K,
                       cutoffs,
                       adjust0bin = TRUE, ...){
   X <- X[!is.na(X)]
   data_temp <- data.frame(X=X)
   range_temp <- range(X)
   cutoffsOrig <- cutoffs
   if(min(cutoffs)>=range_temp[1]){
         cutoffs <- c(range_temp[1]-stats::sd(X)*0.1,cutoffs)
   }
   if(max(cutoffs)<=range_temp[2]){
         cutoffs <- c(cutoffs,range_temp[2]+stats::sd(X)*0.1)
   }

   cut_temp <- cut(X, breaks=cutoffs)
   if(adjust0bin){
     if(sum(table(cut_temp)==0)!=0){
        cutoffs <- cutoffs[-which(table(cut_temp)==0)]
        cut_temp <- cut(X, breaks=cutoffs)
     }
   }

   temp <- tapply(X,cut_temp,stats::median)
   tab <- data.frame(cut_temp =factor(names(temp),levels=levels(cut_temp)),
                     mid_COV = round(temp,2))
   
   cutoffs[c(1,length(cutoffs))] <- range_temp
   LU_temp <- cbind(cutoffs[-length(cutoffs)],cutoffs[-1])
   colnames(LU_temp) <- c("lower_COV","upper_COV")
   mid_LU <- apply(LU_temp,1,mean)
   tab <- data.frame(tab,n_bin=c(table(cut_temp)),LU_temp,mid_LU=mid_LU)
   return(list(COV_bin=cut_temp,COVbin_summary=tab))
}

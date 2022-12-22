#' @title Find quantiles of the original data.
#' @usage findQuantile(Y,
#'             X,
#'             X_bin,
#'             probs=c(0.1,0.5,0.9),...)
#' @param Y A numeric vector whose sample quantiles are wanted.
#' @param X A numeric vector corresponding to Y.
#' @param X_bin Binning result from makeCOVbin function.
#' @param probs A numeric vector of probabilities.
#' @param ... Arguments to be passed to methods.
#' @return quantiles of Y using X_bin
#' @export
#' @examples
#' data(origdata)
#' CUT = FindBestCut(origdata$TIME,8)$cutoffs
#' time_bin = makeCOVbin(origdata$TIME,K=8,cutoffs = CUT)
#' findQuantile(origdata$DV,origdata$TIME,X_bin=time_bin)

findQuantile <-function(Y,
                        X,
                        X_bin,
                        probs = c(0.1,0.5,0.9), ...){

      dataT <- data.frame(X,Y,Xbin=X_bin$COV_bin)
      temp <- matrix(0,ncol=length(probs),nrow=nrow(X_bin$COVbin_summary))
      for(i in 1:length(probs)){
        temp[,i] <- tapply(dataT$Y,dataT$Xbin,
                           function(x) stats::quantile(x,prob=probs[i],
                                                na.rm=TRUE,names=FALSE))
      }
      nameT <- strsplit(rapply(list(probs),sprintf,
                              fmt="%.2f",how="replace")[[1]],"\\.")
      nameT <- unlist(lapply(nameT,function(x) x[2]))
      colnames(temp) <- paste0("Q",nameT,"th")
      DV_q <- data.frame(X_bin$COVbin_summary,temp)
      return(DV_q)
}

#' @title Find quantiles of the simulated data.
#' @importFrom stats  quantile
#' @usage findSIMQuantile(sim_data,
#'                X,
#'                X_bin,
#'                probs = c(0.1,0.5,0.9),
#'                conf.level = 0.95,
#'                approx = FALSE, ...)
#' @param sim_data A matrix of simulated data with only Y values collected.
#' @param X A numeric vector corresponding to Y.
#' @param X_bin Binning result from makeCOVbin function.
#' @param probs A numeric vector of probabilities.
#' @param conf.level Confidence level of the interval.
#' @param approx Arguments to be passed to methods
#' @param ... Arguments to be passed to methods
#' @return quantiles of sim_data using X_bin
#' @export
#' @examples
#' data(origdata)
#' data(simdata)
#' CUT = FindBestCut(origdata$TIME,8)$cutoffs
#' time_bin = makeCOVbin(origdata$TIME,K=8,cutoffs = CUT)
#' findSIMQuantile(simdata,origdata$TIME,X_bin=time_bin)

findSIMQuantile <- function(sim_data,
                            X,
                            X_bin,
                            probs = c(0.1,0.5,0.9),
                            conf.level = 0.95,
                            approx = FALSE, ...){
   
   tempSIMQ <- findSIMQ(sim_data,X,X_bin,
                       probs,conf.level,approx)
   nameT <- strsplit(rapply(list(probs),sprintf,
                           fmt="%.2f",how="replace")[[1]],"\\.")
   nameT <- paste0("Q",unlist(lapply(nameT,function(x) x[2])),"th")
   nameT <- strsplit(rapply(list(tempSIMQ$quantile),sprintf,
                           fmt="%.2f",how="replace")[[1]],"\\.")
   nameT <- paste0("Q",unlist(lapply(nameT,function(x) x[2])),"th")
   
   SIM_quantile <- data.frame(cut_temp=
                                X_bin$COVbin_summary$cut_temp[tempSIMQ$cuttemp],
                              quantile=nameT,
                              tempSIMQ$quantileV)
               
   nameT <- as.character(c((1-conf.level)/2,
                          0.5, 
                          1-(1-conf.level)/2)*100)
   colnames(SIM_quantile)[-c(1:2)] <- paste0(nameT,"%")
   SIM_quantile <- SIM_quantile
   return(SIM_quantile)
}

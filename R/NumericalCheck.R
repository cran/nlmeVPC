#' @title The numerical predictive checks
#' @description This function calculates the numerical predictive checks for 
#' each prediction level. For a given level of prediction, 
#' the predicted interval is calculated using the simulated data, 
#' and the number of observed data below the predicted interval is counted. 
#' The expected number of points below the predicted interval 
#' is also calculated and compared to the observed number.
#' @usage NumericalCheck(orig_data,
#'                sim_data,
#'                N_xbin = NULL,
#'                pred.level = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
#'                conf.level = 0.95,
#'                X_name = "TIME",
#'                Y_name = "DV",
#'                MissingDV = NULL,
#'                Kmethod = "cluster",                
#'                maxK = NULL,
#'                beta = 0.2,
#'                lambda = 0.3,
#'                R = 4,
#'                C1 = 2.5,
#'                C2 = 7.8, ...)
#' @param orig_data A data frame of original data with X and Y variable.
#' @param sim_data A matrix of simulated data with only Y values collected.
#' @param N_xbin Number of bins in X variable. If NULL, optimal number of bins are automatically calcuated using optK function.
#' @param pred.level Numeric vector of probabilities.
#' @param conf.level Confidence level of the interval.
#' @param X_name Name of X variable in orig_data (usually "TIME" in pharmacokinetic data).
#' @param Y_name Name of Y variable in orig_data (usually "DV" in pharmacokinetic data).
#' @param MissingDV Name of missing indicator variable in orig_data, which have value 1 if missing, value 0 otherwise. (usually "MDV" in pharmacokinetic data).
#' @param maxK The maximum number of bins.
#' @param Kmethod The way to calculate the penalty in automatic binning."cluster" or "kernel". 
#' @param beta Additional parameter for automatic binning, used in optK function.  
#' @param lambda Additional parameter for automatic binning, used in optK function.  
#' @param R Additional parameter for automatic binning, used in optK function. 
#' @param C1 Additional parameter for automatic binning, used in optK function.  
#' @param C2 Additional parameter for automatic binning, used in optK function. 
#' @param ... Arguments to be passed to methods.
#' @return The result of numerical predictive check
#' @references Holford N, & Karlsson M. (2008). "A tutorial on visual predictive checks,
#'  abstr 1434." Annual Meeting of the Populations Approach Group in Europe. www.page-meeting.org. 2008.
#' @references Harling, Uekcert, K. 2018. VPC and NPC User Guide. ICON plc. 
#' @references https://github.com/UUPharmacometrics/PsN/releases/download/4.9.0/vpc_npc_userguide.pdf.
#' @export
#' @examples
#' \donttest{
#' data(origdata)
#' data(simdata)
#' NumericalCheck(origdata,simdata,N_xbin=8)$NPC
#' }

NumericalCheck <- function(orig_data,
                           sim_data,
                           N_xbin = NULL,
                           pred.level = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                           conf.level = 0.95,
                           X_name = "TIME",
                           Y_name = "DV",
                           MissingDV = NULL,
                           Kmethod = "cluster",                
                           maxK = NULL,
                           beta = 0.2,
                           lambda = 0.3,
                           R = 4,
                           C1 = 2.5,
                           C2 = 7.8, ...){

   sel.id <- !is.na(orig_data[,Y_name]) 
   if(!is.null(MissingDV))
      sel.id <- sel.id & orig_data[,MissingDV]==0
   sim_data <- sim_data[sel.id,]  
   orig_data <- orig_data[sel.id,]
   
   N_sim <- ncol(sim_data)
   if(is.null(N_xbin))
      N_xbin <- optK(orig_data[,X_name],...)$K
   
   plot_data <- data.frame(X=orig_data[,X_name],Y=orig_data[,Y_name])

   if(N_xbin < length(unique(plot_data$X))){
      CUT <- FindBestCut(plot_data$X,K=N_xbin,...)$cutoffs
      time_bin <- makeCOVbin(plot_data$X,N_xbin,cutoffs=CUT,...)
   } else{
      cut_temp.id <- as.numeric(sort(unique(plot_data$X)))
      temp_diff <- diff(cut_temp.id)/2
      CUT <- c(cut_temp.id[1]-temp_diff[1],
             cut_temp.id[-length(cut_temp.id)]+temp_diff,
             cut_temp.id[length(cut_temp.id)]+
                temp_diff[length(temp_diff)])
      CUT <- CUT[CUT<max(plot_data$X) & CUT>min(plot_data$X)]
      time_bin <- makeCOVbin(plot_data$X,N_xbin,cutoffs=CUT,...)
   }

   plot_data <- data.frame(plot_data, cut_temp = time_bin$COV_bin,
                           tempID = 1:nrow(plot_data))
   keepAll <- matrix(0,nrow=length(pred.level),ncol=9)
   probkeep <- NULL
   keepAll2 <- data.frame(belowPI = rep(0,length(pred.level)*N_xbin),
                          abovePI =rep(0,length(pred.level)*N_xbin),
                          n = rep(0,length(pred.level)*N_xbin),
                          pred.level =rep(0,length(pred.level)*N_xbin),
                          belowE = rep(0,length(pred.level)*N_xbin),
                          aboveE = rep(0,length(pred.level)*N_xbin),
                          cut_temp = rep("0",length(pred.level)*N_xbin))
   CIkeep <- matrix(0,nrow=N_xbin,ncol=2*length(pred.level))

   for(i in 1:length(pred.level)){
      probs <- c((1-pred.level[i])/2,1-(1-pred.level[i])/2)
      probkeep <- c(probkeep,probs)
      simAll <- data.frame(TIME =rep(plot_data$X,ncol(sim_data)),
                          DV = c(sim_data),
                          cut_temp = rep(time_bin$COV_bin,ncol(sim_data)))
      temp_simQ <- matrix(unlist(tapply(simAll$DV,simAll$cut_temp,
                                        quantile,probs)),ncol=2,byrow=TRUE)
      colnames(temp_simQ) <- c("LB","UB")
      temp_simQ <- data.frame(cut_temp = time_bin$COVbin_summary$cut_temp,
                              temp_simQ)      

      temp_data <- merge(plot_data,temp_simQ,by="cut_temp",all.x=TRUE,sort=FALSE)
      temp_data <- temp_data[sort.list(temp_data$tempID),]
       
      temp_data$belowPIflag <- temp_data$Y<temp_data$LB
      temp_data$abovePIflag <- temp_data$Y>temp_data$UB
      temp1 <- data.frame(belowPI = tapply(temp_data$belowPIflag,
                                          temp_data$cut_temp,sum),
                          abovePI = tapply(temp_data$abovePIflag,
                                          temp_data$cut_temp,sum),
                          n = tapply(temp_data$belowPIflag,
                                     temp_data$cut_temp,length),
                         pred.level = pred.level[i])
      temp1$belowE <- temp1$n*(1-temp1$pred.level)/2
      temp1$aboveE <-temp1$belowE
      temp1$cut_temp <- rownames(temp1)
      keepAll2$belowPI[((i-1)*N_xbin+1):(i*N_xbin)] <- temp1$belowPI
      keepAll2$abovePI[((i-1)*N_xbin+1):(i*N_xbin)] <- temp1$abovePI
      keepAll2$n[((i-1)*N_xbin+1):(i*N_xbin)] <- temp1$n
      keepAll2$pred.level[((i-1)*N_xbin+1):(i*N_xbin)] <- temp1$pred.level
      keepAll2$belowE[((i-1)*N_xbin+1):(i*N_xbin)] <- temp1$belowE
      keepAll2$aboveE[((i-1)*N_xbin+1):(i*N_xbin)] <- temp1$aboveE
      keepAll2$cut_temp[((i-1)*N_xbin+1):(i*N_xbin)] <- temp1$cut_temp
      temp1 <- apply(temp1[,c(3,5,6,1,2)],2,sum)

      tempAll <- t(apply(sim_data,2,
                        function(x) return(c(sum(x<temp_data$LB),
                                             sum(x>temp_data$UB)))))
      temp2 <- c(apply(tempAll,2,
                    function(x) stats::quantile(x,
                                    probs=c((1-conf.level)/2,
                                            1-(1-conf.level)/2))))
      names(temp2)<- c("belowPI_L","belowPI_U","abovePI_L","abovePI_U")
      keepAll[i,] <- c(temp1,temp2)
      colnames(temp_simQ)[2:3] <- paste0(as.character(pred.level[i]*100),
                                        "%",colnames(temp_simQ)[2:3])
  
      CIkeep[,(2*i-1)] <- temp_simQ[,2]
      CIkeep[,2*i] <- temp_simQ[,3]      
   }
   CIkeep <- data.frame(X_bin=time_bin$COVbin_summary$cut_temp,CIkeep)
   colnames(keepAll) <- c("n","belowE","aboveE","belowPI","abovePI",
                        "belowPI_L","belowPI_U","abovePI_L","abovePI_U")
   keepAll <- data.frame(keepAll)
   n <- keepAll$n[1]
   NPC <- data.frame(PI = pred.level*100,
                     n = keepAll$n,
                     belowE = keepAll$belowE,
                     belowPI = keepAll$belowPI,
                     belowPI_p = keepAll$belowPI/n*100,
                     belowPI_L = keepAll$belowPI_L/n*100,
                     belowPI_U = keepAll$belowPI_U/n*100,
                     aboveE = keepAll$aboveE,
                     abovePI = keepAll$abovePI,
                     abovePI_p = keepAll$abovePI/n*100,
                     abovePI_L = keepAll$abovePI_L/n*100,
                     abovePI_U = keepAll$abovePI_U/n*100)

   colnames(NPC) <- c("PI","n",
                      "Expected points below PI",
                      "points below PI","points below PI(%)",
                      paste0(round(conf.level*100),"%CIBelowFrom(%)"),
                      paste0(round(conf.level*100),"%CIBelowTo(%)"),
                      "Expected points above PI",
                      "points above PI","points above PI(%)",
                      paste0(round(conf.level*100),"%CIAboveFrom(%)"),
                      paste0(round(conf.level*100),"%CIAboveTo(%)")) 

   return(list(NPC = NPC,
               DV_point = plot_data,
               sim_quant = CIkeep,
               NPCcut = keepAll2))
}








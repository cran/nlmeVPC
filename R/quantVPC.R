#' @title the quantified visual predictive check plot (QVPC)
#' @description  The quantified visual predictive check visually represents 
#' actual and unavailable observations around predicted medians, regardless 
#' of the density or shape of the observed data distribution,  through the 
#' form of a percent. 
#' @import ggplot2 quantreg optimx
#' @usage quantVPC(orig_data,
#'          sim_data,
#'          N_xbin = NULL,
#'          prob = 0.5,
#'          X_name = "TIME",
#'          Y_name = "DV",
#'          MissingDV = NULL,
#'          Kmethod = "cluster",                
#'          maxK = NULL,
#'          beta = 0.2,
#'          lambda = 0.3,
#'          R = 4,
#'          C1 = 2.5,
#'          C2 = 7.8, ...)
#' @param orig_data A data frame of original data with X and Y variable.
#' @param sim_data A matrix of simulated data with only Y values collected.
#' @param N_xbin Number of bins in X variable. If NULL, optimal number of bins are automatically calcuated using optK function.
#' @param prob Scalar of probability.
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
#' @return quantVPC plot 
#' @references Post, T.M., et al. (2008)
#' Extensions to the visual predictive check for facilitate model performance
#' evaluation, Journal of pharmacokinetics and pharmacodynamics,
#' 35(2), 185-202
#' @export
#' @examples
#' \donttest{
#' data(origdata)
#' data(simdata)
#' quantVPC(origdata,simdata,prob=0.5,N_xbin=8)
#' }

quantVPC<-function(orig_data,
                   sim_data,
                   N_xbin = NULL,
                   prob = 0.5,
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
  
   main_title <- "Quantified VPC"
   unavailable <- TIME <- percent <- type <- NA
   sel.id <- !is.na(orig_data[,Y_name]) 
   if(!is.null(MissingDV))
     sel.id <- sel.id & orig_data[,MissingDV]==0
   sim_data <- sim_data[sel.id,]  
   orig_data <- orig_data[sel.id,]
   plot_data <- data.frame(X=orig_data[,X_name],Y=orig_data[,Y_name])
   if(length(prob)>1){
      prob <- prob[1]
      warning(paste("prob should be a scalar. Use the first element",
                    prob[1]))
   }

  if(is.null(N_xbin))
    N_xbin <- optK(orig_data[,X_name],...)$K


   time_bin <- makeCOVbin(plot_data$X,
                        K = N_xbin,
                        cutoffs=FindBestCut(plot_data$X,N_xbin)$cutoffs,...)

   temp.sim.Q <- findQuantile(plot_data$Y,plot_data$X,time_bin,probs=prob, ...)

   temp.sim.Q2 <- findSIMQuantile(sim_data,plot_data$X,time_bin,probs=prob,...)
   temp.sim.Q$Q50th <- unlist(temp.sim.Q2[,"50%"])
   QVPC.temp <- matrix(0,nrow=nrow(temp.sim.Q),ncol=7)
   expN <- nrow(plot_data)/N_xbin

   for(i in 1:nrow(temp.sim.Q)){
     time.cut <- temp.sim.Q[i,1]
     N <- temp.sim.Q[i,3]
     midtime <- temp.sim.Q[i,6]
     Mt <- temp.sim.Q[i,7]
     temp <- plot_data[which(time_bin$COV_bin==time.cut),]
     above <- (sum(temp$Y>Mt)+sum(temp$Y==Mt)/2)/N*100
     below <- (sum(temp$Y<Mt)+sum(temp$Y==Mt)/2)/N*100
     unavail <- 100-(above+below)
     obsmid <- prob*100+unavail/2
     QVPC.temp[i,] <- c(midtime,N,Mt,obsmid,above,below,unavail)
   }

  QVPC.temp <- as.data.frame(QVPC.temp)
  colnames(QVPC.temp) <- c("TIME","N","Mt","obsmid","above","below","unavailable")
    
  plot.data <- data.frame(TIME = rep(QVPC.temp$TIME,3),
                         N = rep(QVPC.temp$N,3),
                         type = rep(c("above","below","unavailable"),
                                    each=nrow(QVPC.temp)),
                         obsmid = rep(QVPC.temp$obsmid,3),
                         percent = unlist(QVPC.temp[,5:7]))
   
  ggplot(data=plot.data) +
       geom_bar(aes(as.factor(TIME),percent,fill=type),
                stat="identity",position = "stack") +
       geom_point(aes(as.factor(TIME),obsmid),col="white") +
       scale_y_continuous(expand=c(0,0)) +
       scale_fill_manual("Type",values = c("grey45","grey20","grey75")) +
       theme_bw() +
       labs(x=X_name,y="Percentage Observations(%)",title=main_title) +
       scale_x_discrete(labels=levels(temp.sim.Q$cut_temp)) +
       theme(axis.text.x=element_text(angle=45,hjust=1))

}





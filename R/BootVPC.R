#' @title the bootstrap visual predicrive checks. 
#' @description This function draws the visual predictive check plot with bootstrapped data.
#'   It compares
#'  the distribution of quantiles obtained from the bootstrapped data to 
#'  the distribution of quantiles from simulated data drawn from the fitted model. 
#' @import ggplot2 quantreg optimx
#' @title Visual predictive checks
#' @usage bootVPC(orig_data,
#'         sim_data,
#'         B = 1000,                  
#'         N_xbin = NULL,
#'         conf.level = 0.95,
#'         X_name = "TIME",
#'         Y_name = "DV",
#'         subject_name = "ID",                  
#'         MissingDV = NULL,
#'         DV_point = TRUE,                  
#'         plot_caption = TRUE,                  
#'         plot_flag = TRUE,
#'         linesize = 0.7,
#'         pointsize = 0.7,
#'         Kmethod = "cluster",                
#'         maxK = NULL,
#'         beta = 0.2,
#'         lambda = 0.3,
#'         R = 4,                
#'         C1 = 2.5,
#'         C2 = 7.8, ...)
#' @param orig_data A data frame of original data with X and Y variable.
#' @param sim_data A matrix of simulated data with only Y values collected.
#' @param B Number of bootstrap samples.
#' @param N_xbin Number of bins in X variable. If NULL, optimal number of bins are automatically calcuated using optK function.
#' @param conf.level Confidence level of the interval.
#' @param X_name Name of X variable in orig_data (usually "TIME" in pharmacokinetic data).
#' @param Y_name Name of Y variable in orig_data (usually "DV" in pharmacokinetic data).
#' @param subject_name Name of subject variable in orig_data (usually "ID" in pharmacokinetic data).
#' @param MissingDV Name of missing indicator variable in orig_data, which have value 1 if missing, value 0 otherwise. (usually "MDV" in pharmacokinetic data).
#' @param DV_point Draw point (X, Y) in the plot if TRUE; omit if FALSE.
#' @param plot_caption Put caption with additional information if TRUE; omit if FALSE.
#' @param plot_flag Draw plot if TRUE; generate data for drawing plot if FALSE.
#' @param linesize Size of line in the plot.
#' @param pointsize Size of point in the plot.
#' @param maxK The maximum number of bins.
#' @param Kmethod The way to calculate the penalty in automatic binning."cluster" or "kernel". 
#' @param beta Additional parameter for automatic binning, used in optK function.  
#' @param lambda Additional parameter for automatic binning, used in optK function.  
#' @param R Additional parameter for automatic binning, used in optK function. 
#' @param C1 Additional parameter for automatic binning, used in optK function.  
#' @param C2 Additional parameter for automatic binning, used in optK function. 
#' @param ... Arguments to be passed to methods.
#' @return bootVPC plot or the values to draw bootVPC plot.
#' @references Post, T. M., et al. (2008)
#' Extensions to the visual predictive check for facilitate model performance
#' evaluation, Journal of pharmacokinetics and pharmacodynamics,
#' 35(2), 185-202
#' @export
#' @examples
#' \donttest{
#' data(origdata)
#' data(simdata)
#' bootVPC(origdata,simdata,N_xbin=8)
#' }

bootVPC <- function(orig_data,
                    sim_data,
                    B = 1000,                  
                    N_xbin = NULL,
                    conf.level = 0.95,
                    X_name = "TIME",
                    Y_name = "DV",
                    subject_name = "ID",                  
                    MissingDV = NULL,
                    DV_point = TRUE,                  
                    plot_caption = TRUE,                  
                    plot_flag = TRUE,
                    linesize = 0.7,
                    pointsize = 0.7,                
                    Kmethod = "cluster",                
                    maxK = NULL,
                    beta = 0.2,
                    lambda = 0.3,
                    R = 4,
                    C1 = 2.5,
                    C2 = 7.8, ...){
   main_title <- "Bootstrap VPC"
   quantile <- `50%` <- cut_temp <- mid_LU <- X <- Y <- LB <- UB <- Mid <- NULL
   sel.id <- !is.na(orig_data[,Y_name]) 
   if(!is.null(MissingDV))
     sel.id <- sel.id & orig_data[,MissingDV]==0
   sim_data <- sim_data[sel.id,]  
   orig_data <- orig_data[sel.id,]
   
   probs <- c( (1-conf.level)/2, 0.5, 1-(1-conf.level)/2)

   plot_data <- data.frame(X=orig_data[,X_name],Y=orig_data[,Y_name],
                         Subject = orig_data[,subject_name])

   if(is.null(N_xbin))
      N_xbin <- optK(orig_data[,X_name],...)$K

   caption <- paste(paste("No. of bins = ",N_xbin," / ",
                          "Percentile=",round(probs[1]*100),"th,",
                          round(probs[2]*100),"th,",
                          round(probs[3]*100),"th",sep=""),"/",
                          round(conf.level*100),"% CI")

   time_bin <- makeCOVbin(plot_data$X,K = N_xbin,
                          cutoffs=FindBestCut(plot_data$X,N_xbin)$cutoffs,...)
   Y_min <- min(plot_data$Y,na.rm=T)
   Y_max <- max(plot_data$Y,na.rm=T)

   bootTemp <- matrix(0,nrow=B,ncol=N_xbin)
   IDlist <- sort(unique(plot_data$Subject))
   time.cut <- time_bin$COV_bin
   plot_data <- data.frame(plot_data,time.cut)
   for (j in 1:B){
     IDboot <- sample(IDlist,replace=TRUE)
     BootSamp <- plot_data[unlist(lapply(IDboot,function(x)
        (1:nrow(plot_data))[plot_data$Subject==x])),]
     bootTemp[j,] <- tapply(BootSamp$Y,BootSamp$time.cut,
                             function(x) stats::median(x,na.rm=TRUE))
   }

   BootTemp <- t(apply(bootTemp,2,
                    function(x)stats::quantile(x,probs,na.rm=TRUE)))
   Boot_CI <- data.frame(time_bin$COVbin_summary[,c(1,6)],
                       quantile="Q50th",BootTemp)
   colnames(Boot_CI) <- c("X_bin","X","quantile",colnames(BootTemp))


   SIM_quant <- findSIMQuantile(sim_data,plot_data$X,time_bin,probs=probs, ...)
   
   sel.id <- (1:(nrow(SIM_quant)/3))*3
   TEMP <- data.frame(cut_temp = SIM_quant$cut_temp[sel.id],
                     X1 = SIM_quant[sel.id-2,"50%"],
                     X2 = SIM_quant[sel.id-1,"50%"],
                     X3 = SIM_quant[sel.id,"50%"])
   colnames(TEMP)[-1] <- SIM_quant$quantile[1:3]
   SIM_quant <- merge(TEMP,time_bin$COVbin_summary[,c(1,6)],
                     by="cut_temp",all.x=TRUE,sort=FALSE)
   SIM_quant <- SIM_quant[,c(1,ncol(SIM_quant),2:(ncol(SIM_quant)-1))]
   colnames(SIM_quant)[1:2] <- c("X_bin","X")

   if(plot_flag){
      Y_min <- as.numeric(min(c(unlist(SIM_quant[,5]),plot_data$Y),na.rm=T))
      Y_max <- as.numeric(max(c(unlist(SIM_quant[,3]),plot_data$Y),na.rm=T))

      Ptemp <- ggplot()+
                 theme_bw()+
                 theme(panel.grid.minor = element_blank(),
                       panel.grid.major=element_blank())#+ylim(Y_min,Y_max)
  
      if(DV_point){
        Ptemp <- Ptemp+geom_point(data = plot_data,aes(X,Y),
                              color="grey30",size=pointsize,alpha=0.5)
      } 
      D1 <- Boot_CI
      colnames(D1)[c(4:6)] <- c("LB","Mid","UB")
      D2 <- SIM_quant
      colnames(D2)[c(3:5)] <- c("LB","Mid","UB")
      if(D1$X[1]>min(plot_data$X)){
        temp <- D1[1,]; temp$X <- min(plot_data$X)
        D1 <- rbind(temp,D1)
        temp <- D2[1,]; temp$X <- min(plot_data$X)
        D2 <- rbind(temp,D2)
      }
      if(D1$X[nrow(D1)]<max(plot_data$X)){
        temp <- D1[nrow(D1),]; temp$X <- max(plot_data$X)
        D1 <- rbind(D1,temp)
        temp <- D2[nrow(D2),]; temp$X <- max(plot_data$X)
        D2 <- rbind(D2,temp)
      }       
      Ptemp <- Ptemp +
            geom_ribbon(data=D1,aes(x=X,ymin=LB,ymax = UB),
                        fill="lightpink",alpha=0.5)+
            geom_line(data = D1,aes(x=X,y=Mid),
                         color="red",size=linesize) +      
            geom_line(data = D2,aes(x=X,y=LB),linetype="dashed",
                         color="blue",size=linesize)+
            geom_line(data = D2,aes(x=X,y=Mid),
                         color="blue",size=linesize)+
            geom_line(data = D2,aes(x=X,y=UB),linetype="dashed",
                         color="blue",size=linesize)
  

     Ptemp <- Ptemp +
              theme(panel.grid.major=element_line(colour = "white"),
                    panel.grid.minor=element_line(colour="white"),
                    plot.caption=element_text(hjust=0.5))+
             labs(x=X_name,y=Y_name,title=main_title)
     Ptemp

  } else{
    DV_point <- orig_data[,c(X_name,Y_name)]
    return(list(DV_point=DV_point,
                Boot_CI=Boot_CI,
                SIM_quant=SIM_quant))
  }
}

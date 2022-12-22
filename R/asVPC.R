#' @title The average shifted visual predictive checks (asVPC)
#' @description  This function draws the average shifted visual predictive check
#' (asVPC) plot. It calculates original and simulated data percentiles using 
#' the average shifted histogram method. After calculating percentiles with 
#' bin-related or distance-related weights, draw the VPC type plot.
#'
#' @import ggplot2
#' @importFrom Hmisc wtd.quantile
#' @usage asVPC(orig_data,
#'       sim_data,
#'       type = "CI",                
#'       weight_method = "bin",                  
#'       N_xbin = NULL,
#'       N_hist = NULL,
#'       probs = c(0.1,0.5,0.9),
#'       conf.level = 0.95,
#'       X_name = "TIME",
#'       Y_name = "DV",
#'       MissingDV = NULL,
#'       DV_point = TRUE,                
#'       CIvpc_type = "line",
#'       bin_grid = TRUE,                   
#'       plot_caption = TRUE,
#'       plot_flag = TRUE,
#'       linesize = 0.7,
#'       pointsize = 0.7,
#'       captionsize = 10,
#'       Kmethod = "cluster",                
#'       maxK = NULL,
#'       beta = 0.2,
#'       lambda = 0.3,
#'       R = 4,
#'       C1 = 2.5,
#'       C2 = 7.8,...)
#' @param orig_data A data frame of original data with X and Y variable.
#' @param sim_data A matrix of simulated data with only Y values collected.
#' @param type Type of VPC graph; "CI", "percentile", or "scatter".
#' @param weight_method The way to put weight when the average shifted values are calculated. "bin" or "distance".
#' @param N_xbin Number of bins in X variable. If NULL, optimal number of bins are automatically calcuated using optK function.
#' @param N_hist The number of shifted histograms.
#' @param probs A numeric vector of probabilities.
#' @param conf.level Confidence level of the interval.
#' @param X_name Name of X variable in orig_data (usually "TIME" in pharmacokinetic data).
#' @param Y_name Name of Y variable in orig_data (usually "DV" in pharmacokinetic data).
#' @param MissingDV Name of missing indicator variable in orig_data, which have value 1 if missing, value 0 otherwise. (usually "MDV" in pharmacokinetic data).
#' @param DV_point Draw point (X, Y) in the plot if TRUE; omit if FALSE.
#' @param CIvpc_type Type of CI area in VPC graph; "line" or "segment".
#' @param bin_grid Draw grid lines for binning in X variable if TRUE; omit if FALSE.
#' @param plot_caption Put caption with additional information if TRUE; omit if FALSE.
#' @param plot_flag Draw plot if TRUE; generate data for drawing plot if FALSE.
#' @param linesize Size of line in the plot.
#' @param pointsize Size of point in the plot.
#' @param captionsize Size of caption.
#' @param maxK The maximum number of bins.
#' @param Kmethod The way to calculate the penalty in automatic binning."cluster" or "kernel". 
#' @param beta Additional parameter for automatic binning, used in optK function.  
#' @param lambda Additional parameter for automatic binning, used in optK function.  
#' @param R Additional parameter for automatic binning, used in optK function. 
#' @param C1 Additional parameter for automatic binning, used in optK function.  
#' @param C2 Additional parameter for automatic binning, used in optK function. 
#' @param ... Arguments to be passed to methods.
#' @return asVPC plot or the values to draw asVPC plot.
#' @export
#' @examples
#' \donttest{
#' data(origdata)
#' data(simdata)
#' asVPC(origdata,simdata,type="CI",N_hist=3,weight_method="distance",N_xbin=8)
#' asVPC(origdata,simdata,type="CI",N_hist=3,weight_method="bin",N_xbin=8)
#' }
#' 
asVPC <- function(orig_data,
                  sim_data,
                  type = "CI",                
                  weight_method = "bin",                  
                  N_xbin = NULL,
                  N_hist = NULL,
                  probs = c(0.1,0.5,0.9),
                  conf.level = 0.95,
                  X_name = "TIME",
                  Y_name = "DV",
                  MissingDV = NULL,
                  DV_point = TRUE,                
                  CIvpc_type = "line",
                  bin_grid = TRUE,                   
                  plot_caption = TRUE,
                  plot_flag = TRUE,
                  linesize = 0.7,
                  pointsize = 0.7,
                  captionsize = 10,
                  Kmethod = "cluster",                
                  maxK = NULL,
                  beta = 0.2,
                  lambda = 0.3,
                  R = 4,
                  C1 = 2.5,
                  C2 = 7.8, ...){
   
   main_title <- paste0("average shifted VPC : ",
                        weight_method,"-related weights")

   sel.id <- !is.na(orig_data[,Y_name]) 
   if(!is.null(MissingDV))
      sel.id <- sel.id & orig_data[,MissingDV]==0
   sim_data <- sim_data[sel.id,]  
   orig_data <- orig_data[sel.id,]
   plot_data <- data.frame(X=orig_data[,X_name],Y=orig_data[,Y_name])
   approxi <- FALSE
   if(!is.null(N_xbin))
      approxi <- ifelse(N_xbin >= length(unique(plot_data$X)) |
                        min(table(plot_data$X)<5),TRUE,FALSE)

   if(type=="scatter"){
      DV_qline <- FALSE; SIM_qline <- TRUE; SIM_qCI <- FALSE;
   } else if(type=="percentile"){
      DV_qline <- TRUE; SIM_qline <- TRUE; SIM_qCI <- FALSE;
   } else if(type=="CI"){
      DV_qline <- TRUE; SIM_qline <- FALSE; SIM_qCI <- TRUE ;
   }
   if(!plot_flag)
      DV_qline <- SIM_qline <- SIM_qCI <- TRUE
   if(approxi) {
      DV_qline <- FALSE
      DV_point <- TRUE
      SIM_qCI <- FALSE
      SIM_qline <- TRUE
   }

   if(is.null(N_xbin)){
      optK_t <- optK(orig_data[,X_name],...)
      N_xbin <- optK_t$K
   }

   if(is.null(N_hist))
     N_hist <- ifelse(N_xbin==length(unique(orig_data[,X_name])),
                    2,round(N_xbin/2))

   bintot.N <- N_xbin*N_hist
   timebin_flag <- bintot.N < length(unique(plot_data$X))
   if(timebin_flag){
      CUTtemp <- c(min(plot_data$X),
                  FindBestCut(plot_data$X,K=N_xbin,...)$cutoffs,
                  max(plot_data$X))
      CUT <- rep(CUTtemp[-length(CUTtemp)],each=N_hist)+
            rep(diff(CUTtemp)/N_hist, each=N_hist)*rep((1:N_hist)-1,N_xbin)
      time_bin <- makeCOVbin(plot_data$X,K=bintot.N,cutoffs=CUT[-1],
                           adjust0bin=FALSE,...)
   } else{
      cut_temp.id <- as.numeric(sort(unique(plot_data$X)))
      temp_diff <- diff(cut_temp.id)/2
      CUT <- c(cut_temp.id[1]-temp_diff[1],
               cut_temp.id[-length(cut_temp.id)] + temp_diff,
               cut_temp.id[length(cut_temp.id)] + temp_diff[length(temp_diff)])
      time_bin <- makeCOVbin(plot_data$X,K=bintot.N,
                           cutoffs=CUT,adjust0bin=FALSE,...)
   }

   caption <- paste(paste("No. of bins = ",N_xbin," / ",
                          "No. of Nhist = ", N_hist,"/",
                          "Percentile=",round(probs[1]*100),"th,",
                          round(probs[2]*100),"th,",
                          round(probs[3]*100),"th",sep=""),"/",
                          round(conf.level*100),"% CI")

   alpha <- 1-conf.level
   bintot.N <- nrow(time_bin$COVbin_summary)
  
   SQY <- data.frame(cut_temp = rep("NA",bintot.N*3),
                    quantile = paste0("Q",
                                  as.character(rep(probs,bintot.N)*100),"th"), 
                    V1 = rep(NA,bintot.N*3),
                    V2 = rep(NA,bintot.N*3),
                    V3 = rep(NA,bintot.N*3))                    
   FQY <- matrix(NA,nrow=bintot.N,ncol=3)

   for(i in 1:bintot.N){
      if(timebin_flag){
        if(i<N_hist){
           sel.id <- which(as.numeric(time_bin$COV_bin)<=i+N_hist-1)
        } else if(i>(bintot.N-N_hist+1)){
           sel.id <- which(as.numeric(time_bin$COV_bin)>=i-(N_hist-1))
        } else{
           sel.id <- which(as.numeric(time_bin$COV_bin)>i-N_hist &
                           as.numeric(time_bin$COV_bin)<i+N_hist)
        }
      } else{
         sel.id <- which(as.numeric(time_bin$COV_bin)==i)
      }
 
      low.point <- time_bin$COVbin_summary$lower_COV[i]
      upper.point <- time_bin$COVbin_summary$upper_COV[i]
      mid.point <- time_bin$COVbin_summary$mid_LU[i]      
      A <- as.numeric(time_bin$COV_bin[sel.id])
      if(weight_method=="bin"){
         temp <- N_hist-abs(A-i)
         temp.weight <- temp/N_hist
      } else{
         dist.temp <- abs(orig_data$TIME[sel.id]-mid.point)
         if(diff(range(dist.temp))!=0){
            temp.weight <- (max(dist.temp)-dist.temp)/diff(range(dist.temp))
         } else{
            temp.weight <-rep(1,length(dist.temp))
         }
      }

      suppressWarnings(
         temp.quantile <- t(apply(sim_data[sel.id,],2,function(x)
                                   wtd.quantile(x,weights=temp.weight,
                                   probs=probs,na.rm=TRUE))))
      suppressWarnings(
         temp.orig.q <- wtd.quantile(orig_data[,Y_name][sel.id],
                                weights=temp.weight,
                                probs=probs,na.rm=TRUE))

      temp <- t(apply(temp.quantile,2,
                      function(x) 
                       stats::quantile(x,prob=c(alpha/2,0.5,1-alpha/2),
                                         na.rm=TRUE)))

      SQY[(1:3)+(i-1)*3,3:5] <- temp
      SQY[(1:3)+(i-1)*3,1] <- as.character(time_bin$COVbin_summary$cut_temp[i])
      FQY[i,] <- c(temp.orig.q)
   }
   
   
   colnames(SQY)[3:5] <- paste0(as.character(c(alpha/2,0.5,1-alpha/2)*100),"%")
   
   FQY <- cbind(time_bin$COVbin_summary, FQY)
   nameT <- strsplit(rapply(list(probs),sprintf,
                              fmt="%.2f",how="replace")[[1]],"\\.")
   nameT <- unlist(lapply(nameT,function(x) x[2]))
   colnames(FQY)[7:9] <- paste0("Q",nameT,"th")

   Y_min <- as.numeric(min(c(unlist(FQY[,7]),
                           unlist(SQY[,3]),plot_data$Y),na.rm=T))
   Y_max <- as.numeric(max(c(unlist(FQY[,9]),
                           unlist(SQY[,5]),plot_data$Y),na.rm=T))

   Ptemp <- ggplot() + 
              ylim(Y_min,Y_max) +
              theme_bw() +
              theme(panel.grid.minor = element_blank(),
                    panel.grid.major=element_blank())

   if(bin_grid){
      temp_tick <- unique(unlist(time_bin$COVbin_summary[,c("lower_COV",
                                                            "upper_COV")]))
      mid_tick <- unlist(time_bin$COVbin_summary$mid_LU)
      Ptemp <- Ptemp +
         geom_vline(xintercept=temp_tick[2:(length(temp_tick)-1)],
                    linetype=3, color="grey70") +
         geom_rug(aes(mid_tick),color="grey70")
   }
   SQY$ID <- 1:nrow(SQY)
   TT <- sort(unique(SQY$quantile))
   D1 <- merge(SQY[SQY$quantile==TT[1],],FQY,by="cut_temp",
               all.x=TRUE,sort=FALSE) 
   D1 <- D1[order(D1$ID),c(1,3:5,11,12,9:10)]; 
   colnames(D1)[c(2:4,6)]<- c("LB","Mid","UB","Quant")
   D2 <- merge(SQY[SQY$quantile==TT[2],],FQY,by="cut_temp",
               all.x=TRUE,sort=FALSE) 
   D2 <- D2[order(D2$ID),c(1,3:5,11,13,9:10)]; 
   colnames(D2)[c(2:4,6)]<- c("LB","Mid","UB","Quant")
   D3 <- merge(SQY[SQY$quantile==TT[3],],FQY,by="cut_temp",
               ll.x=TRUE,sort=FALSE) 
   D3 <- D3[order(D3$ID),c(1,3:5,11,14,9:10)]; 
   colnames(D3)[c(2:4,6)]<- c("LB","Mid","UB","Quant")
      
   if(SIM_qCI){
 
      if(D1$mid_LU[1]>min(plot_data$X)){
        temp <- D1[1,]; temp$mid_LU <- min(plot_data$X)
        D1 <- rbind(temp,D1)
        temp <- D2[1,]; temp$mid_LU <- min(plot_data$X)
        D2 <- rbind(temp,D2)
        temp <- D3[1,]; temp$mid_LU <- min(plot_data$X)
        D3 <- rbind(temp,D3)  
      }
      if(D1$mid_LU[nrow(D1)]<max(plot_data$X)){
        temp <- D1[nrow(D1),]; temp$mid_LU <- max(plot_data$X)
        D1 <- rbind(D1,temp)
        temp <- D2[nrow(D2),]; temp$mid_LU <- max(plot_data$X)
        D2 <- rbind(D2,temp)
        temp <- D3[nrow(D3),]; temp$mid_LU <- max(plot_data$X)
        D3 <- rbind(D3,temp)  
      }     

      if(CIvpc_type=="line"){
         Ptemp <- Ptemp +
            geom_ribbon(data=D1,aes(x=mid_LU,ymin=LB,ymax = UB),
                        fill="lightblue",alpha=0.5)+
            geom_ribbon(data=D3,aes(x=mid_LU,ymin=LB,ymax = UB),
                        fill="lightblue",alpha=0.5)+
            geom_ribbon(data=D2,aes(x=mid_LU,ymin=LB,ymax = UB),
                        fill="lightpink",alpha=0.5)
      } else if(CIvpc_type=="segment"){
         D11 <- data.frame(LB = rep(D1$LB,each=2),UB = rep(D1$UB,each=2),
                          Quant = rep(D1$Quant,each=2),X = c(t(D1[,7:8])))  
         D22 <- data.frame(LB = rep(D2$LB,each=2),UB = rep(D2$UB,each=2),
                          Quant = rep(D2$Quant,each=2),X = c(t(D2[,7:8])))  
         D33 <- data.frame(LB = rep(D3$LB,each=2),UB = rep(D3$UB,each=2),
                          Quant = rep(D3$Quant,each=2),X = c(t(D3[,7:8])))  
         Ptemp <- Ptemp +
            geom_ribbon(data = D11,aes(x=X,ymin=LB,ymax=UB),
                        fill="lightblue",alpha=0.5)+
            geom_ribbon(data = D33,aes(x=X,ymin=LB,ymax=UB),
                        fill="lightblue",alpha=0.5)+
            geom_ribbon(data = D22,aes(x=X,ymin=LB,ymax=UB),
                        fill="lightpink",alpha=0.5)
      }
   }

   if(DV_point){
      Ptemp <- Ptemp+geom_point(data = plot_data,aes(X,Y),
                              color="grey30",size=pointsize,alpha=0.5)
   }

   if(DV_qline){
      if(CIvpc_type=="line"){
         Ptemp <- Ptemp +
               geom_line(data = D1,aes(x=mid_LU,y=Quant),linetype="dashed",
                         color="red",size=linesize)+
               geom_line(data = D2,aes(x=mid_LU,y=Quant),
                         color="red",size=linesize)+
               geom_line(data = D3,aes(x=mid_LU,y=Quant),linetype="dashed",
                         color="red",size=linesize)
      } else if(CIvpc_type=="segment"){
         Ptemp <- Ptemp +
               geom_line(data = D11,aes(x=X,y=Quant),linetype="dashed",
                      color="red",size=linesize)+
               geom_line(data = D22,aes(x=X,y=Quant),
                         color="red",size=linesize)+
               geom_line(data = D33,aes(x=X,y=Quant),linetype="dashed",
                      color="red",size=linesize)  
      }
   }

   if(SIM_qline){
      Ptemp <- Ptemp +
          geom_line(data = D1,aes(x=mid_LU,y=Mid),linetype="dashed",
                    color="blue",size=linesize)+
          geom_line(data = D2,aes(x=mid_LU,y=Mid),
                    color="blue",size=linesize)+
          geom_line(data = D3,aes(x=mid_LU,y=Mid),linetype="dashed",
                    color="blue",size=linesize) 
   }

   if(plot_caption){
      Ptemp <- Ptemp +
                   labs(caption=caption) +
                   theme(plot.caption = element_text(size=captionsize))
   }
   if(plot_flag){
      Ptemp + theme(panel.grid.major=element_line(colour = "white"),
                    panel.grid.minor=element_line(colour="white"),
                    plot.caption=element_text(hjust=0.5))+
                    labs(x=X_name,y=Y_name,title=main_title)
   } else{
     
      QTemp <- FQY[,c(1,6,7:9)] 
      colnames(QTemp)[1:2] <- c("X_bin","X_mid")

      STemp <- merge(SQY,FQY[,c(1,6)],by="cut_temp",all.x=TRUE,sort=FALSE)
      STemp <- STemp[,c(1,ncol(STemp),2:(ncol(STemp)-1))]
      colnames(STemp)[1:2] <- c("X_bin","X_mid")    

      DV_point <- orig_data[,c(X_name,Y_name)]
      return(list(DV_point=DV_point,
                  DV_quant=QTemp,
                  SIM_quantCI=STemp))      
  }
}





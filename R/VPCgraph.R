#' @title the original visual predictive check plot (VPC)
#' @description This function draws the original visual predictive check plot 
#' proposed by Holford & Karlsson (2008). The visual predictive check plot is 
#' a graphical comparison of observations and simulated predictions. It 
#' compares the distribution of quantiles obtained from the observed data 
#' to the distribution of quantiles from simulated data drawn 
#' from the fitted model. 
#' @usage VPCgraph(orig_data,
#'          sim_data,
#'          type = "CI",                   
#'          N_xbin = NULL,
#'          probs = c(0.1,0.5,0.9),
#'          conf.level = 0.95,
#'          X_name = "TIME",
#'          Y_name = "DV",
#'          MissingDV = NULL,
#'          DV_point = TRUE,
#'          CIvpc_type = "line",
#'          bin_grid = TRUE,
#'          plot_caption = TRUE,
#'          plot_flag = TRUE,
#'          linesize = 0.7,
#'          pointsize = 0.7,
#'          captionsize = 10,
#'          Kmethod = "cluster",                   
#'          maxK = NULL,
#'          beta = 0.2,
#'          lambda = 0.3,
#'          R = 4,
#'          C1 = 2.5,
#'          C2 = 7.8, ...)
#' @param orig_data A data frame of original data with X and Y variable.
#' @param sim_data A matrix of simulated data with only Y values collected.
#' @param type Type of VPC graph; "CI", "percentile", or "scatter".
#' @param N_xbin Number of bins in X variable. If NULL, optimal number of bins are automatically calcuated using optK function.
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
#' @param captionsize Size of caption .
#' @param maxK The maximum number of bins.
#' @param Kmethod The way to calculate the penalty in automatic binning."cluster" or "kernel". 
#' @param beta Additional parameter for automatic binning, used in optK function.  
#' @param lambda Additional parameter for automatic binning, used in optK function.  
#' @param R Additional parameter for automatic binning, used in optK function. 
#' @param C1 Additional parameter for automatic binning, used in optK function.  
#' @param C2 Additional parameter for automatic binning, used in optK function.   
#' @param ... Arguments to be passed to methods.
#' @return Visual predictive check plot or the values to draw VPC plot.
#' @references Holford N, & Karlsson M. (2008). "A tutorial on visual predictive checks,
#'  abstr 1434." Annual Meeting of the Populations Approach Group in Europe. www.page-meeting.org. 2008.
#' @references Harling, Uekcert, K. 2018. VPC and NPC User Guide. ICON plc. 
#' @references https://github.com/UUPharmacometrics/PsN/releases/download/4.9.0/vpc_npc_userguide.pdf.
#' @export
#' @examples
#' \donttest{
#' data(origdata)
#' data(simdata)
#' VPCgraph(origdata,simdata,type="CI",X_name="TIME",Y_name="DV",N_xbin=8)
#' }

# VPC graph ----------------------------------------------------------------

VPCgraph <- function(orig_data,
                     sim_data,
                     type = "CI",                   
                     N_xbin = NULL,
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
   
   
   main_title <- "Visual Predictive Check"   
   lower <- upper <- NA
   quantile<- mid_COV <- LB <- UB <- X <- Y <- Quant <- Mid <- NULL
   cut_temp <- mid_LU <- X_bin <- X_mid<-NULL
   sel.id <- !is.na(orig_data[,Y_name]) 
   if(!is.null(MissingDV))
      sel.id <- sel.id & orig_data[,MissingDV]==0
   sim_data <- sim_data[sel.id,]  
   orig_data <- orig_data[sel.id,]
   plot_data <- data.frame(X=orig_data[,X_name],Y=orig_data[,Y_name])
   
   if(is.null(N_xbin)){
      optK_t <- optK(plot_data$X,...)
      CUT <- optK_t$cutoffs
      time_bin <- optK_t$time_bin     
      N_xbin <- optK_t$K
   } else{
     if(N_xbin < length(table(plot_data$X))){
        CUT <- FindBestCut(plot_data$X,K=N_xbin,...)$cutoffs
        time_bin <- makeCOVbin(plot_data$X,K=N_xbin,cutoffs=CUT,...)
     } else{
        cut_temp.id <- as.numeric(sort(unique(plot_data$X)))
        temp_diff <- diff(cut_temp.id)/2
        nn <- length(cut_temp.id)
        CUT <- cut_temp.id[-nn]+temp_diff
        time_bin <- makeCOVbin(plot_data$X,K=N_xbin,cutoffs=CUT,...)
     }
   }   
   
   approxi <- ifelse(N_xbin >= length(table(plot_data$X)),TRUE,FALSE)   
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
   
   caption <- paste(paste("No. of bins = ",N_xbin," / ",
                          "Percentile=",round(probs[1]*100),"th,",
                          round(probs[2]*100),"th,",
                          round(probs[3]*100),"th",sep=""),"/",
                    round(conf.level*100),"% PI")
   
   FQY <- findQuantile(plot_data$Y,plot_data$X,
                     time_bin,probs=probs)
   
   SQY <- findSIMQuantile(sim_data,plot_data$X,
                        time_bin,probs=probs,conf.level,approxi)  
   
   Y_min <- as.numeric(min(c(unlist(FQY[,7]),
                           unlist(SQY[,3]),plot_data$Y),na.rm=T))
   Y_max <- as.numeric(max(c(unlist(FQY[,9]),
                           unlist(SQY[,5]),plot_data$Y),na.rm=T))
   
   Ptemp <- ggplot()+
      theme_bw()+
      theme(panel.grid.minor = element_blank(),
            panel.grid.major=element_blank())+ylim(Y_min,Y_max)
   
   if(bin_grid){
      temp_tick <- unique(unlist(time_bin$COVbin_summary[,c("lower_COV",
                                                          "upper_COV")]))
      mid_tick <- unlist(time_bin$COVbin_summary$mid_LU)
      Ptemp <- Ptemp +
         geom_vline(xintercept=temp_tick[2:(length(temp_tick)-1)],
                    linetype=3, color="grey70")+
         geom_rug(aes(mid_tick),color="grey70")
   }
   
   TT <- sort(unique(SQY$quantile))
   D1 <- merge(SQY[SQY$quantile==TT[1],],FQY,
               by = "cut_temp",all.x=TRUE,sort=FALSE)
   D1 <- D1[,c(1,3:5,10,11,8:9)];    
   colnames(D1)[c(2:4,6)] <- c("LB","Mid","UB","Quant")
   D2 <- merge(SQY[SQY$quantile==TT[2],],FQY,
               by = "cut_temp",all.x=TRUE,sort=FALSE)
   D2 <- D2[,c(1,3:5,10,12,8:9)];  
   colnames(D2)[c(2:4,6)] <- c("LB","Mid","UB","Quant")
   D3 <- merge(SQY[SQY$quantile==TT[3],],FQY,
               by = "cut_temp",all.x=TRUE,sort=FALSE)
   D3 <- D3[,c(1,3:5,10,13,8:9)];  
   colnames(D3)[c(2:4,6)] <- c("LB","Mid","UB","Quant")
   
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
                        fill="lightblue",alpha=0.5) +
            geom_ribbon(data=D3,aes(x=mid_LU,ymin=LB,ymax = UB),
                        fill="lightblue",alpha=0.5) +
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
                        fill="lightblue",alpha=0.5) +
            geom_ribbon(data = D33,aes(x=X,ymin=LB,ymax=UB),
                        fill="lightblue",alpha=0.5) +
            geom_ribbon(data = D22,aes(x=X,ymin=LB,ymax=UB),
                        fill="lightpink",alpha=0.5)
      }
   }
   if(DV_point){
      Ptemp <- Ptemp + geom_point(data = plot_data,aes(X,Y),
                              color="grey30",size=pointsize,alpha=0.5)
   }
   
   if(DV_qline){
      if(CIvpc_type=="line"){
         Ptemp <- Ptemp +
            geom_line(data = D1,aes(x=mid_LU,y=Quant),linetype="dashed",
                      color="red",size=linesize) +
            geom_line(data = D2,aes(x=mid_LU,y=Quant),
                      color="red",size=linesize) +
            geom_line(data = D3,aes(x=mid_LU,y=Quant),linetype="dashed",
                      color="red",size=linesize)
      } else if(CIvpc_type=="segment"){
         Ptemp <- Ptemp +
            geom_line(data = D11,aes(x=X,y=Quant),linetype="dashed",
                      color="red",size=linesize) +
            geom_line(data = D22,aes(x=X,y=Quant),
                      color="red",size=linesize) +
            geom_line(data = D33,aes(x=X,y=Quant),linetype="dashed",
                      color="red",size=linesize)  
      }
   }
   
   if(SIM_qline){
      Ptemp <- Ptemp +
         geom_line(data = D1,aes(x=mid_LU,y=Mid),linetype="dashed",
                   color="blue",size=linesize) +
         geom_line(data = D2,aes(x=mid_LU,y=Mid),
                   color="blue",size=linesize) +
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
      QTemp <- data.frame(X_bin=FQY$cut_temp,
                          X_mid = FQY$mid_LU,
                          FQY[,7:9])

      STemp <- merge(SQY,FQY[,c(1,6)],by = "cut_temp",all.x = TRUE,sort=FALSE)
      STemp <- STemp[,c(1,6,2:5)]
      colnames(STemp)[1:2] <- c("X_bin","X_mid")
      
      DV_point <- orig_data[,c(X_name,Y_name)]
      return(list(DV_point=DV_point,
                  DV_quant=QTemp,
                  SIM_quantCI=STemp))
   }
}








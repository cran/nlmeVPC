#' @title The visual predictive checks using the additive quantile regression (aqrVPC)
#' @description  This function draws the visual predictive check (VPC) plot 
#' using additive quantile regression. The quantile 
#' regression methods are used to calculate quantiles.
#' @usage aqrVPC(orig_data, 
#'        sim_data, 
#'        probs = c(0.1,0.5,0.9),
#'        conf.level = 0.95, 
#'        X_name = "TIME", 
#'        Y_name = "DV",
#'        MissingDV = NULL, 
#'        plot_caption = TRUE,
#'        DV_point = TRUE,
#'        plot_flag = TRUE,
#'        linesize = 0.7, 
#'        pointsize = 0.7, 
#'        captionsize = 10,
#'        qss_lambda = NULL, ...)
#' @param orig_data A data frame of original data with X and Y variable.
#' @param sim_data A matrix of simulated data with only Y values collected.
#' @param probs A numeric vector of probabilities.
#' @param conf.level Confidence level of the interval.
#' @param X_name Name of X variable in orig_data (usually "TIME" in pharmacokinetic data).
#' @param Y_name Name of Y variable in orig_data (usually "DV" in pharmacokinetic data).
#' @param MissingDV Name of missing indicator variable in orig_data, which have value 1 if missing, value 0 otherwise. (usually "MDV" in pharmacokinetic data).
#' @param DV_point Draw point (X, Y) in the plot if TRUE; omit if FALSE.
#' @param plot_caption Put caption with additional information if TRUE; omit if FALSE.
#' @param plot_flag Draw plot if TRUE; generate data for drawing plot if FALSE.
#' @param linesize Size of line in the plot.
#' @param pointsize Size of point in the plot.
#' @param captionsize Size of caption. 
#' @param qss_lambda Smoothing parameter in quantreg::qss function. 
#' Larger lambda produces a smoother fit.
#' @param ... Arguments to be passed to methods.
#' @return aqrVPC plot or the values to draw aqrVPC plot.
#' @references Koenker, Roger, and Kevin F. Hallock. "Quantile regression." 
#'             Journal of economic perspectives 15.4 (2001): 143-156.
#' @references Jamsen, K. M., Patel, K., Nieforth, K., & Kirkpatrick, C. M. 
#'             (2018). A regression approach to visual predictive checks for
#'              population pharmacometric models. CPT: pharmacometrics & 
#'              systems pharmacology, 7(10), 678-686.
#' @import ggplot2 quantreg optimx
#' @export
#' @examples
#'  \donttest{
#' data(origdata)
#' data(simdata)
#' aqrVPC(origdata,simdata)
#' }

aqrVPC <- function(orig_data,
                   sim_data,
                   probs = c(0.1,0.5,0.9),
                   conf.level = 0.95,
                   X_name = "TIME",
                   Y_name = "DV",
                   MissingDV = NULL,
                   plot_caption = TRUE, 
                   DV_point = TRUE,
                   plot_flag = TRUE,
                   linesize = 0.7,
                   pointsize = 0.7,
                   captionsize = 10,
                   qss_lambda = NULL, ...){
  
   main_title <- "Additive Quantile Regression VPC"  
   G <- X <- Y <- NA
   lower <- upper <- mid <- NA
   quantile <- LB <- UB <- Quant <- NULL
   sel.id <- !is.na(orig_data[,Y_name]) 
   if(!is.null(MissingDV))
      sel.id <- sel.id & orig_data[,MissingDV]==0
   sim_data <- sim_data[sel.id,]  
   orig_data <- orig_data[sel.id,]
   Lp <- (1-conf.level)/2
   Mp <- 0.5
   Up <- 1 - Lp
   caption <- paste(paste("Percentile=",round(probs[1]*100),"th,",
                          round(probs[2]*100),"th,",
                          round(probs[3]*100),"th",sep=""),"/",
                          round(conf.level*100),"% CI")   
   orig.aqr <- function(DV.data,TIME.data,q,lambda){
      data <- data.frame(X=TIME.data,Y=DV.data)
      orig.rqss.mod <- vector(mode="list",length=length(q))
      for(i in 1:length(q)){
         orig.rqss.mod[[i]] <- quantreg::rqss(Y~qss(X,lambda=lambda[i]),
                                            data=data,tau=q[i])
      }
      df <- data.frame(X=data$X)      
      p <- sapply(1:length(q),function(i){
        quantreg::predict.rqss(orig.rqss.mod[[i]],newdata=df)})
      
      df <- data.frame(X=rep(data$X,3), 
                       G=rep(paste("Q",round(q*100),"th",sep = ""),
                             each=nrow(data)),
                       Y = c(p))
      return(df)
   }

   N_sim <- ncol(sim_data)

   plot_data <- data.frame(orig_data,
                           X=orig_data[,X_name],
                           Y=orig_data[,Y_name])

   model.optim <- function(data,tau,lambda){
      mod.rqss <- rqss(Y~qss(X,lambda=lambda),data=data,tau=tau)
      stats::AIC(mod.rqss)
   }
 
   if(is.null(qss_lambda)){
      opt.lambda <- NULL
      for (i in 1:length(probs)){
         lambda <- stats::optimize(model.optim,c(1,10),
                    data=plot_data,tau=probs[i])$minimum
         opt.lambda[i] <- lambda
      }
   } else{
      opt.lambda <- qss_lambda
   }
   names(opt.lambda) <- paste("Q",round(probs*100,2),sep="")

   DV_quant <- orig.aqr(plot_data$Y,plot_data$X,q=probs,lambda=opt.lambda)
   DV_quant$temp <- paste(DV_quant$G,DV_quant$X)
   DV_quant <- DV_quant[!duplicated(DV_quant$temp),]
   DV_quant <- DV_quant[order(DV_quant$G,DV_quant$X),c("X","Y","G")]
   ntemp <- nrow(DV_quant)/3
   DV_quant <- data.frame(X= DV_quant$X[1:ntemp],
                           matrix(DV_quant$Y,ncol=3))   
   colnames(DV_quant)[-1] <- paste("Q",round(probs*100,2),sep="")
   
   SIM_quant <- orig.aqr(c(unlist(sim_data)),
                        rep(plot_data$X,N_sim),
                        q=probs,lambda=opt.lambda)
   SIM_quant$temp <- paste(SIM_quant$G,SIM_quant$X)
   SIM_quant <- SIM_quant[!duplicated(SIM_quant$temp),]
   SIM_quant <- SIM_quant[order(SIM_quant$G,SIM_quant$X),c("X","Y","G")]
   
   
   sim.Q.temp <- apply(sim_data,2,function(x)
                       orig.aqr(x,plot_data$X,q=probs,lambda=opt.lambda)[,3])

   sim.Q.temp <- array(sim.Q.temp,dim=c(nrow(sim_data),length(probs),N_sim))

   N <- nrow(orig_data)
   P <- length(probs)

   SIM_quant <- data.frame(X=rep(plot_data$X,P),
                           quantile=paste("Q",round(rep(probs,
                                     each=N)*100),"th",sep=""),
                           V1 = rep(NA,N*P),
                           V2 = rep(NA,N*P),
                           V3 = rep(NA,N*P))
                           
   for(i in 1:P){
      temp <-   t(apply(sim.Q.temp[,i,],1,function(x)
                            stats::quantile(x,probs=c(Lp,Mp,Up),na.rm=TRUE)))
      SIM_quant[(1:N)+(i-1)*N,3:5] <- temp
   }
   colnames(SIM_quant)[-(1:2)] <- paste0(as.character(c(Lp,Mp,Up)*100),"%")
   options(warn=0)

   if(plot_flag){
      Y_min <- as.numeric(min(c(unlist(DV_quant[,2]),
                            unlist(SIM_quant[,3]),plot_data$Y),na.rm=T))
      Y_max <- as.numeric(max(c(unlist(DV_quant[,4]),
                            unlist(SIM_quant[,5]),plot_data$Y),na.rm=T))

      Ptemp <- ggplot()+ylim(Y_min,Y_max)+
                theme_bw()+
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major=element_blank())
  
      TT <- sort(unique(SIM_quant$quantile))
      SIM_quant$ID <- 1:nrow(SIM_quant)
      D1 <- merge(SIM_quant[SIM_quant$quantile==TT[1],],
                  DV_quant,by="X",all.x=TRUE,sort = FALSE)
      D1 <- D1[order(D1$ID),c(1,3:5,7)]; 
      colnames(D1)[c(2:5)] <- c("LB","Mid","UB","Quant")
      
      D2 <- merge(SIM_quant[SIM_quant$quantile==TT[2],],
                  DV_quant,by="X",all.x=TRUE,sort = FALSE)
      D2 <- D2[order(D2$ID),c(1,3:5,8)]; 
      colnames(D2)[c(2:5)] <- c("LB","Mid","UB","Quant")
      
      D3 <- merge(SIM_quant[SIM_quant$quantile==TT[3],],
                  DV_quant,by="X",all.x=TRUE,sort = FALSE)
      D3 <- D3[order(D3$ID),c(1,3:5,9)];
      colnames(D3)[c(2:5)] <- c("LB","Mid","UB","Quant")
 
      Ptemp <- Ptemp+
                geom_ribbon(data=D1,aes(x=X,ymin=LB,ymax = UB),
                            fill="lightblue",alpha=0.5)+
                geom_ribbon(data=D3,aes(x=X,ymin=LB,ymax = UB),
                            fill="lightblue",alpha=0.5)+
                geom_ribbon(data=D2,aes(x=X,ymin=LB,ymax = UB),
                            fill="lightpink",alpha=0.5)+  
                geom_line(data = D1,aes(x=X,y=Quant),linetype="dashed",
                          color="red",size=linesize)+
                geom_line(data = D2,aes(x=X,y=Quant),
                          color="red",size=linesize)+
                geom_line(data = D3,aes(x=X,y=Quant),linetype="dashed",
                          color="red",size=linesize)   
 
      if(DV_point){
         Ptemp <- Ptemp+geom_point(aes(X,Y),data = plot_data,
                                   color="grey30",size=pointsize,alpha=0.5)
      }
      if(plot_caption){
         Ptemp <- Ptemp +
                     labs(caption=caption) +
                     theme(plot.caption = element_text(size=captionsize))
      } 
      Ptemp <- Ptemp+
                theme(panel.grid.major=element_line(colour = "white"),
                      panel.grid.minor=element_line(colour="white"),
                      plot.caption=element_text(hjust=0.5))+
               labs(x=X_name,y=Y_name,title=main_title)
      Ptemp
   } else{
      DV_point <- orig_data[,c(X_name,Y_name)]
      return(list(DV_point=DV_point,
                  DV_quant=DV_quant,
                  SIM_quantCI=SIM_quant))
   }
}





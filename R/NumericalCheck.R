#' @title The numerical predictive checks
#' @description This function calculates the numerical predictive checks for 
#' each prediction level.
#' @importFrom magrittr %>%
#' @usage NumericalCheck(orig_data,
#'                sim_data,
#'                N_xbin=NULL,
#'                pred.level=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
#'                conf.level=0.95,
#'                X_name="TIME",
#'                Y_name="DV",
#'                MissingDV = NULL,
#'                Kmethod="cluster",                
#'                maxK=NULL,
#'                beta=0.2,
#'                lambda=0.3,
#'                R=4,
#'                C1=2.5,
#'                C2=7.8,...)
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
#' NumericalCheck(origdata,simdata)$NPC
#' }

NumericalCheck<-function(orig_data,
                         sim_data,
                         N_xbin=NULL,
                         pred.level=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                         conf.level=0.95,
                         X_name="TIME",
                         Y_name="DV",
                         MissingDV = NULL,
                         Kmethod="cluster",                
                         maxK=NULL,
                         beta=0.2,
                         lambda=0.3,
                         R=4,
                         C1=2.5,
                         C2=7.8,...){

   Y = LB = UB = belowPIflag = abovePIflag = NA
   n = belowE = aboveE = belowPI = abovePI = NA
   belowPI_L = belowPI_U=abovePI_L = abovePI_U = NA
   PI = belowPI_p = abovePI_p = cut_temp=NA
   quantile=value=NA
   sel.id = !is.na(orig_data[,Y_name]) 
   if(!is.null(MissingDV))
      sel.id = sel.id & orig_data[,MissingDV]==0
   sim_data = sim_data[sel.id,]  
   orig_data = orig_data[sel.id,]
   
   DV_quant<-NULL
   SIM_quant<-NULL
   temp_simCI <- NULL
   ID<-NULL;G<-NULL
   N_sim = ncol(sim_data)
   if(is.null(N_xbin))
      N_xbin<-optK(orig_data[,X_name],...)$K
   
   plot_data<-data.frame(X=orig_data[,X_name],Y=orig_data[,Y_name])

   if(N_xbin < length(table(plot_data$X))){
      CUT = FindBestCut(plot_data$X,K=N_xbin,...)$cutoffs
      time_bin<-makeCOVbin(plot_data$X,N_xbin,cutoffs=CUT,...)
   } else{
      cut_temp.id<-as.numeric(names(table(plot_data$X)))
      temp_diff<-diff(cut_temp.id)/2
      CUT<-c(cut_temp.id[1]-temp_diff[1],
             cut_temp.id[-length(cut_temp.id)]+temp_diff,
             cut_temp.id[length(cut_temp.id)]+
                temp_diff[length(temp_diff)])
      CUT =CUT[CUT<max(plot_data$X) & CUT>min(plot_data$X)]
      time_bin<-makeCOVbin(plot_data$X,N_xbin,cutoffs=CUT,...)
   }

   plot_data = data.frame(plot_data, cut_temp = time_bin$COV_bin)
   keepAll = NULL
   probkeep = NULL
   keepAll2 = NULL
   CIkeep = data.frame(X_bin=time_bin$COVbin_summary$cut_temp)
   for(i in 1:length(pred.level)){
      probs = c((1-pred.level[i])/2,1-(1-pred.level[i])/2)
      probkeep = c(probkeep,probs)
      probs = probs[!duplicated(probs)]
      temp_simQ<-findSIMQuantile(sim_data,plot_data$X,time_bin,probs=probs)
      if(length(probs)==1){
         temp_simQ = temp_simQ[,c(1,4,4)]
      } else{
         temp_simQ = temp_simQ[,c(1,2,4)]
         colnames(temp_simQ)[3]="value"
         temp_simQ=temp_simQ %>% spread(quantile,value)
      }
      colnames(temp_simQ)[2:3] = c("LB","UB")
      temp_data = plot_data %>% dplyr::left_join(temp_simQ,by="cut_temp")
      temp1 = temp_data %>%
         group_by(cut_temp) %>%
         dplyr::mutate(belowPIflag= Y<LB,
                       abovePIflag = Y>UB) %>%
         dplyr::summarize(belowPI = sum(belowPIflag),
                          abovePI = sum(abovePIflag),
                          n=n())%>%
         dplyr::mutate(pred.level=pred.level[i],
                       belowE = n*(1-pred.level)/2,
                       aboveE = n*(1-pred.level)/2)
      
      keepAll2 = rbind(keepAll2, temp1)
      temp1 = unlist(temp1 %>% summarize(n=sum(n),
                                  belowE=sum(belowE),
                                  aboveE=sum(aboveE),
                                  belowPI=sum(belowPI),
                                  abovePI=sum(abovePI)))
      tempAll = t(apply(sim_data,2,
                        function(x) return(c(sum(x<temp_data$LB),
                                             sum(x>temp_data$UB)))))
      temp2 = apply(tempAll,2,
                    function(x) stats::quantile(x,
                                    probs=c((1-conf.level)/2,
                                            1-(1-conf.level)/2)))
      names(temp2)= c("belowPI_L","belowPI_U","abovePI_L","abovePI_U")
      keepAll = rbind(keepAll,c(temp1,temp2))
      colnames(temp_simQ)[2:3] = paste0(as.character(pred.level[i]*100),
                                        "%",colnames(temp_simQ)[2:3])
      CIkeep = cbind(CIkeep,temp_simQ[,2:3])
   }

   NPC = keepAll %>% dplyr::as_tibble() %>%
      dplyr::mutate(PI = pred.level*100,
             belowPI_p=belowPI/n*100,
             abovePI_p=abovePI/n*100,
             belowPI_L=belowPI_L/n*100,
             belowPI_U=belowPI_U/n*100,
             abovePI_L=abovePI_L/n*100,
             abovePI_U=abovePI_U/n*100) %>%
      dplyr::select(PI,n,belowE,belowPI,belowPI_p,belowPI_L,belowPI_U,
             aboveE,abovePI,abovePI_p,abovePI_L,abovePI_U)


   colnames(NPC) = c("PI","n",
                     "Expected points below PI",
                     "points below PI","points below PI(%)",
                     paste0(round(conf.level*100),"%CIBelowFrom(%)"),
                     paste0(round(conf.level*100),"%CIBelowTo(%)"),
                     "Expected points above PI",
                     "points above PI","points above PI(%)",
                     paste0(round(conf.level*100),"%CIAboveFrom(%)"),
                     paste0(round(conf.level*100),"%CIAboveTo(%)"))

   return(list(NPC=NPC,
               DV_point=plot_data,
               sim_quant=CIkeep,
               NPCcut = keepAll2))
}








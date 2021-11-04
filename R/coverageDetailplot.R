#' @title draw the detailed coverage plot
#' @description This function draws the detailed coverage plot for 
#' the specific prediction level.
#' @import ggplot2 quantreg optimx
#' @usage coverageDetailplot(orig_data,
#'                   sim_data,
#'                   N_xbin=NULL,
#'                   predL=0.5,
#'                   conf.level=0.95,
#'                   X_name="TIME",
#'                   Y_name="DV",
#'                   MissingDV = NULL,
#'                   subject_name="ID",
#'                   Kmethod="cluster",                
#'                   maxK=NULL,
#'                   beta=0.2,
#'                   lambda=0.3,
#'                   R=4,
#'                   C1=2.5,
#'                   C2=7.8,...)
#' @title Visual predictive checks
#' @param orig_data A data frame of original data with X and Y variable.
#' @param sim_data A matrix of simulated data with only Y values collected.
#' @param N_xbin Number of bins in X variable. If NULL, optimal number of bins are automatically calcuated using optK function.
#' @param predL Scalar of probability
#' @param conf.level Confidence level of the interval.
#' @param X_name Name of X variable in orig_data (usually "TIME" in pharmacokinetic data).
#' @param Y_name Name of Y variable in orig_data (usually "DV" in pharmacokinetic data)
#' @param MissingDV Name of missing indicator variable in orig_data, which have value 1 if missing, value 0 otherwise. (usually "MDV" in pharmacokinetic data).
#' @param subject_name Name of subject variable in orig_data (usually "ID" in pharmacokinetic data).
#' @param maxK The maximum number of bins
#' @param Kmethod The way to calculate the penalty in automatic binning."cluster" or "kernel". 
#' @param beta Additional parameter for automatic binning, used in optK function.  
#' @param lambda Additional parameter for automatic binning, used in optK function.  
#' @param R Additional parameter for automatic binning, used in optK function. 
#' @param C1 Additional parameter for automatic binning, used in optK function.  
#' @param C2 Additional parameter for automatic binning, used in optK function. 
#' @param ... Arguments to be passed to methods.
#' @return the detailed coverage plot 
#' @references Post, T. M., et al. (2008)
#' Extensions to the visual predictive check for facilitate model performance
#' evaluation, Journal of pharmacokinetics and pharmacodynamics,
#' 35(2), 185-202
#' @export
#' @examples
#' \donttest{
#' data(origdata)
#' data(simdata)
#' coverageDetailplot(origdata,simdata,predL=0.5)
#' }

coverageDetailplot<-function(orig_data,
                   sim_data,
                   N_xbin=NULL,
                   predL=0.5,
                   conf.level=0.95,
                   X_name="TIME",
                   Y_name="DV",
                   MissingDV = NULL,
                   subject_name="ID",
                   Kmethod="cluster",                
                   maxK=NULL,
                   beta=0.2,
                   lambda=0.3,
                   R=4,
                   C1=2.5,
                   C2=7.8,...){
   
   belowPI= abovePI=Lower=Upper=cut_temp=Middle=type=percent=LE=UE=NA
   main_title=paste("Coverage Detailed plot : PI = ",predL*100)
   probs = c((1-predL)/2,1-(1-predL)/2)
   pred.level=predL
   NPCcut=NumericalCheck(orig_data,sim_data,pred.level = pred.level,
                         conf.level=conf.level,
                         X_name=X_name,Y_name=Y_name,...)$NPCcut
   NPCcut %>% 
      filter(pred.level==predL) %>% 
      mutate(Lower = belowPI/n,Upper = abovePI/n,Middle=1-Lower-Upper) %>% 
      select(cut_temp,Lower,Upper,Middle) %>%
      gather(key="type",value="percent",-cut_temp) %>%
      mutate(type=factor(type,levels=c("Upper","Middle","Lower")),
             LE = probs[1],UE=probs[2]) %>%
      ggplot() +
      geom_bar(aes(cut_temp,percent,fill=type),stat="identity",
               position="stack")+
      geom_point(aes(cut_temp,LE),col="white")+
      geom_point(aes(cut_temp,UE),col="white")+
      scale_y_continuous(expand=c(0,0))+
      scale_fill_manual("Type",values = c("grey45","grey75","grey20"))+
      theme_bw()+
      labs(x=X_name,y="Percentage Observations(%)",title=main_title)+
      scale_x_discrete(labels=levels(NPCcut$cut_temp))+
      theme(axis.text.x=element_text(angle=45,hjust=1))
   
 
}





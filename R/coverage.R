#' @title coverage plot
#' @description draw the coverage plot that visualize the result 
#' of the numerical predictive checks.
#' @import ggplot2 tidyr dplyr
#' @usage coverageplot(orig_data,
#'              sim_data,
#'              N_xbin=NULL,
#'              pred.level=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
#'              conf.level=0.95,                         
#'              X_name="TIME",
#'              Y_name="DV",
#'              MissingDV = NULL,
#'              plot_flag=TRUE,
#'              linesize=0.7,
#'              pointsize=1.5,
#'              Kmethod="cluster",                
#'              maxK=NULL,
#'              beta=0.2,
#'              lambda=0.3,
#'              R=4,
#'              C1=2.5,
#'              C2=7.8,...)
#' @param orig_data A data frame of original data with X and Y variable.
#' @param sim_data A matrix of simulated data with only Y values collected.
#' @param N_xbin Number of bins in X variable. If NULL, optimal number of bins are automatically calculated using optK function.
#' @param pred.level Numeric vector of probabilities.
#' @param conf.level Confidence level of the interval.
#' @param X_name Name of X variable in orig_data (usually "TIME" in pharmacokinetic data).
#' @param Y_name Name of Y variable in orig_data (usually "DV" in pharmacokinetic data).
#' @param MissingDV Name of missing indicator variable in orig_data, which have value 1 if missing, value 0 otherwise. (usually "MDV" in pharmacokinetic data).
#' @param plot_flag Draw plot if TRUE; generate data for drawing plot if FALSE.
#' @param linesize Size of line in the plot.
#' @param pointsize Size of point in the plot.
#' @param maxK Yhe maximum number of bins.
#' @param Kmethod The way to calculate the penalty in automatic binning."cluster" or "kernel". 
#' @param beta Additional parameter for automatic binning, used in optK function.  
#' @param lambda Additional parameter for automatic binning, used in optK function.  
#' @param R Additional parameter for automatic binning, used in optK function. 
#' @param C1 Additional parameter for automatic binning, used in optK function.  
#' @param C2 Additional parameter for automatic binning, used in optK function. 
#' @param ... arguments to be passed to methods
#' @return coverage plot
#' @references Holford N, & Karlsson M. (2008). "A tutorial on visual predictive checks,
#'  abstr 1434." Annual Meeting of the Populations Approach Group in Europe. www.page-meeting.org. 2008.
#' @export
#' @examples
#' \donttest{
#' data(origdata)
#' data(simdata)
#' coverageplot(origdata,simdata)
#' }
#'

coverageplot <- function(orig_data,
                         sim_data,
                         N_xbin=NULL,
                         pred.level=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                         conf.level=0.95,                         
                         X_name="TIME",
                         Y_name="DV",
                         MissingDV = NULL,
                         plot_flag=TRUE,
                         linesize=0.7,
                         pointsize=1.5,
                         Kmethod="cluster",                
                         maxK=NULL,
                         beta=0.2,
                         lambda=0.3,
                         R=4,
                         C1=2.5,
                         C2=7.8,...){
   
   sel.id = !is.na(orig_data[,Y_name]) 
   if(!is.null(MissingDV))
     sel.id = sel.id & orig_data[,MissingDV]==0
   sim_data = sim_data[sel.id,]  
   orig_data = orig_data[sel.id,]
  
   PI =  PBp=prob=PAp=BF=BT=AF=AT=NA
   L =U= Ratio=NA
   RatioLower=RatioUpper=G=LLower=LUpper=ULower=UUpper=NA

   if(is.null(N_xbin))
      N_xbin<-optK(orig_data[,X_name],...)$K
   NPC = NumericalCheck(orig_data,sim_data,N_xbin=N_xbin,
                        X_name=X_name,Y_name=Y_name,
                        pred.level=pred.level,conf.level=conf.level,...)$NPC
   colnames(NPC)[-c(1:2)]= c("EPB","PB","PBp","BF","BT",
                             "EPA","PA","PAp","AF","AT")
   NPCtemp = NPC %>%
     dplyr::transmute(PI = PI,
               prob = (100-PI)/2,
               RatioLower=PBp/prob,
               RatioUpper=PAp/prob,
               LLower=BF/prob,
               ULower=BT/prob,
               LUpper=AF/prob,
               UUpper=AT/prob)
  
   T1 = NPCtemp %>% 
      dplyr::select(PI,RatioLower,RatioUpper) %>%
      tidyr::gather(key="G",value="Ratio",-PI) %>% 
      tidyr::separate(G,c("t","G"),sep=5) %>%
      dplyr::select(-t)
   T2 = NPCtemp %>% 
      dplyr::select(PI,LLower,LUpper) %>%
      tidyr::gather(key="G",value="L",-PI) %>% 
      tidyr::separate(G,c("t","G"),sep=1) %>%
      dplyr::select(-t)
   T3 = NPCtemp %>% 
      dplyr::select(PI,ULower,UUpper) %>%
      tidyr::gather(key="G",value="U",-PI) %>% 
      tidyr::separate(G,c("t","G"),sep=1) %>%
      dplyr::select(-t)
   PI.temp= T1 %>% 
      dplyr::left_join(T2, by = c("PI", "G")) %>% 
      dplyr::left_join(T3, by = c("PI", "G"))
   if(plot_flag){
     Y_min <- 0
     Y_max <- max(PI.temp$U)+0.2

     PI.temp %>%
        ggplot(aes(x=PI,y=Ratio)) +
        geom_ribbon(aes(ymin=L,ymax=U),alpha=0.3) +
        geom_hline(yintercept=1,color="white",size=linesize*1.3)+
        geom_point(size=pointsize)+geom_line(size=linesize) +
        ylim(c(Y_min,Y_max)) + 
        facet_wrap(~G,ncol=1)+theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major=element_blank())+
        labs(y="Obs/Exp",title="Coverage Plot")
   } else{
      NPCtemp = NPCtemp %>% dplyr::select(PI,RatioLower,LLower,ULower,
                                    RatioUpper,LUpper,UUpper)
      colnames(NPCtemp) = c("PI","Lower",
                            paste0("Lower",round(conf.level*100),"%CI_LB"),
                            paste0("Lower",round(conf.level*100),"%CI_UB"),
                            "Upper",
                            paste0("Upper",round(conf.level*100),"%CI_LB"),
                            paste0("Upper",round(conf.level*100),"%CI_UB"))                           
     return(list(NPC=as_tibble(NPC),coverage=as_tibble(NPCtemp)))
   }
}

















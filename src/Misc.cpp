#include "RcppArmadillo.h"
#include "Rcpp.h"

using namespace Rcpp; 
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]



//' @title Find quantiles of the simulated data using Rcpp
//' @usage findSIMQ(SIM,
//'         X,
//'         Xbin,
//'         probs,
//'         confLevel,
//'         approx)
//' @param SIM A matrix of simulated data with only Y values collected.
//' @param X A numeric vector corresponding to Y
//' @param Xbin Binning result from makeCOVbin function
//' @param probs A numeric vector of probabilities
//' @param confLevel Confidence level of the interval.
//' @param approx Arguments to be passed to methods
//' @return quantiles of SIM using xbin
//' @export
//' @examples
//' data(origdata)
//' data(simdata)
//' CUT = FindBestCut(origdata$TIME,8)$cutoffs
//' time_bin = makeCOVbin(origdata$TIME,K=8,cutoffs = CUT)
//' findSIMQ(simdata,origdata$TIME,Xbin=time_bin,probs=c(0.1,0.5,0.9),
//' confLevel=0.95,approx=FALSE)
//' @exportPattern "^[[:alpha:]]+"
//' @importFrom Rcpp evalCpp
//' @useDynLib nlmeVPC
// [[Rcpp::export]]
extern "C" SEXP findSIMQ(SEXP SIM,SEXP X, SEXP Xbin, SEXP probs, 
                               SEXP confLevel, SEXP approx){
      Rcpp::NumericMatrix SIM_(SIM);
      Rcpp::NumericVector X_ (X);
      Rcpp::List Xbin_ (Xbin);
      Rcpp::IntegerVector COV_bin = Xbin_["COV_bin"];
      Rcpp::List COVbin_summary = Xbin_["COVbin_summary"];
      Rcpp::NumericVector probs_(probs);
      double confLevel_ =Rcpp::as<double>(confLevel);
      bool approx_ = Rcpp::as<bool>(approx);
    
      Environment pkg1=Environment::namespace_env("stats");
      Function STATQUANT = pkg1["quantile"];     
      Function QNORM = pkg1["qnorm"];  
      Rcpp::IntegerVector CUTTEMP=COVbin_summary["cut_temp"];       
      int Nsim = SIM_.ncol();
      int n = SIM_.nrow();
    
      Rcpp::IntegerVector nbinN = COVbin_summary["n_bin"];
      int Nbin= nbinN.length();
      int lq  = probs_.length();        
      Rcpp::NumericMatrix SIM_summary(lq*Nbin,3);  
      Rcpp::IntegerVector cut_temp(lq*Nbin);
      Rcpp::NumericVector quantile(lq*Nbin);   

        
      int m = max(nbinN);       
      if(approx_){
        Rcpp::NumericMatrix Vkeep(m*Nsim,Nbin);   
        std::fill( Vkeep.begin(), Vkeep.end(), NumericVector::get_na() ) ;  
 
        for(int i=0; i<Nbin; i++){ 
          int count=0;
          for(int j=0; j<n; j++){
            for(int k=0; k<Nsim; k++){
              if(COV_bin(j)== (i+1)){
                Vkeep(count,i) = SIM_(j,k);
                count++;
              }
            }
          }
        }
        Rcpp::NumericMatrix confP(lq,3);
        Rcpp::NumericVector tt(1); tt(0)=1-(1-confLevel_)/2;
        Rcpp::NumericVector Qtt = QNORM(tt);
        
        confP.column(0)=probs_-Qtt(0)*sqrt(probs_*(1-probs_)/n);
        confP.column(1)=probs_;
        confP.column(2)=probs_+Qtt(0)*sqrt(probs_*(1-probs_)/n);        

        for(int i=0; i<Nbin; i++){
          for(int j=0; j<lq; j++){
             Rcpp::NumericVector AA = STATQUANT(Vkeep.column(i),
                       Rcpp::Named("probs")=confP.row(j),
                        Rcpp::Named("na.rm")=1);
             SIM_summary.row(i*lq+j) = AA;
             cut_temp(i*lq+j)=CUTTEMP(i);
             quantile(i*lq+j)=probs_(j);
          }
        }   
      } else{
       arma::cube SIM_Q(Nsim,lq,Nbin) ;
       for(int II=0; II<Nsim; II++){
          Rcpp::IntegerVector countbin(Nbin);
          Rcpp::NumericVector simtemp = SIM_.column(II);
          Rcpp::NumericMatrix transM(m,Nbin);
          std::fill( transM.begin(), transM.end(), NumericVector::get_na() ) ;
          
          for(int i=0; i<n; i++){
             transM(countbin(COV_bin(i)-1),COV_bin(i)-1)=simtemp(i) ;  
             countbin(COV_bin(i)-1) ++;
          }
          for(int i=0; i<Nbin; i++){
           
            Rcpp::NumericVector AA = STATQUANT(transM.column(i),
                        Rcpp::Named("probs")=probs_,
                        Rcpp::Named("na.rm")=1);
           
            for(int j=0; j<lq; j++)
               SIM_Q(II,j,i) = AA(j);
             
          }
             
       }
      
       Rcpp::NumericVector confP(3); 
       confP(0)=(1-confLevel_)/2.0;
       confP(1)=0.5;
       confP(2)=1-(1-confLevel_)/2.0;

       for(int i=0; i<Nbin; i++){
          for(int j=0; j<lq; j++){
             Rcpp::NumericVector AA = STATQUANT(SIM_Q.slice(i).col(j),
                       Rcpp::Named("probs")=confP,
                        Rcpp::Named("na.rm")=1);
             SIM_summary.row(i*lq+j) = AA;
             cut_temp(i*lq+j)=CUTTEMP(i);
             quantile(i*lq+j)=probs_(j);
          }
       }

    }
    
    return Rcpp::List::create(
    Rcpp::Named("cuttemp") = cut_temp,
    Rcpp::Named("quantile") = quantile,
    Rcpp::Named("quantileV") =SIM_summary );
   
} 



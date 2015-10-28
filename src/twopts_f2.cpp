#include <Rcpp.h>
#include <R_ext/PrtUtil.h>
#include "f2_est.h"
#include "utils.h"
using namespace Rcpp;
using namespace std;

RcppExport SEXP est_rf_f2(SEXP geno_R, 
			  SEXP segreg_type_R, 
			  SEXP n_ind_R, 
			  SEXP bar_width_R)
{
  int n_ind = Rcpp::as<int>(n_ind_R);
  int bar_width = Rcpp::as<int>(bar_width_R);
  Rcpp::NumericVector segreg_type = Rcpp::as<Rcpp::NumericVector>(segreg_type_R);
  Rcpp::NumericVector geno = Rcpp::as<Rcpp::NumericVector>(geno_R);\
  int n_mar=((int)geno.size()/n_ind), k1, k2;
  NumericMatrix r(n_mar, n_mar);
  Rcpp::NumericVector rtemp(2);   
  for(int i=0; i < n_mar-1; i++)
    {
      R_CheckUserInterrupt(); // check for ^C 
      for(int j=(i+1); j  < n_mar; j++)
        {
	  std::vector<int> k_sub(&geno[i*n_ind],&geno[i*n_ind+n_ind]);
	  std::vector<int> k1_sub(&geno[j*n_ind],&geno[j*n_ind+n_ind]);
	  k1=segreg_type(i); k2=segreg_type(j);
	  // Rcpp::Rcout << k1 << "--" << k2 << "\n";
	  switch(k1){
	  case 1: //C
	    switch(k2){
	    case 1: 
	      rtemp = est_rf_C_C(k_sub, k1_sub, n_ind);  //Markers: Codominant - Codominant
	      break;
	    case 2:
	      rtemp = est_rf_C_D_43(k_sub, k1_sub, n_ind);  //Markers: Codominant - Dominant( (12)=4, 3 )
	      break;
	    case 3:
	      rtemp = est_rf_C_D_51(k_sub, k1_sub, n_ind);  //Markers Codominant - Dominant( (1, (23)=5) 
	      break;
	    case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind);  //Any type of markers - slower function
	      break;
	    }
	    break;
	  case 2: //D_43
	    switch(k2){
	    case 1: 
	      rtemp = est_rf_C_D_43(k1_sub, k_sub, n_ind);  //Markers: Dominant( (12)=4, 3 ) - Codominant
	      break;  
	    case 2:
	      rtemp = est_rf_D_D_43(k_sub, k1_sub, n_ind);  //Markers: Dominant( (12)=4, 3 ) - Dominant( (12)=4, 3 )
	      break;
	    case 3:
	      rtemp = est_rf_D_D_43_51(k_sub, k1_sub, n_ind);  //Markers: Dominant( (12)=4, 3 ) - Dominant( (1, (23)=5) 
	      break;
	    case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind);  //Any type of markers - slower function
	      break;
	    }
	    break;
	  case 3: //D51
	    switch(k2){
	    case 1: 
	      rtemp = est_rf_C_D_51(k1_sub, k_sub, n_ind); //Markers  Dominant( (1, (23)=5) - Codominant
	      break;
	    case 2:
	      rtemp = est_rf_D_D_43_51(k1_sub, k_sub, n_ind); //Markers:  Dominant( (1, (23)=5) - Dominant( (12)=4, 3 )
	      break;
	    case 3:
	      rtemp = est_rf_D_D_51(k_sub, k1_sub, n_ind); //Markers:  Dominant( (1, (23)=5) - Dominant( (1, (23)=5)
	      break;
	    case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind);  //Any type of markers - slower function
	      break;
	    }
	    break;
	  case 4:
	    switch(k2){
	    case 1: case 2: case 3: case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind); //Any type of markers - slower function
	      break;
	    }
	    break;
	  }
	  r(j,i)=rtemp(0);
	  r(i,j)=rtemp(1);
	  rtemp(0)=rtemp(1)=0;
	}
    }
  return wrap(r);
}

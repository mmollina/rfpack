#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.00001
#define rf_TOL 1e-50
#define LN_5 log(1.0/2.0)

RcppExport SEXP est_rf_bc(SEXP geno_R, 
			  SEXP n_ind_R, 
			  SEXP bar_width_R)
{
  Rcpp::NumericVector geno = Rcpp::as<Rcpp::NumericVector>(geno_R);
  int n_ind = Rcpp::as<int>(n_ind_R);
  int bar_width = Rcpp::as<int>(bar_width_R);
  int n_mar=((int)geno.size()/n_ind);
  double rtemp, l, l0, mis=0, nr=0;
  Rcpp::NumericMatrix r(n_mar, n_mar);
  for(int i=0; i < (int)(geno.size()/n_ind)-1; i++)
    {
      R_CheckUserInterrupt(); /* check for ^C */
      for(int j=(i+1); j  < (int)geno.size()/n_ind; j++)
	{
	  nr=mis=0;
	  std::vector<int> k_sub(&geno[i*n_ind],&geno[i*n_ind+n_ind]);
	  std::vector<int> k1_sub(&geno[j*n_ind],&geno[j*n_ind+n_ind]);
	  for(int k=0; k < n_ind; k++)
	    {
	      if((k_sub[k] + k1_sub[k])==3) nr++;
	      else if(k_sub[k]==0 || k1_sub[k]==0) mis++;
	    }
	  rtemp=nr/(n_ind-mis);
	  if(rtemp < rf_TOL) rtemp=rf_TOL;
	  l=(n_ind-(mis+nr))*log(1-rtemp)+nr*log(rtemp);
	  l0=LN_5*(n_ind-mis);
	  r(j,i)=rtemp;
	  r(i,j)=(l-l0)/log(10.0);
	}
    }
  return wrap(r);
}

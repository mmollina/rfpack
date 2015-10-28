#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define TOL 0.00001
#define rf_TOL 1e-50
#define LN_5 log(1.0/2.0)

RcppExport SEXP est_rf_bc(SEXP geno_R, 
			  SEXP n_ind_R, 
			  SEXP bar_width_R);

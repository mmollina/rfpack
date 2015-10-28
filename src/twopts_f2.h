#include <Rcpp.h>
#include <R_ext/PrtUtil.h>
#include "f2_est.h"
#include "utils.h"
using namespace Rcpp;
using namespace std;

RcppExport SEXP est_rf_f2(SEXP geno_R, 
			  SEXP segreg_type_R, 
			  SEXP n_ind_R, 
			  SEXP bar_width_R);

#include <Rcpp.h>
#include <R_ext/PrtUtil.h>
#include "out_est.h"
#include "utils.h"
using namespace Rcpp;
using namespace std;

RcppExport SEXP est_rf_out(SEXP geno_R, 
			   SEXP segreg_type_R, 
			   SEXP n_ind_R, 
			   SEXP bar_width_R);

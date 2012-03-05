#ifndef _UNMARKED_NLL_DISTSAMP_H
#define _UNMARKED_NLL_DISTSAMP_H

#include "R_ext/Applic.h"
#include "RcppArmadillo.h"


RcppExport SEXP nll_distsamp( SEXP y_, SEXP X_, SEXP V_, SEXP beta_lam_, SEXP beta_sig_, SEXP X_offset_, SEXP V_offset_, SEXP A_, SEXP a_, SEXP u_, SEXP output_, SEXP db_ ) ;

#endif

#ifndef _unmarked_NLL_DISTSAMP_H
#define _unmarked_NLL_DISTSAMP_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_distsamp( SEXP y_, SEXP X_, SEXP V_, SEXP beta_lam_, SEXP beta_sig_, SEXP X_offset_, SEXP V_offset_, SEXP A_, SEXP a_, SEXP u_, SEXP output_, SEXP db_ ) ;

#endif

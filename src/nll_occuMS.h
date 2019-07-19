#ifndef _unmarked_NLL_OCCUMS_H
#define _unmarked_NLL_OCCUMS_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_occuMS( SEXP beta_, SEXP y_, 
    SEXP dm_state_, SEXP dm_phi_, SEXP dm_det_, 
    SEXP sind_, SEXP pind_, SEXP dind_, SEXP prm_, 
    SEXP S_, SEXP T_, SEXP J_, SEXP N_,
    SEXP naflag_, SEXP guide_) ;

#endif

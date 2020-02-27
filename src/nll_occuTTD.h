#ifndef _unmarked_NLL_OCCUTTD_H
#define _unmarked_NLL_OCCUTTD_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_occuTTD( SEXP beta_, SEXP y_, SEXP delta_,
    SEXP W_, SEXP V_, SEXP Xgam_, SEXP Xeps_, 
    SEXP pind_, SEXP dind_, SEXP cind_, SEXP eind_, 
    SEXP lpsi_, SEXP tdist_,
    SEXP N_, SEXP T_, SEXP J_,
    SEXP naflag_) ;

#endif

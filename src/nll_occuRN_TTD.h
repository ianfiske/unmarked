#ifndef _unmarked_NLL_OCCURN_TTD_H
#define _unmarked_NLL_OCCURN_TTD_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_occuRN_TTD( SEXP beta_, SEXP y_, SEXP delta_,
    SEXP W_, SEXP V_, SEXP pind_, SEXP dind_, SEXP tdist_,
    SEXP N_, SEXP J_, SEXP K_, SEXP naflag_) ;

#endif

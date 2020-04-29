#ifndef _unmarked_NLL_NMIXTTD_H
#define _unmarked_NLL_NMIXTTD_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_nmixTTD( SEXP beta_, SEXP y_, SEXP delta_,
    SEXP W_, SEXP V_, SEXP pind_, SEXP dind_, SEXP tdist_,
    SEXP N_, SEXP J_, SEXP K_, SEXP naflag_) ;

#endif

#ifndef _unmarked_NLL_MULTINOMPOIS_H
#define _unmarked_NLL_MULTINOMPOIS_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_multinomPois(SEXP betaR, SEXP pi_funR, 
    SEXP XlamR, SEXP Xlam_offsetR, SEXP XdetR, SEXP Xdet_offsetR,  
    SEXP yR, SEXP navecR, SEXP nPr, SEXP nAPr) ;

#endif

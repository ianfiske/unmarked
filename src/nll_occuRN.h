#ifndef _unmarked_NLL_OCCURN_H
#define _unmarked_NLL_OCCURN_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_occuRN(SEXP betaR, 
    SEXP XlamR, SEXP Xlam_offsetR, SEXP XdetR, SEXP Xdet_offsetR,  
    SEXP KR, SEXP yR, SEXP navecR, SEXP nPr, SEXP nOPr) ;

#endif

#ifndef _unmarked_TRANPROBS_H
#define _unmarked_TRANPROBS_H

#include <RcppArmadillo.h>

RcppExport SEXP tranProbs( SEXP Ni1, SEXP Ni2, SEXP omegaR, SEXP gammaR, 
    SEXP deltaR, SEXP dynamicsR) ;

#endif

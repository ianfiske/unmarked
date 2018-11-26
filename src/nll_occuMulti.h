#ifndef _unmarked_NLL_OCCUMULTI_H
#define _unmarked_NLL_OCCUMULTI_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_occuMulti( SEXP fStartR, SEXP fStopR, SEXP dmFr, SEXP dmOccR, 
    SEXP betaR, SEXP dmDetR, SEXP dStartR, SEXP dStopR, SEXP yR, SEXP yStartR, 
    SEXP yStopR, SEXP Iy0r, SEXP zR, SEXP fixed0r) ;

#endif

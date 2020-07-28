#ifndef _unmarked_NLL_GMULTMIX_H
#define _unmarked_NLL_GMULTMIX_H

#include <RcppArmadillo.h>
#include "utils.h"

RcppExport SEXP nll_gmultmix( SEXP betaR, SEXP mixtureR, SEXP pi_funR, 
    SEXP XlamR, SEXP Xlam_offsetR, SEXP XphiR, SEXP Xphi_offsetR, SEXP XdetR, 
    SEXP Xdet_offsetR, SEXP kR, SEXP lfac_kR, SEXP lfac_kmytR, SEXP kmytR, 
    SEXP yR, SEXP naflagR, SEXP finR, SEXP nPr, SEXP nLPr, SEXP nPPr, SEXP nDPr) ;

#endif

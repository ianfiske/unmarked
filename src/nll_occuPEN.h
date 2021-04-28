#ifndef _unmarked_NLL_OCCUPEN_H
#define _unmarked_NLL_OCCUPEN_H

#include <RcppArmadillo.h>
#include <float.h>

RcppExport SEXP nll_occuPEN( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_psiR, SEXP beta_pR, SEXP ndR, SEXP knownOccR, SEXP navecR, SEXP X_offsetR, SEXP V_offsetR, SEXP penaltyR ) ;

#endif

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP getDetVecs(SEXP y_arr, SEXP mp_arr, SEXP J_i, SEXP tin, SEXP K_) ;
extern SEXP getSingleDetVec(SEXP y_, SEXP mp_, SEXP K_);
extern SEXP nll_distsamp( SEXP y_, SEXP lam_, SEXP sig_, SEXP scale_, SEXP a_, SEXP u_, SEXP w_, SEXP db_, SEXP keyfun_, SEXP survey_, SEXP reltol_ );
extern SEXP nll_gpcount( SEXP y_, SEXP Xlam_, SEXP Xphi_, SEXP Xp_, SEXP beta_lam_, SEXP beta_phi_, SEXP beta_p_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xphi_offset_, SEXP Xp_offset_, SEXP M_, SEXP mixture_, SEXP numPrimary_ ) ;
extern SEXP nll_occu( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_psiR, SEXP beta_pR, SEXP ndR, SEXP knownOccR, SEXP navecR, SEXP X_offsetR, SEXP V_offsetR );
extern SEXP nll_occuPEN( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_psiR, SEXP beta_pR, SEXP ndR, SEXP knownOccR, SEXP navecR, SEXP X_offsetR, SEXP V_offsetR, SEXP penaltyR );
extern SEXP nll_pcount( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_lamR, SEXP beta_pR, SEXP log_alphaR, SEXP X_offsetR, SEXP V_offsetR, SEXP naMatR, SEXP lkR, SEXP mixtureR ) ;
extern SEXP nll_pcountOpen( SEXP y_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xp_, SEXP Xiota_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_p_, SEXP beta_iota_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xp_offset_, SEXP Xiota_offset_, SEXP ytna_, SEXP yna_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, SEXP M_, SEXP J_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP fix_, SEXP go_dims_, SEXP immigration_, SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_) ;

static const R_CallMethodDef CallEntries[] = {
    {"getDetVecs",      (DL_FUNC) &getDetVecs,       5},
    {"getSingleDetVec", (DL_FUNC) &getSingleDetVec,  3},
    {"nll_distsamp",    (DL_FUNC) &nll_distsamp,    11},
    {"nll_gpcount",     (DL_FUNC) &nll_gpcount,     14},
    {"nll_occu",        (DL_FUNC) &nll_occu,        10},
    {"nll_occuPEN",     (DL_FUNC) &nll_occuPEN,     11},
    {"nll_pcount",      (DL_FUNC) &nll_pcount,      11},
    {"nll_pcountOpen",  (DL_FUNC) &nll_pcountOpen,  35},
    {NULL, NULL, 0}
};

void R_init_unmarked(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

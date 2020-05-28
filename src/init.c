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
extern SEXP nll_gmultmix( SEXP betaR, SEXP mixtureR, SEXP pi_funR, SEXP XlamR, SEXP Xlam_offsetR, SEXP XphiR, SEXP Xphi_offsetR, SEXP XdetR, SEXP Xdet_offsetR, SEXP kR, SEXP lfac_kR, SEXP lfac_kmytR, SEXP kmytR, SEXP yR, SEXP naflagR, SEXP finR, SEXP nPr, SEXP nLPr, SEXP nPPr, SEXP nDPr ) ;
extern SEXP nll_gdistsamp(SEXP beta_, SEXP mixture_, SEXP keyfun_, SEXP survey_, SEXP Xlam_, SEXP Xlam_offset_, SEXP A_, SEXP Xphi_, SEXP Xphi_offset_, SEXP Xdet_, SEXP Xdet_offset_, SEXP db_, SEXP a_, SEXP u_, SEXP w_, SEXP k_, SEXP lfac_k_, SEXP lfac_kmyt_, SEXP kmyt_, SEXP y_, SEXP naflag_, SEXP fin_, SEXP nP_, SEXP nLP_, SEXP nPP_, SEXP nDP_, SEXP rel_tol_) ;
extern SEXP nll_gpcount( SEXP y_, SEXP Xlam_, SEXP Xphi_, SEXP Xp_, SEXP beta_lam_, SEXP beta_phi_, SEXP beta_p_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xphi_offset_, SEXP Xp_offset_, SEXP M_, SEXP mixture_, SEXP numPrimary_ ) ;
extern SEXP nll_multinomPois( SEXP betaR, SEXP pi_funR, SEXP XlamR, SEXP Xlam_offsetR, SEXP XdetR, SEXP Xdet_offsetR, SEXP yR, SEXP navecR, SEXP nPr, SEXP nAPr ) ;
extern SEXP nll_occu( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_psiR, SEXP beta_pR, SEXP ndR, SEXP knownOccR, SEXP navecR, SEXP X_offsetR, SEXP V_offsetR, SEXP link_psiR);
extern SEXP nll_occuPEN( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_psiR, SEXP beta_pR, SEXP ndR, SEXP knownOccR, SEXP navecR, SEXP X_offsetR, SEXP V_offsetR, SEXP penaltyR );
extern SEXP nll_occuMulti( SEXP fStartR, SEXP fStopR, SEXP dmFr, SEXP dmOccR, SEXP betaR, SEXP dmDetR, SEXP dStartR, SEXP dStopR, SEXP yR, SEXP yStartR, SEXP yStopR, SEXP Iy0r, SEXP zR, SEXP fixed0r);
extern SEXP nll_occuMS( SEXP beta_, SEXP y_, SEXP dm_state_, SEXP dm_phi_, SEXP dm_det_, SEXP sind_, SEXP pind_, SEXP dind_, SEXP prm_, SEXP S_, SEXP T_, SEXP J_, SEXP N_, SEXP naflag_, SEXP guide_);
extern SEXP nll_occuTTD( SEXP beta_, SEXP y_, SEXP delta_, SEXP W_, SEXP V_, SEXP Xgam_, SEXP Xeps_, SEXP pind_, SEXP dind_, SEXP cind_, SEXP eind_, SEXP lpsi_, SEXP tdist_, SEXP N_, SEXP T_, SEXP J_, SEXP naflag_);
extern SEXP get_mlogit(SEXP lp_mat_, SEXP type_, SEXP S_, SEXP guide_);
extern SEXP nll_occuRN( SEXP betaR, SEXP XlamR, SEXP Xlam_offsetR, SEXP XdetR, SEXP Xdet_offsetR, SEXP KR, SEXP yR, SEXP navecR, SEXP nPr, SEXP nOPr );
extern SEXP nll_pcount( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_lamR, SEXP beta_pR, SEXP log_alphaR, SEXP X_offsetR, SEXP V_offsetR, SEXP naMatR, SEXP lkR, SEXP mixtureR ) ;
extern SEXP nll_pcountOpen( SEXP y_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xp_, SEXP Xiota_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_p_, SEXP beta_iota_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xp_offset_, SEXP Xiota_offset_, SEXP ytna_, SEXP yna_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, SEXP M_, SEXP J_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP fix_, SEXP go_dims_, SEXP immigration_, SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_) ;
extern SEXP nll_distsampOpen( SEXP y_, SEXP yt_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xsig_, SEXP Xiota_, SEXP beta_, SEXP beta_ind_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xsig_offset_, SEXP Xiota_offset_, SEXP ytna_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, SEXP first1_, SEXP M_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP survey_, SEXP fix_, SEXP go_dims_, SEXP immigration_, SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_, SEXP a_, SEXP u_, SEXP w_, SEXP db_, SEXP keyfun_, SEXP lfac_k_, SEXP kmyt_, SEXP lfac_kmyt_, SEXP fin_, SEXP A_ ) ;
extern SEXP get_lik_trans(SEXP I_, SEXP I1_);
extern SEXP nll_multmixOpen( SEXP y_, SEXP yt_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_,SEXP Xp_, SEXP Xiota_, SEXP beta_, SEXP beta_ind_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xp_offset_, SEXP Xiota_offset_, SEXP ytna_, SEXP yna_, SEXP lk_, SEXP mixture_,  SEXP first_, SEXP last_, SEXP first1_, SEXP M_, SEXP T_, SEXP J_, SEXP delta_, SEXP dynamics_, SEXP fix_, SEXP go_dims_, SEXP immigration_,  SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_, SEXP pi_fun_, SEXP lfac_k_, SEXP kmyt_, SEXP lfac_kmyt_, SEXP fin_) ;
static const R_CallMethodDef CallEntries[] = {
    {"getDetVecs",      (DL_FUNC) &getDetVecs,       5},
    {"getSingleDetVec", (DL_FUNC) &getSingleDetVec,  3},
    {"nll_distsamp",    (DL_FUNC) &nll_distsamp,    11},
    {"nll_gmultmix",    (DL_FUNC) &nll_gmultmix,    20},
    {"nll_gdistsamp",   (DL_FUNC) &nll_gdistsamp,   27},
    {"nll_gpcount",     (DL_FUNC) &nll_gpcount,     14},
    {"nll_multinomPois",(DL_FUNC) &nll_multinomPois,10},
    {"nll_occu",        (DL_FUNC) &nll_occu,        11},
    {"nll_occuPEN",     (DL_FUNC) &nll_occuPEN,     11},
    {"nll_occuMulti",   (DL_FUNC) &nll_occuMulti,   14},
    {"nll_occuMS",      (DL_FUNC) &nll_occuMS,      15},
    {"nll_occuTTD",     (DL_FUNC) &nll_occuTTD,     17},
    {"get_mlogit",      (DL_FUNC) &get_mlogit,       4},
    {"nll_occuRN",      (DL_FUNC) &nll_occuRN,      10},
    {"nll_pcount",      (DL_FUNC) &nll_pcount,      11},
    {"nll_pcountOpen",  (DL_FUNC) &nll_pcountOpen,  35},
    {"nll_distsampOpen", (DL_FUNC) &nll_distsampOpen, 42},
    {"get_lik_trans", (DL_FUNC) &get_lik_trans,     2},
    {"nll_multmixOpen", (DL_FUNC) &nll_multmixOpen, 38}, 
    {NULL, NULL, 0}
};

void R_init_unmarked(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

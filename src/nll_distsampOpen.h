#ifndef _unmarked_NLL_DISTSAMPOPEN_H
#define _unmarked_NLL_DISTSAMPOPEN_H

#include <RcppArmadillo.h>
#include "tranprobs.h"
#include "distprob.h"
#include "distr.h"

RcppExport SEXP nll_distsampOpen( SEXP y_, SEXP yt_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, 
    SEXP Xsig_, SEXP Xiota_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, 
    SEXP beta_sig_, SEXP beta_iota_, SEXP log_alpha_, SEXP Xlam_offset_, 
    SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xsig_offset_, SEXP Xiota_offset_, 
    SEXP ytna_, SEXP yna_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, 
    SEXP M_, SEXP J_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP survey_, 
    SEXP fix_, SEXP go_dims_, SEXP immigration_, SEXP I_, SEXP I1_, SEXP Ib_, 
    SEXP Ip_, SEXP scale_, SEXP a_, SEXP u_, SEXP w_, SEXP db_, SEXP keyfun_,
    SEXP lfac_k_, SEXP lfac_kmyt_, SEXP kmyt_, SEXP lgy1_, SEXP fin_ ) ;

#endif

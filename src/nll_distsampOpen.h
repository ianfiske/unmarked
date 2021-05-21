#ifndef _unmarked_NLL_DISTSAMPOPEN_H
#define _unmarked_NLL_DISTSAMPOPEN_H

#include <RcppArmadillo.h>
#include <float.h>
#include "tranprobs.h"
#include "distprob.h"
#include "distr.h"

RcppExport SEXP nll_distsampOpen( SEXP y_, SEXP yt_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, 
    SEXP Xsig_, SEXP Xiota_, SEXP beta_, SEXP beta_ind_,
    SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xsig_offset_, 
    SEXP Xiota_offset_, SEXP ytna_, SEXP lk_, SEXP mixture_, SEXP first_, 
    SEXP last_, SEXP first1_, SEXP M_, SEXP T_, SEXP delta_, SEXP dynamics_, 
    SEXP survey_, SEXP fix_, SEXP go_dims_, SEXP immigration_, SEXP I_, 
    SEXP I1_, SEXP Ib_, SEXP Ip_, SEXP a_, SEXP u_, SEXP w_, SEXP db_, 
    SEXP keyfun_, SEXP lfac_k_, SEXP kmyt_, SEXP lfac_kmyt_, SEXP fin_, SEXP A_ ) ;

#endif

#ifndef _unmarked_NLL_PCOUNTOPEN_H
#define _unmarked_NLL_PCOUNTOPEN_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_pcountOpen( SEXP y_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xp_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_p_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xp_offset_, SEXP ytr_, SEXP yr_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, SEXP M_, SEXP J_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP fix_, SEXP go_dims_  ) ;



// constant model
void tp1(arma::mat& g3, int lk, double gam, double om) {
    int Nmin=0;
    for(int n1=0; n1<lk; n1++) {
	for(int n2=0; n2<lk; n2++) {
	    Nmin = std::min(n1, n2);
	    for(int c=0; c<=Nmin; c++) {
		g3.at(n1, n2) += exp(Rf_dbinom(c, n1, om, true) +
				  Rf_dpois(n2-c, gam, true));
	    }
	}
    }
}




// autoregressive model
void tp2(arma::mat& g3, int lk, double gam, double om) {
    int Nmin=0;
    for(int n1=0; n1<lk; n1++) {
	for(int n2=0; n2<lk; n2++) {
	    Nmin = std::min(n1, n2);
	    for(int c=0; c<=Nmin; c++) {
		g3.at(n1, n2) += exp(Rf_dbinom(c, n1, om, true) +
				  Rf_dpois(n2-c, gam*n1, true));
	    }
	}
    }
}



#endif

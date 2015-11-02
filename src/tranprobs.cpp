#include "tranprobs.h"

// constant model
void tp1(arma::mat& g3, int nrI, int nrI1, Rcpp::IntegerVector N, arma::imat I, arma::imat I1, Rcpp::List Ib, Rcpp::List Ip, double gam, double om) {
  Rcpp::NumericVector pois1 = dpois(N, gam, true);
  arma::vec pois = Rcpp::as<arma::vec>(pois1);
  arma::vec bin = arma::zeros<arma::vec>(nrI1);
  for(int i=0; i<nrI1; i++) {
    bin(i) = Rf_dbinom(I1(i,0), I1(i,1), om, true);
  }
  for(int s=0; s<nrI; s++) {
    arma::uvec indB = Rcpp::as<arma::uvec>(Ib[s]);
    arma::uvec indP = Rcpp::as<arma::uvec>(Ip[s]);
    int nc = indB.n_elem;
    for(int q=0; q<nc; q++) {
      g3(s) += exp(bin(indB(q)) + pois(indP(q)));
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







// trend model (exponential growth)
void tp3(arma::mat& g3, int lk, double gam) {
    for(int n1=0; n1<lk; n1++) {
	for(int n2=0; n2<lk; n2++) {
	  g3.at(n1, n2) = Rf_dpois(n2, gam*n1, false);
	}
    }
}


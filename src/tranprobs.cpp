#include "tranprobs.h"


using namespace Rcpp ;



// constant model
void tp1(arma::mat& g3, int nrI, int nrI1, Rcpp::IntegerVector N, arma::imat I, arma::imat I1, Rcpp::List Ib, Rcpp::List Ip, double gam, double om) {
  Rcpp::NumericVector pois1 = dpois(N, gam, true);
  arma::vec pois = as<arma::vec>(pois1);
  arma::vec bin = arma::zeros<arma::vec>(nrI1);
  for(int i=0; i<nrI1; i++) {
    bin(i) = Rf_dbinom(I1(i,0), I1(i,1), om, true);
  }
  for(int s=0; s<nrI; s++) {
    arma::uvec indB = as<arma::uvec>(Ib[s]);
    arma::uvec indP = as<arma::uvec>(Ip[s]);
    int nc = indB.n_elem;
    for(int q=0; q<nc; q++) {
      g3(s) += exp(bin(indB(q)) + pois(indP(q)));
    }
  }
}


// autoregressive + immigration model
void tp2(arma::mat& g3, int lk, double gam, double om, double imm) {
    int Nmin=0;
    for(int n1=0; n1<lk; n1++) {
	for(int n2=0; n2<lk; n2++) {
	    Nmin = std::min(n1, n2);
	    for(int c=0; c<=Nmin; c++) {
		g3.at(n1, n2) += exp(Rf_dbinom(c, n1, om, true) +
				  Rf_dpois(n2-c, gam*n1 + imm, true));
	    }
	}
    }
}
 
// trend + immigration model 
void tp3(arma::mat& g3, int lk, double gam, double imm) {
    for(int n1=0; n1<lk; n1++) {
      for(int n2=0; n2<lk; n2++) {
        g3.at(n1, n2) = Rf_dpois(n2, n1*gam+imm, false);
      }
    }
}





// Ricker + immigration model
void tp4(arma::mat& g3, int lk, double gam, double om, double imm) {
    for(int n1=0; n1<lk; n1++) {
	for(int n2=0; n2<lk; n2++) {
	  g3.at(n1, n2) = Rf_dpois(n2, n1*exp(gam*(1-n1/om)) + imm, false);
	}
    }
}

// Gompertz + immigration model
void tp5(arma::mat& g3, int lk, double gam, double om, double imm) {
    for(int n1=0; n1<lk; n1++) {
	   for(int n2=0; n2<lk; n2++) {
	     g3.at(n1, n2) = Rf_dpois(n2, n1*exp(gam * (1 - log(double (n1) + 1)/log(om + 1))) + imm, false);
	   }
    }
}


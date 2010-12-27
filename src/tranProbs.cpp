#include "tranProbs.h"

using namespace Rcpp ;

SEXP tranProbs( SEXP Kr, SEXP omegaR, SEXP gammaR, SEXP deltaR )
{
    IntegerVector K(Kr);
    int lk = K.size();
    double omega = as<double>(omegaR);
    double gamma = as<double>(gammaR);
    int delta = as<int>(deltaR);
    arma::mat bpsum(lk, lk);
    
    for(int i=0; i<lk; i++) {
        for(int j=0; j<lk; j++) {
            int Km = std::min(K[i], K[j]);
            IntegerVector cmin0 = seq(0, Km);
            bpsum(i, j) = sum(
                dbinom(cmin0, K[j], omega, false) *
                dpois(K[i]-cmin0, gamma, false)
                );
            }
        }
    if(delta > 1) {    
        for(int i=1; i<delta; i++) {
            bpsum *= bpsum;
            }
        }
    return wrap(bpsum);
}


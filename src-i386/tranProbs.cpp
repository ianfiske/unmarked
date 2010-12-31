#include "tranProbs.h"

using namespace Rcpp ;

SEXP tranProbs( SEXP Ni1, SEXP Ni2, SEXP omegaR, SEXP gammaR, SEXP deltaR, 
    SEXP dynamicsR)
{
    IntegerVector N_i1(Ni1);
    IntegerVector N_i2(Ni2);
    int lj = N_i1.size();
    int lk = N_i2.size();
    int N_i2min = N_i2[0];
    double omega = as<double>(omegaR);
    double gamma = as<double>(gammaR);
    int delta = as<int>(deltaR);
    std::string dynamics = as<std::string>(dynamicsR) ;
    
    arma::mat bpsum = arma::zeros(lj, lj);
    for(int j=0; j<lj; j++) {
        for(int k=0; k<lk; k++) {
            int Nm = std::min(N_i1[j], N_i2[k]);
            int counter = k+N_i2min;
            IntegerVector cmin0 = seq(0, Nm);
            if(dynamics=="autoreg")
                bpsum(counter, j) = sum(dbinom(cmin0, N_i1[j], omega, false) *
                    dpois(N_i2[k]-cmin0, gamma*N_i1[j], false));
            else
                bpsum(counter, j) = sum(dbinom(cmin0, N_i1[j], omega, false) *
                    dpois(N_i2[k]-cmin0, gamma, false));
            }
        }
    if(delta > 1) {    
        for(int d=1; d<delta; d++) {
            bpsum *= bpsum;
            arma::mat cs = sum(bpsum, 0);
            arma::mat csm = arma::repmat(cs, lj, 1);
            bpsum = bpsum / csm;
            }
        }
    bpsum = bpsum.rows(N_i2min, lj-1);
    return wrap(bpsum);
}


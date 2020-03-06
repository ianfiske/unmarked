#include "get_lik_trans.h"

using namespace Rcpp;
using namespace arma;

SEXP get_lik_trans(SEXP I_, SEXP I1_){

  umat I = as<umat>(I_);
  umat I1 = as<umat>(I1_);

  List Ib(I.n_rows);
  List Ip(I.n_rows);
  for (unsigned i=0; i<I.n_rows; i++){

    unsigned minI = min(I.row(i));

    IntegerVector Ztmp = seq(0, minI);
    uvec Z = as<uvec>(Ztmp);
   
    uvec Ib_el = find( I1.col(0) <= minI && I1.col(1) == I(i,0) );

    uvec Ip_el = I(i, 1) - Z; 

    Ib[i] = Ib_el;
    Ip[i] = Ip_el;
  }

  List out = List::create(Named("Ib") = Ib , _["Ip"] = Ip);
  return(wrap(out));
}

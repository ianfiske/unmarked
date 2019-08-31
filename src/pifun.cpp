#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

vec removalPiFun ( vec p ){
  int J = p.size();
  vec pi(J);
  pi(0) = p(0);
  for(int j=1; j<J; j++){
    pi(j) = pi(j-1) / p(j-1) * (1-p(j-1)) * p(j);
  }
  return(pi);
}

vec doublePiFun( vec p ){
  //p must have 2 columns
  vec pi(3);
  pi(0) = p(0) * (1 - p(1));
  pi(1) = p(1) * (1 - p(0));
  pi(2) = p(0) * p(1);
  return(pi);
}

vec depDoublePiFun( vec p ){
  //p must have 2 columns
  vec pi(2);
  pi(0) = p(0);
  pi(1) = p(1) * (1 - p(0));
  return(pi);
}

vec piFun( vec p , std::string pi_fun ){
  if(pi_fun == "removalPiFun"){
    return(removalPiFun(p));
  } else if(pi_fun == "doublePiFun"){
    return(doublePiFun(p));
  } else if(pi_fun == "depDoublePiFun"){
    return(depDoublePiFun(p));
  } else {
    stop("Invalid pifun type");
  }
}

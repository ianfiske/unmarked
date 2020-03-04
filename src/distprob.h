#ifndef _unmarked_DISTPROB_H
#define _unmarked_DISTPROB_H

#include <RcppArmadillo.h>
#include "detfuns.h"

using namespace arma;

vec distprob(const std::string& keyfun, const double param1, 
             const double param2, const std::string& survey, const vec& db,
             const vec& w, const rowvec& a) ;

#endif

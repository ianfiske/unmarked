#include "utils.h"

arma::mat inv_logit( arma::mat inp ){
  return(1 / (1 + exp(-1 * inp)));
}

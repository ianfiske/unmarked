# TODO:  update Mh so that it is in the flavor of unmarked.

#' Fits model Mh - closed population model with heterogeneity
#'
#' This fits the logit-normal version of Model Mh in which the logit of
#' individual detection probability is assumed to have a normal distribution with
#' parameters  mu and sigma. These parameters and the unobserved number of
#' individuals, n0, are estimated by integrated likelihood.
#'
#' @title Fit model Mh.
#' @param nx the detection frequency vector, of length J where J is the number of sample occasions
#' @return optim output for now
#' @references
#' Coull and Agresit (1999) Dorazio and Royle (2003)
#' @author
#' Andy Royle \email{aroyle@@usgs.gov}
#' @examples
#' nx<-c(34, 16, 10, 4, 2, 2,0,0,0,0,0,0,0,0)
#' Mh(nx)
#' @keywords models
#' @export
#' @note
#' Assumes balanced data right now.
#' this will be the engine for a poisson reg model for spec richness
Mh <- function(nx) {

  nind<-sum(nx)
  J<-length(nx)

  Mhlik<-function(parms){
    mu<-parms[1]
    sigma<-exp(parms[2])
    n0<-exp(parms[3])

    il<-rep(NA,J+1)
    for(k in 0:J){
      il[k+1]<-integrate(function(x){
                           dbinom(k,J,exp(x)/(1+exp(x)))*dnorm(x,mu,sigma)
                         },
                         lower=-30,upper=30)$value
    }

    -1*(lgamma(n0+nind+1) - lgamma(n0+1) + sum(c(n0,nx)*log(il)))
  }

  fm <- optim(rep(0,3), Mhlik, method="BFGS", hessian = TRUE)

  return(fm)
}



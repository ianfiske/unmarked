# TODO:  update Mh so that it is in the flavor of unmarked.

# Fits model Mh - closed population model with heterogeneity
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



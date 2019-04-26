\name{occuMS}

\alias{occuMS}

\title{Fit Single-Season Multi-State Occupancy Models}

\usage{occuMS(detformulas, stateformulas, data, parameterization="multinomial",
    starts, method="BFGS", se=TRUE, engine=c("C","R"), silent=FALSE, ...)}

\arguments{
    \item{detformulas}{Character vector of formulas for detection probabilities.
      See details for a description of how to order these formulas.}
    \item{stateformulas}{Character vector of formulas for occupancy probabilities. 
      See details for a description of how to order these formulas.}
    \item{data}{An \code{\link{unmarkedFrameOccuMS}} object}
    \item{parameterization}{Either \code{"multinomial"} for the multinomial 
      parameterization (MacKenzie et al. 2009) which allows an arbitrary
      number of occupancy states, or \code{"condbinom"} for the conditional 
      binomial parameterization (Nichols et al. 2007) which requires exactly
      3 occupancy states. See details.}
    \item{starts}{Vector of parameter starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{Logical specifying whether or not to compute standard
      errors.}
    \item{engine}{Either "C" to use fast C++ code or "R" to use native R
      code during the optimization.}
    \item{silent}{Boolean; if \code{TRUE}, suppress warnings.}
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{This function fits single-season multi-state occupancy models with both the multinomial and conditional binomial parameterizations.}

\details{

See \code{\link{unmarkedFrame}} and \code{\link{unmarkedFrameOccuMS}} for a
description of how to supply data to the \code{data} argument.


}

\value{unmarkedFitOccuMS object describing the model fit.}

\references{
MacKenzie, D. I., Nichols, J. D., Seamans, M. E., and R. J. Gutierrez, 2009. 
  Modeling species occurrence dynamics with multiple states and imperfect 
  detection. Ecology 90: 823-835.

Nichols, J. D., Hines, J. E., Mackenzie, D. I., Seamans, M. E., and 
  R. J. Gutierrez. 2007. Occupancy estimation and modeling with multiple states 
  and state uncertainty. Ecology 88: 1395-1400.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuMS}}}


\examples{

\dontrun{

#Multinomial parameterization
#Simulate data

#Parameters
N <- 500; J <- 5; S <- 3
site_covs <- matrix(rnorm(N*2),ncol=2)
obs_covs <- matrix(rnorm(N*J*2),ncol=2)
a1 <- -0.5; b1 <- 1; a2 <- -0.6; b2 <- -0.7
p11 <- -0.4; p12 <- -1.09; p22 <- -0.84
truth <- c(a1,b1,a2,b2,p11,0,p12,p22)

#State process
lp <- matrix(NA,ncol=S,nrow=N)
for (n in 1:N){
  lp[n,2] <- exp(a1+b1*site_covs[n,1])
  lp[n,3] <- exp(a2+b2*site_covs[n,2])
  lp[n,1] <- 1  
}
psi_mat <- lp/rowSums(lp)

z <- rep(NA,N)
for (n in 1:N){
  z[n] <- sample(0:2, 1, replace=T, prob=psi_mat[n,])
}

probs_raw <- matrix(c(1,0,0,1,exp(p11),0,1,exp(p12),exp(p22)),nrow=3,byrow=T)
probs_raw <- probs_raw/rowSums(probs_raw)
  
y <- matrix(0,nrow=N,ncol=J)
for (n in 1:N){

  probs <- switch(z[n]+1,
                  probs_raw[1,],
                  probs_raw[2,],
                  probs_raw[3,])
  if(z[n]>0){
    y[n,] <- sample(0:2, J, replace=T, probs)
  }
}

#Construct unmarkedFrame
umf <- unmarkedFrameOccuMS(y=y,siteCovs=as.data.frame(site_covs),
                           obsCovs=as.data.frame(obs_covs))

#Formulas

#3 states, so detformulas is a character vector of formulas of 
#length 3 in following order:
#1) p[11]: prob of detecting state 1 given true state 1
#2) p[12]: prob of detecting state 1 given true state 2
#3) p[22]: prob of detecting state 2 given true state 2
detformulas <- c('~V1','~1','~1')
#If you had 4 states, it would be p[11],p[12],p[13],p[22],p[23],p[33] and so on

#3 states, so stateformulas is a character vector of length 2 in following order:
#1) psi[1]: probability of state 1
#2) psi[2]: probability of state 2
#You can get probability of state 0 (unoccupied) as 1 - psi[1] - psi[2]
stateformulas <- c('~V1','~V2')

#Fit model
fit <- occuMS(detformulas, stateformulas, umf, parameterization="multinomial")

#Look at results
fit
#Compare with truth
cbind(truth=truth,estimate=coef(fit))

#Generate predicted values
lapply(predict(fit,type='state'),head)
lapply(predict(fit,type='det'),head)

#Fit a null model
detformulas <- rep('~1',3)
stateformulas <- rep('~1',2)
fit_null <- occuMS(detformulas, stateformulas, umf, parameterization="multinomial")

#Compare fits
modSel(fitList(fit,fit_null))
}

}

\keyword{models}

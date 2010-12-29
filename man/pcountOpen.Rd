\name{pcountOpen}
\alias{pcountOpen}
\title{Fit the open N-mixture model of Dail and Madsen}
\description{Fit the model of Dail and Madsen, which is 
a generalized form of the Royle (2004) N-mixture model allowing for
estimation of recruitment (births/imigrants) and 
apparent survival (mortality/emigration).}
\usage{
pcountOpen(lambdaformula, gammaformula, omegaformula, pformula, data, 
	mixture = c("P", "NB"), K, fix=c("none", "gamma", "omega"), starts, 
    method = "BFGS", se = TRUE, ...)
}
\arguments{
  \item{lambdaformula}{Right-hand sided formula for initial abundance}
  \item{gammaformula}{Right-hand sided formula for recruitment rate}
  \item{omegaformula}{Right-hand sided formula for survival probability}
  \item{pformula}{Right-hand sided formula for detection probability}
  \item{data}{An object of class unmarkedFramePCountOpen. See details}
  \item{mixture}{character specifying mixture: either "P" or "NB".}
  \item{K}{Integer upper index of integration for N-mixture.}
  \item{fix}{If "omega", omega is fixed at 1. If "gamma", gamma is fixed at 0.}
  \item{starts}{vector of starting values}
  \item{method}{Optimization method used by \code{\link{optim}}.}
  \item{se}{logical specifying whether or not to compute standard errors.}
  \item{\dots}{additional arguments to be passed to \code{\link{optim}}.}
}
\details{
This model generalizes the Royle (2004) N-mixture model by relaxing the 
closure assumption.

The latent initial abundance distribution, 
\eqn{f(N | \mathbf{\theta})}{f(N | theta)} can be set as either a 
Poisson or a negative binomial random variable, depending on the setting of the 
\code{mixture} argument. \code{mixture = "P"} or \code{mixture = "NB"} select 
the Poisson or negative binomial distribution respectively.  The mean of 
\eqn{N_i} is \eqn{\lambda_i}{lambda_i}.  If \eqn{N_i \sim NB}{N_i ~ NB}, then an
additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
\eqn{\alpha}{alpha} implies higher variance).

The latent abundance state following the initial sampling period arises from a
Markovian process in which survivors are modeled as \eqn{S_it \sim 
Binomial(N_it-1, omega_it)}{S_it ~ Binomial(N_it-1, omega_it)}, and recruits
follow \eqn{G_it \sim Poisson(gamma_it)}{G_it ~ Poisson(gamma_it)}.

The detection process is modeled as binomial: \eqn{y_{it} \sim
Binomial(N_it, p_it)}{y_it ~ Binomial(N_it, p_it)}.

Covariates of \eqn{\lambda_i}{lamdba_i} use the log link and
covariates of \eqn{p_it}{p_it} use the logit link.
 
}
\value{
An unmarkedFitPCO object.
}
\references{
Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Dail, D. and L. Madsen (In press) Models for Estimating Abundance from Repeated Counts of an Open Metapopulation. \emph{Biometrics}.

}

\author{Richard Chandler \email{richard.chandlers@gmail.com}}

\note{
When gamma or omega are modeled using observation-level covariates, the covariate data for the final survey occasion will be ignored; however, they must be supplied.

If time intervals vary among survey occasions (sampling periods), an M by T matrix of dates (integer values) should be supplied to unmarkedFramePCO using the dates argument. Currently, the handling of non-constant survey intervals differs from the methods described in the original paper. 

At this time, the robust design is not available. However, if you have multiple seasons of data, you can model within vs. among rates by creating an MxT observation-level covariate identifying the transition periods and include it the gamma and omega formulas.
}

\section{Warning}{This function can be extremely slow. Consider testing the timing on a small subset of the data, perhaps with se=FALSE. Finding the lowest value of K that does not affect estimates will also help with speed. }

\seealso{
\code{\link{pcount}}
}

\examples{

## Simulate no covariates, constant sampling period intervals	
set.seed(3)
M <- 50
T <- 5
lambda <- 1
gamma <- 0.5
omega <- 0.8
p <- 0.7
y <- N <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], omega)
	G[,t] <- rpois(M, gamma)
	N[,t+1] <- S[,t] + G[,t]
	}
y[] <- rbinom(M*T, N, p)
# y[M-5, 1:4] <- y[M-4, 2:5] <- y[M-3, 2:4] <- y[M-2, c(2,4)] <- y[M-1, 1] <- y[M, T] <- NA

                            
# Prepare data                               
umf <- unmarkedFramePCO(y = y)

summary(umf)

# Fit model and backtransform
(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10)) # 14s on 64bit.
backTransform(m1, "lambda")
backTransform(m1, "gamma")
backTransform(m1, "omega")
backTransform(m1, "det")


}

\keyword{models}

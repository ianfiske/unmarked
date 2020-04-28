\name{pcountOpen}
\alias{pcountOpen}
\title{Fit the open N-mixture models of Dail and Madsen and extensions}
\description{Fit the models of Dail and Madsen (2011) and Hostetler and
Chandler (in press), which are
generalized forms of the Royle (2004) N-mixture model for open
populations.}
\usage{
pcountOpen(lambdaformula, gammaformula, omegaformula, pformula,
  data, mixture = c("P", "NB", "ZIP"), K, dynamics=c("constant", "autoreg",
  "notrend", "trend", "ricker", "gompertz"), fix=c("none", "gamma", "omega"),
  starts, method = "BFGS", se = TRUE, immigration = FALSE,
  iotaformula = ~1, ...)
}
\arguments{
  \item{lambdaformula}{
    Right-hand sided formula for initial abundance
  }
  \item{gammaformula}{
    Right-hand sided formula for recruitment rate (when dynamics is "constant",
    "autoreg", or "notrend") or population growth rate (when dynamics is
    "trend", "ricker", or "gompertz")
  }
  \item{omegaformula}{
    Right-hand sided formula for apparent survival probability
    (when dynamics is "constant", "autoreg", or "notrend") or equilibrium
    abundance (when dynamics is "ricker" or "gompertz")
  }
  \item{pformula}{
    Right-hand sided formula for detection probability
  }
  \item{data}{
    An object of class \code{\link{unmarkedFramePCO}}. See details
  }
  \item{mixture}{
    character specifying mixture: "P", "NB", or "ZIP" for
    the Poisson, negative binomial, and zero-inflated Poisson
    distributions.
  }
  \item{K}{
    Integer defining upper bound of discrete integration. This
    should be higher than the maximum observed count and high enough
    that it does not affect the parameter estimates. However, the higher
    the value the slower the compuatation.
  }
  \item{dynamics}{
    Character string describing the type of population
    dynamics. "constant" indicates that there is no relationship between
    omega and gamma. "autoreg" is an auto-regressive model in which
    recruitment is modeled as gamma*N[i,t-1]. "notrend" model gamma as
    lambda*(1-omega) such that there is no temporal trend. "trend" is
    a model for exponential growth, N[i,t] = N[i,t-1]*gamma, where gamma
    in this case is finite rate of increase (normally referred to as
    lambda). "ricker" and "gompertz" are models for density-dependent
    population growth.  "ricker" is the Ricker-logistic model, N[i,t] =
    N[i,t-1]*exp(gamma*(1-N[i,t-1]/omega)), where gamma is the maximum
    instantaneous population growth rate (normally referred to as r) and
    omega is the equilibrium abundance (normally referred to as K).  "gompertz"
    is a modified version of the Gompertz-logistic model, N[i,t] =
    N[i,t-1]*exp(gamma*(1-log(N[i,t-1]+1)/log(omega+1))), where the
    interpretations of gamma and omega are similar to in the Ricker model.
  }
  \item{fix}{
    If "omega", omega is fixed at 1. If "gamma", gamma is fixed at 0.
  }
  \item{starts}{
    vector of starting values
  }
  \item{method}{
    Optimization method used by \code{\link{optim}}.
  }
  \item{se}{
    logical specifying whether or not to compute standard errors.
  }
  \item{immigration}{
    logical specifying whether or not to include an immigration term (iota) in
    population dynamics.
  }
  \item{iotaformula}{
    Right-hand sided formula for average number of immigrants to a site
    per time step
  }
  \item{\dots}{
    additional arguments to be passed to \code{\link{optim}}.
  }
}
\details{
These models generalize the Royle (2004) N-mixture model by relaxing the
closure assumption. The models include two or three additional parameters:
gamma, either the recruitment rate (births and immigrations), the
finite rate of increase, or the maximum instantaneous rate of increase;
omega, either the apparent survival rate (deaths and emigrations) or the
equilibrium abundance (carrying capacity); and iota, the number of immigrants
per site and year. Estimates of
population size at each time period can be derived from these
parameters, and thus so can trend estimates. Or, trend can be estimated
directly using dynamics="trend".

When immigration is set to FALSE (the default), iota is not modeled.
When immigration is set to TRUE and dynamics is set to "autoreg", the model
will separately estimate birth rate (gamma) and number of immigrants (iota).
When immigration is set to TRUE and dynamics is set to "trend", "ricker", or
"gompertz", the model will separately estimate local contributions to
population growth (gamma and omega) and number of immigrants (iota).

The latent abundance distribution, \eqn{f(N | \mathbf{\theta})}{f(N |
theta)} can be set as a Poisson, negative binomial, or zero-inflated
Poisson random
variable, depending on the setting of the \code{mixture} argument,
\code{mixture = "P"}, \code{mixture = "NB"}, \code{mixture = "ZIP"}
respectively.  For the first two distributions, the mean of \eqn{N_i} is
\eqn{\lambda_i}{lambda_i}.  If \eqn{N_i \sim NB}{N_i ~ NB}, then an
additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
\eqn{\alpha}{alpha} implies higher variance). For the ZIP distribution,
the mean is \eqn{\lambda_i(1-\psi)}{lambda_i*(1-psi)}, where psi is the
zero-inflation parameter.

For "constant", "autoreg", or "notrend" dynamics, the latent abundance state
following the initial sampling period arises
from a
Markovian process in which survivors are modeled as \eqn{S_{it} \sim
Binomial(N_{it-1}, \omega_{it})}{S(i,t) ~ Binomial(N(i,t-1),
omega(i,t))}, and recruits
follow \eqn{G_{it} \sim Poisson(\gamma_{it})}{G(i,t) ~
  Poisson(gamma(i,t))}.
Alternative population dynamics can be specified
using the \code{dynamics} and \code{immigration} arguments.

The detection process is modeled as binomial: \eqn{y_{ijt} \sim
Binomial(N_{it}, p_{ijt})}{y(i,j,t) ~ Binomial(N(i,t), p(i,j,t))}.

\eqn{\lambda_i}{lambda_i}, \eqn{\gamma_{it}}{gamma_it}, and
\eqn{\iota_{it}}{iota_it} are modeled 
using the the log link.
\eqn{p_{ijt}}{p_ijt} is modeled using
the logit link.
\eqn{\omega_{it}}{omega_it} is either modeled using the logit link (for
"constant", "autoreg", or "notrend" dynamics) or the log link (for "ricker"
or "gompertz" dynamics).  For "trend" dynamics, \eqn{\omega_{it}}{omega_it}
is not modeled.
}
\value{
An object of class unmarkedFitPCO.
}
\references{
Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from
Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Dail, D. and L. Madsen (2011) Models for Estimating Abundance from
Repeated Counts of an Open Metapopulation. \emph{Biometrics}. 67, pp
577-587.

Hostetler, J. A. and R. B. Chandler (2015) Improved State-space Models for
Inference about Spatial and Temporal Variation in Abundance from Count Data.
\emph{Ecology}
96:1713-1723.

}

\author{Richard Chandler \email{rbchan@uga.edu} and Jeff Hostetler}

\note{
When gamma or omega are modeled using year-specific covariates, the
covariate data for the final year will be ignored; however,
they must be supplied.

If the time gap between primary periods is not constant, an M by T
matrix of integers should be supplied to \code{\link{unmarkedFramePCO}}
using the \code{primaryPeriod} argument.

Secondary sampling periods are optional, but can greatly improve the
precision of the estimates.
}

\section{Warning}{This function can be extremely slow, especially if
  there are covariates of gamma or omega. Consider testing the timing on
  a small subset of the data, perhaps with se=FALSE. Finding the lowest
  value of K that does not affect estimates will also help with speed. }

\seealso{
\code{\link{pcount}, \link{unmarkedFramePCO}}
}

\examples{

## Simulation
## No covariates, constant time intervals between primary periods, and
## no secondary sampling periods

set.seed(3)
M <- 50
T <- 5
lambda <- 4
gamma <- 1.5
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


# Prepare data
umf <- unmarkedFramePCO(y = y, numPrimary=T)
summary(umf)


# Fit model and backtransform
(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20)) # Typically, K should be higher

(lam <- coef(backTransform(m1, "lambda"))) # or
lam <- exp(coef(m1, type="lambda"))
gam <- exp(coef(m1, type="gamma"))
om <- plogis(coef(m1, type="omega"))
p <- plogis(coef(m1, type="det"))

\dontrun{
# Finite sample inference. Abundance at site i, year t
re <- ranef(m1)
devAskNewPage(TRUE)
plot(re, layout=c(5,5), subset = site \%in\% 1:25 & year \%in\% 1:2,
     xlim=c(-1,15))
devAskNewPage(FALSE)

(N.hat1 <- colSums(bup(re)))

# Expected values of N[i,t]
N.hat2 <- matrix(NA, M, T)
N.hat2[,1] <- lam
for(t in 2:T) {
    N.hat2[,t] <- om*N.hat2[,t-1] + gam
    }

rbind(N=colSums(N), N.hat1=N.hat1, N.hat2=colSums(N.hat2))


}

}

\keyword{models}

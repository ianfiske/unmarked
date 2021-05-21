\name{multmixOpen}
\alias{multmixOpen}
\title{
  Open population multinomial N-mixture model
}

\description{
  Fit the model of Dail and Madsen (2011) and Hostetler and Chandler
  (2015) for designs involving repeated counts that yield multinomial outcomes.
  Possible data collection methods include repeated removal sampling and 
  double observer sampling.
}

\usage{
multmixOpen(lambdaformula, gammaformula, omegaformula, pformula,
    data, mixture=c("P", "NB", "ZIP"), K,
    dynamics=c("constant", "autoreg", "notrend", "trend", "ricker", "gompertz"),
    fix=c("none", "gamma", "omega"), immigration=FALSE, iotaformula = ~1,
    starts, method="BFGS", se=TRUE, ...)
}

\arguments{
  \item{lambdaformula}{Right-hand sided formula for initial abundance}
  \item{gammaformula}{Right-hand sided formula for recruitment rate (when 
    dynamics is "constant", "autoreg", or "notrend") or population growth rate 
    (when dynamics is "trend", "ricker", or "gompertz")}
  \item{omegaformula}{Right-hand sided formula for apparent survival probability
    (when dynamics is "constant", "autoreg", or "notrend") or equilibrium
    abundance (when dynamics is "ricker" or "gompertz")}
  \item{pformula}{A right-hand side formula describing the detection 
    function covariates}
  \item{data}{An object of class \code{\link{unmarkedFrameMMO}}}
  \item{mixture}{String specifying mixture: "P", "NB", or "ZIP" for
    the Poisson, negative binomial, or zero-inflated Poisson
    distributions respectively}
  \item{K}{Integer defining upper bound of discrete integration. This
    should be higher than the maximum observed count and high enough
    that it does not affect the parameter estimates. However, the higher
    the value the slower the computation}
  \item{dynamics}{Character string describing the type of population
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
    interpretations of gamma and omega are similar to in the Ricker model}
  \item{fix}{If "omega", omega is fixed at 1. If "gamma", gamma is fixed at 0}
  \item{immigration}{Logical specifying whether or not to include an immigration 
    term (iota) in population dynamics}
  \item{iotaformula}{Right-hand sided formula for average number of immigrants 
    to a site per time step}
  \item{starts}{Vector of starting values}
  \item{method}{Optimization method used by \code{\link{optim}}}
  \item{se}{Logical specifying whether or not to compute standard errors}
  \item{\dots}{Additional arguments to optim, such as lower and upper bounds}
}

\details{

These models generalize multinomial N-mixture models (Royle et al. 2004) by
relaxing the closure assumption (Dail and Madsen 2011, Hostetler and Chandler
2015, Sollmann et al. 2015).

The models include two or three additional parameters:
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

\eqn{\lambda_i}{lambda_i}, \eqn{\gamma_{it}}{gamma_it}, and
\eqn{\iota_{it}}{iota_it} are modeled 
using the the log link.
\eqn{p_{ijt}}{p_ijt} is modeled using
the logit link.
\eqn{\omega_{it}}{omega_it} is either modeled using the logit link (for
"constant", "autoreg", or "notrend" dynamics) or the log link (for "ricker"
or "gompertz" dynamics).  For "trend" dynamics, \eqn{\omega_{it}}{omega_it}
is not modeled.

The detection process is modeled as multinomial:
\eqn{\mathbf{y_{it}} \sim
Multinomial(N_{it}, \pi_{it})}{y(i,1:J,t) ~
Multinomial(N(i,t), pi(i,1,t), pi(i,2,t), ..., pi(i,J,t))},
where \eqn{\pi_{ijt}}{pi(ijt)} is the multinomial cell probability for
plot i at time t on occasion j.

Options for the detection process include equal-interval removal sampling
(\code{"removal"}), double observer sampling (\code{"double"}), or 
dependent double-observer sampling (\code{"depDouble"}). This option is
specified when setting up the data using \code{\link{unmarkedFrameMMO}}. Note
that unlike the related functions \code{\link{multinomPois}} and
\code{\link{gmultmix}}, custom functions for the detection process (i.e., 
\code{piFun}s) are not supported. To request additional options contact the author.
}

\value{An object of class unmarkedFitMMO}

\references{
  Dail, D. and L. Madsen (2011) Models for Estimating Abundance from
  Repeated Counts of an Open Metapopulation. \emph{Biometrics}. 67: 577-587.

  Hostetler, J. A. and R. B. Chandler (2015) Improved State-space Models for
  Inference about Spatial and Temporal Variation in Abundance from Count Data.
  \emph{Ecology} 96: 1713-1723.

  Royle, J. A. (2004). Generalized estimators of avian abundance from
  count survey data. \emph{Animal Biodiversity and Conservation} 27(1), 375-386.
}

\author{Ken Kellner \email{contact@kenkellner.com}, Richard Chandler}

\note{
  When gamma or omega are modeled using year-specific covariates, the
  covariate data for the final year will be ignored; however,
  they must be supplied.

  If the time gap between primary periods is not constant, an M by T
  matrix of integers should be supplied to \code{\link{unmarkedFrameMMO}}
  using the \code{primaryPeriod} argument.

  Secondary sampling periods are optional, but can greatly improve the
  precision of the estimates.
}

\section{Warning}{This function can be extremely slow, especially if
  there are covariates of gamma or omega. Consider testing the timing on
  a small subset of the data, perhaps with se=FALSE. Finding the lowest
  value of K that does not affect estimates will also help with speed. }

\seealso{
\code{\link{multinomPois}, \link{gmultmix}, \link{unmarkedFrameMMO}}
}

\examples{
   
  #Generate some data 
  set.seed(123)
  lambda=4; gamma=0.5; omega=0.8; p=0.5
  M <- 100; T <- 5
  y <- array(NA, c(M, 3, T))
  N <- matrix(NA, M, T)
  S <- G <- matrix(NA, M, T-1)

  for(i in 1:M) {
    N[i,1] <- rpois(1, lambda)
    y[i,1,1] <- rbinom(1, N[i,1], p)    # Observe some
    Nleft1 <- N[i,1] - y[i,1,1]         # Remove them
    y[i,2,1] <- rbinom(1, Nleft1, p)   # ...
    Nleft2 <- Nleft1 - y[i,2,1]
    y[i,3,1] <- rbinom(1, Nleft2, p)

    for(t in 1:(T-1)) {
      S[i,t] <- rbinom(1, N[i,t], omega)
      G[i,t] <- rpois(1, gamma)
      N[i,t+1] <- S[i,t] + G[i,t]
      y[i,1,t+1] <- rbinom(1, N[i,t+1], p)    # Observe some
      Nleft1 <- N[i,t+1] - y[i,1,t+1]         # Remove them
      y[i,2,t+1] <- rbinom(1, Nleft1, p)   # ...
      Nleft2 <- Nleft1 - y[i,2,t+1]
      y[i,3,t+1] <- rbinom(1, Nleft2, p)
    }
  }
  y=matrix(y, M)
  
  #Create some random covariate data
  sc <- data.frame(x1=rnorm(100))

  \dontrun{
  #Create unmarked frame
  umf <- unmarkedFrameMMO(y=y, numPrimary=5, siteCovs=sc, type="removal")

  #Fit model
  (fit <- multmixOpen(~x1, ~1, ~1, ~1, K=30, data=umf))
  
  #Compare to truth
  cf <- coef(fit)
  data.frame(model=c(exp(cf[1]), cf[2], exp(cf[3]), plogis(cf[4]), plogis(cf[5])), 
             truth=c(lambda, 0, gamma, omega, p))

  #Predict
  head(predict(fit, type='lambda'))

  #Check fit with parametric bootstrap
  pb <- parboot(fit, nsims=15)
  plot(pb)

  # Empirical Bayes estimates of abundance for each site / year
  re <- ranef(fit)
  plot(re, layout=c(10,5), xlim=c(-1, 10))
  }

}

\keyword{models}

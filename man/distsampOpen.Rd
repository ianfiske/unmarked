\name{distsampOpen}
\alias{distsampOpen}
\title{
  Open population model for distance sampling data
}

\description{
  Fit the model of Dail and Madsen (2011) and Hostetler and Chandler
  (2015) with a distance sampling observation model (Sollmann et al. 2015). 
}

\usage{
distsampOpen(lambdaformula, gammaformula, omegaformula, pformula,
    data, keyfun=c("halfnorm", "exp", "hazard", "uniform"),
    output=c("abund", "density"), unitsOut=c("ha", "kmsq"),
    mixture=c("P", "NB", "ZIP"), K,
    dynamics=c("constant", "autoreg", "notrend", "trend"),
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
  \item{data}{An object of class \code{\link{unmarkedFrameDSO}}}
  \item{keyfun}{One of the following detection functions: "halfnorm", "hazard", 
    "exp", or "uniform"}
  \item{output}{Model either "density" or "abund"}
  \item{unitsOut}{Units of density. Either "ha" or "kmsq" for hectares and 
    square kilometers, respectively}
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

These models generalize distance sampling models (Buckland et al. 2001) by
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

For the distance sampling detection process, half-normal (\code{"halfnorm"}), 
exponential (\code{"exp"}), hazard (\code{"hazard"}), and uniform
(\code{"uniform"}) key functions are available.
}

\value{An object of class unmarkedFitDSO}

\references{
  Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L. 
  and Thomas, L. (2001) \emph{Introduction to Distance Sampling: Estimating 
  Abundance of Biological Populations}. Oxford University Press, Oxford, UK.
  
  Dail, D. and L. Madsen (2011) Models for Estimating Abundance from
  Repeated Counts of an Open Metapopulation. \emph{Biometrics}. 67: 577-587.

  Hostetler, J. A. and R. B. Chandler (2015) Improved State-space Models for
  Inference about Spatial and Temporal Variation in Abundance from Count Data.
  \emph{Ecology} 96: 1713-1723.

  Sollmann, R., Gardner, B., Chandler, R.B., Royle, J.A. and Sillett, T.S. 
  (2015) An open-population hierarchical distance sampling model. 
  \emph{Ecology} 96: 325-331.
}

\author{Richard Chandler, Jeff Hostetler, Andy Royle, Ken Kellner}

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
\code{\link{distsamp}, \link{gdistsamp}, \link{unmarkedFrameDSO}}
}

\examples{

sim7 <- function(lambda=1, gamma=0.5, omega=0.8, sigma=40, scale=NULL, M=100, T=5, J=4,type="line",
keyfun="halfnorm")
    {
    y <- array(NA, c(M, J, T))
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    db <- c(0, 25, 50, 75, 100)
    if(length(db)-1 != J)
      stop("hey, what")
if(keyfun=="halfnorm"){
    if(type=="point")
       g <- function(x, sig) exp(-x^2/(2*sig^2))*x
    if(type=="line")
       g <- function(x, sig) exp(-x^2/(2*sig^2))
}
if(keyfun=="exp"){
    if(type=="point")
       g <- function(x, sig) exp(-x/sig)*x
    if(type=="line")
       g <- function(x, sig) exp(-x/sig)
}
if(keyfun=="hazard"){
    shape<- sigma
    if(type=="point")
       g <- function(x, shape, scale) (1 - exp(-(r/shape)^-scale)) * x
    if(type=="line")
       g <- function(x, sig) (1 - exp(-(x/shape)^-scale) )
}
if(keyfun=="uniform"){
     if(type=="point")
       g <-  function(x,sig) x
    if(type=="line")
       g <- function(x,sig) 1
}

    cp <- u <- a <- numeric(J)
if(keyfun!="uniform"){
if(type=="point"){
   a[1] <- pi*db[2]^2
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value * 2 * pi
    for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value * 2*pi
    }
    }
if(type=="line"){
   L <-  1
   a[1] <- L*db[2]
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value
    for(j in 2:J) {
      a[j] <-  db[j+1]  - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value
    }
    }
    u <- a / sum(a)
    cp <- cp / a * u
    cp[j+1] <- 1-sum(cp)

}
if(keyfun=="uniform"){
if(type=="point"){
   a[1] <- pi*db[2]^2
   for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])
 
    }
    }
if(type=="line"){
   L <-  1
   a[1] <- L*db[2]
 
    for(j in 2:J) {
      a[j] <-  db[j+1]  - sum(a[1:j])
 
    }
    }
    u <- a / sum(a)
    cp <- a
    cp[j+1] <- 0

}
    for(i in 1:M) {
    N[i,1] <- rpois(1, lambda)

    y[i,1:J,1] <- rmultinom(1, N[i,1], cp)[1:J]

    for(t in 1:(T-1)) {
        S[i,t] <- rbinom(1, N[i,t], omega)
        G[i,t] <- rpois(1, gamma)
        N[i,t+1] <- S[i,t] + G[i,t]
        y[i,1:J,t+1] <- rmultinom(1, N[i,t+1], cp)[1:J]
        }
    }
    cp <- array(cp, c(J, M, T))
    cp <- matrix(aperm(cp, c(2,1,3)), M)
    cat("max(N) =", max(N), "\n")
    return(list(y=matrix(y, M),N=N))
}
library(unmarked)
set.seed(711)
lambda <- 4
gamma <- 2
omega <- 0.5
sigma <- 20
T <- 20
 
# half-normal examples
y.sim7 <- sim7(lambda, gamma, omega, sigma=15, M=200, T=T,type="line", keyfun="halfnorm")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,keyfun="halfnorm",
                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 #starts=c(0,0,0,2.21),
                 se=TRUE)

# exponential example 
y.sim7 <- sim7(lambda, gamma, omega, sigma=25, M=200, T=T,type="line", keyfun="exp")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))

fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,keyfun="exp",
                 #starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 starts=c(0,0,0,2),
                   se=TRUE)

# uniform example 
y.sim7 <- sim7(lambda, gamma, omega, sigma=NULL, M=200, T=T,type="line", keyfun="uniform", scale=1)$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,
                   starts=rep(0,3),
                   keyfun="uniform", se=TRUE)
}

\keyword{models}

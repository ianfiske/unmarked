\name{gdistsamp}
\alias{gdistsamp}
\title{
Fit the generalized distance sampling model of Chandler et al. (2011).
}
\description{
This model can be fit to distance sampling data collected at i = 1, 2, ..., R
plots surveyed on t = 1, 2, ..., T occasions. The detections must be binned
into distance classes.
}
\usage{
gdistsamp(lambdaformula, phiformula, pformula, data, keyfun = c("halfnorm", "exp", "hazard", "uniform"), output = c("abund", "density"), unitsOut = c("ha", "kmsq"), mixture = c("P", "NB"), K, starts, method = "BFGS", control = list(), se = TRUE, rel.tol=1e-4)
}
\arguments{
  \item{lambdaformula}{
      A right-hand side formula describing the abundance covariates.
}
  \item{phiformula}{
      A right-hand side formula describing the availability covariates.
}
  \item{pformula}{
      A right-hand side formula describing the detection function covariates.
}
  \item{data}{
      An object of class \code{unmarkedFrameGDS}
}
  \item{keyfun}{
      One of the following detection functions: "halfnorm", "hazard", "exp",
          or "uniform." See details.
}
  \item{output}{
      Model either "density" or "abund"
}
  \item{unitsOut}{
      Units of density. Either "ha" or "kmsq" for hectares and square
          kilometers, respectively.
}
  \item{mixture}{
      Either "P" or "NB" for the Poisson and negative binomial models of
          abundance.
}
  \item{K}{
      An integer value specifying the upper bound used in the integration.
}
  \item{starts}{
      A numeric vector of starting values for the model parameters.
}
  \item{method}{
      Optimization method used by \code{\link{optim}}.
}
  \item{control}{
      Other arguments passed to \code{\link{optim}}.
}
  \item{se}{
      logical specifying whether or not to compute standard errors.
}
\item{rel.tol}{relative accuracy for the integration of the detection function.
    See \link{integrate}. You might try adjusting this if you get an error
    message related to the integral.}
}
\details{
    Unlike conventional distance sampling, which uses the 'conditional on
    detection' likelihood formulation, this model is based upon the
    unconditional likelihood and thus allows for modeling both abundance and
    detection function parameters.

    The latent transect-level abundance distribution
    \eqn{f(N | \mathbf{\theta})}{f(N | theta)} is currently assumed to be
    Poisson with mean \eqn{\lambda}{lambda}.

    The detection process is modeled as multinomial:
    \eqn{y_{ij} \sim Multinomial(N_i, pi_{ij})}{y_ij ~ Multinomial(N_i, pi_i1, pi_i2, ..., pi_iJ)},
    where \eqn{pi_ij} is the multinomial cell probability for transect i in
    distance class j. These are computed based upon a detection function
    \eqn{g(x | \mathbf{\sigma})}{g(x | sigma)}, such as the half-normal,
    negative exponential, or hazard rate.

    Parameters \eqn{\lambda}{lambda} and \eqn{\sigma}{sigma} can be vectors
    affected by transect-specific covariates using the log link.
    }

\value{
    An object of class unmarkedFitGDS.
    }
\references{
    Royle, J. A., D. K. Dawson, and S. Bates (2004) Modeling
    abundance effects in distance sampling. \emph{Ecology} 85, pp. 1591-1597.
    }
\author{
    Richard Chandler \email{rchandler@usgs.gov}
    }
\note{
    You cannot use obsCovs, but you can use yearlySiteCovs (a confusing name
    since this model isn't for multi-year data. It's just a hold-over
    from the colext methods upon which it is based.)
    }
\seealso{
    \code{\link{distsamp}}
    }
\examples{

set.seed(36837)

R <- 50
T <- 5
strip.width <- 50
transect.length <- 100
breaks <- seq(0, 50, by=10)

lambda <- 5
phi <- 0.6
sigma <- 30

J <- length(breaks)-1
y <- array(0, c(R, J, T))
for(i in 1:R) {
    M <- rpois(1, lambda) # Individuals within the 1-ha strip
    for(t in 1:T) {
        # Distances from point
        d <- runif(M, 0, strip.width)
        # Detection process
        if(length(d)) {
            cp <- phi*exp(-d^2 / (2 * sigma^2)) # half-normal w/ g(0)<1
            d <- d[rbinom(length(d), 1, cp) == 1]
            y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
y <- matrix(y, nrow=R)


umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
    dist.breaks=breaks, tlength=rep(transect.length, R), numPrimary=T)
summary(umf)

m1 <- gdistsamp(~1, ~1, ~1, umf, output="density")
backTransform(m1, type="lambda")
backTransform(m1, type="phi")
backTransform(m1, type="det")

}
\keyword{ models }


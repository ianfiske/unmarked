\name{occuRN.TTD}
\alias{occuRN.TTD}
\title{Fit Royle-Nichols Time-to-detection Models}

\usage{occuRN.TTD(lambdaformula= ~1, detformula = ~ 1, data, K=25,
    ttdDist = c("exp", "weibull"), starts, method="BFGS", se=TRUE,
    engine = c("C", "R"), ...)}

\arguments{
  \item{lambdaformula}{Right-hand sided formula for the abundance at each site.}
  \item{detformula}{Right-hand sided formula for mean time-to-detection.}
  \item{data}{\code{unmarkedFrameOccuTTD} object that supplies the data
    (see \code{\link{unmarkedFrameOccuTTD}}). Note that only single-season
    models are supported by \code{occuRN.TTD}.}
  \item{K}{The upper summation index used to numerically integrate out the
        latent abundance. This should be set high enough so that it does not
        affect the parameter estimates. Computation time will increase with K.}
  \item{ttdDist}{Distribution to use for time-to-detection; either
    \code{"exp"} for the exponential, or \code{"weibull"} for the Weibull,
    which adds an additional shape parameter \eqn{k}.}
  \item{starts}{optionally, initial values for parameters in the optimization.}
  \item{method}{Optimization method used by \code{\link{optim}}.}
  \item{se}{logical specifying whether or not to compute standard errors.}
  \item{engine}{Either "C" or "R" to use fast C++ code or native R
    code during the optimization.}
  \item{\dots}{Additional arguments to optim, such as lower and upper bounds}
}

\description{Fit Royle-Nichols occupancy models with time-to-detection data.}

\value{unmarkedFitOccuRN.TTD object describing model fit.}

\details{

More later.

}

\references{

More later.

Kery, Marc, and J. Andrew Royle. 2016. \emph{Applied Hierarchical Modeling in
  Ecology}, Volume 1. Academic Press.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuTTD}}}

\keyword{models}

\examples{

\dontrun{

# Simulate data
M = 1000 # Number of sites
nrep <- 3 # Number of visits per site
Tmax = 5 # Max duration of a visit
alpha1 = -1 # Covariate on rate
beta1 = 1 # Covariate on density
mu.lambda = 1 # Rate at alpha1 = 0
mu.dens = 1 # Density at beta1 = 0

covDet <- matrix(rnorm(M*nrep),nrow = M,ncol = nrep) #Detection covariate
covDens <- rnorm(M) #Abundance/density covariate
dens <- exp(log(mu.dens) + beta1 * covDens)
sum(N <- rpois(M, dens)) # Realized density per site
lambda <- exp(log(mu.lambda) + alpha1 * covDet) # per-individual detection rate
ttd <- NULL
for(i in 1:nrep) {
  ttd <- cbind(ttd,rexp(M, N*lambda[,i]))  # Simulate time to first detection per visit
}
ttd[N == 0,] <- 5 # Not observed where N = 0; ttd set to Tmax
ttd[ttd >= Tmax] <- 5 # Crop at Tmax

#Build unmarked frame
umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength=5,
                            siteCovs = data.frame(covDens=covDens),
                            obsCovs = data.frame(covDet=as.vector(t(covDet))))

#Fit model
fit <- occuRN.TTD(~covDens, ~covDet, data=umf, K=max(N)+10)

#Compare to truth
cbind(coef(fit), c(log(mu.dens), beta1, log(mu.lambda), alpha1))

#Predict abundance/density values
head(predict(fit, type='abun'))

}
}

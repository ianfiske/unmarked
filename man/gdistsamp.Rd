\name{gdistsamp}
\alias{gdistsamp}
\title{
Fit the generalized distance sampling model of Royle et al. (2004) that 
allows estimation of the probability an individual is available for detection.
}
\description{
This model can be fit to distance sampling data collected at i = 1, 2, ..., R
plots surveyed on t = 1, 2, ..., T occasions. The detections must be binned 
into distance classes. 
}
\usage{
gdistsamp(lambdaformula, phiformula, detformula, data, keyfun = c("halfnorm", "exp", "hazard", "uniform"), output = c("density", "abund"), unitsOut = c("ha", "kmsq"), mixture = c("P", "NB"), K, starts, method = "BFGS", control = list(), se = TRUE, rel.tol=1e-4)
}
\arguments{
  \item{lambdaformula}{
      A right-hand side formula describing the abundance covariates.
}
  \item{phiformula}{
      A right-hand side formula describing the availability covariates.
}
  \item{detformula}{
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
    You cannot use obsCovs.
    }
\seealso{
    \code{\link{distsamp}}
    }
\examples{

R <- 50
T <- 3
radius <- 50
breaks <- seq(0, 50, by=10)
    
lambda <- 5
phi <- 0.5
sigma <- 20
    
J <- length(breaks)-1
A <- (2*radius)^2 / 10000     # Area (ha) of square containing circle
y <- array(0, c(R, J, T))
for(i in 1:R) {
    M <- rpois(1, lambda * A) # Individuals within the square
    N <- rbinom(T, M, phi)    # Individualss available for detection at time t
        
    for(t in 1:T) {
        # coordinates of each individual
        xy <- cbind(x=runif(N[t], -radius, radius), 
            y=runif(N[t], -radius, radius))
        
        # Distances from point
        d <- apply(xy, 1, function(x) sqrt(x[1]^2 + x[2]^2))
        d <- d[d <= radius]

        # Detection process
        if(length(d)) {
            p <- exp(-d^2 / (2 * sigma^2)) # half-normal
            d <- d[rbinom(length(d), 1, p) == 1]
            y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
y <- matrix(y, nrow=R)


umf <- unmarkedFrameGDS(y = y, survey="point", unitsIn="m", dist.breaks=breaks, 
    numPrimary=3)
summary(umf)
    
m1 <- gdistsamp(~1, ~1, ~1, umf)
backTransform(m1, type="lambda")
backTransform(m1, type="phi")
backTransform(m1, type="det")

}
\keyword{ models }


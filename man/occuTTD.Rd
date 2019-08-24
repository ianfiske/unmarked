\name{occuTTD}
\alias{occuTTD}
\title{Fit Single-Season and Dynamic Time-to-detection Occupancy Models}

\usage{occuTTD(psiformula= ~1, gammaformula =  ~ 1, epsilonformula = ~ 1,
    detformula = ~ 1, data, ttdDist = c("exp", "weibull"), 
    linkPsi = c("logit", "cloglog"), starts, method="BFGS", se=TRUE, 
    engine = c("C", "R"), ...)}

\arguments{
  \item{psiformula}{Right-hand sided formula for the initial probability of 
    occupancy at each site.}
  \item{gammaformula}{Right-hand sided formula for colonization probability.}
  \item{epsilonformula}{Right-hand sided formula for extinction probability.}
  \item{detformula}{Right-hand sided formula for mean time-to-detection.}
  \item{data}{\code{unmarkedFrameOccuTTD} object that supplies the data 
    (see \code{\link{unmarkedFrameOccuTTD}}).}
  \item{ttdDist}{Distribution to use for time-to-detection; either
    \code{"exp"} for the exponential, or \code{"weibull"} for the Weibull,
    which adds an additional shape parameter \eqn{k}.}
  \item{linkPsi}{Link function for the occupancy model. Options are  
    \code{"logit"} for the standard occupancy model or \code{"cloglog"} 
    for the complimentary log-log link, which relates occupancy
    to site-level abundance.}
  \item{starts}{optionally, initial values for parameters in the optimization.}
  \item{method}{Optimization method used by \code{\link{optim}}.}
  \item{se}{logical specifying whether or not to compute standard errors.}
  \item{engine}{Either "C" or "R" to use fast C++ code or native R
    code during the optimization.}
  \item{\dots}{Additional arguments to optim, such as lower and upper bounds}
}

\description{Fit time-to-detection occupancy models of Garrard et al. 
  (2008, 2013), either single-season or dynamic. Time-to-detection can be 
  modeled with either an exponential or Weibull distribution.} 

\value{unmarkedFitOccuTTD object describing model fit.}

\details{More later.}

\references{

Garrard, G.E., Bekessy, S.A., McCarthy, M.A. and Wintle, B.A. 2008. When have 
  we looked hard enough? A novel method for setting minimum survey effort 
  protocols for flora surveys. Austral Ecology 33: 986-998.

Garrard, G.E., McCarthy, M.A., Williams, N.S., Bekessy, S.A. and Wintle, 
  B.A. 2013. A general model of detectability using species traits. Methods in 
  Ecology and Evolution 4: 45-52.

}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuTTD}}}

\keyword{models}

\examples{

#More later

}

\name{colext}
\alias{colext}
\title{Fit the colonization-extinction model.}
\usage{colext(psiformula= ~1, gammaformula =  ~ 1, epsilonformula = ~ 1,
pformula = ~ 1, data, starts, B=0, method="BFGS", control=list(),
    se=TRUE)}
\description{Estimate parameters of the colonization-extinction model, including covariate-dependent rates and detection process.}
\details{This function fits the colonization-extinction model of MacKenzie et al (2003).  The colonization and extinction rates
can be modeled with covariates that vary yearly at each site using a logit link.  These covariates
are supplied by special unmarkedMultFrame \code{yearlySiteCovs} slot.
These parameters are specified using the \code{gammaformula} and \code{epsilonformula}
arguments.  The initial probability of occupancy is modeled by
covariates specified in the \code{psiformula}.

The conditional detection rate can also be modeled as a function of covariates that vary at the secondary sampling
period (ie., repeat visits).  These covariates are specified by the first part of the \code{formula} argument and
the data is supplied via the usual \code{obsCovs} slot.

See \code{\link{unmarkedMultFrame}} and \code{\link{formatMult}} for more information on supplying the data.}
\value{unmarkedFitColExt object describing model fit.}
\keyword{models}
\references{MacKenzie, D.I. et al. (2002) Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. Ecology, 83(8), 2248-2255.

MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.  Amsterdam: Academic Press.  Royle, J. A. and R. Dorazio. (2008).}
\arguments{\item{psiformula}{Right-hand sided formula for the initial
    probability of occupancy at each site.}
  \item{gammaformula}{Right-hand sided formula for colonization probability.}
  \item{epsilonformula}{Right-hand sided formula for extinction probability.}
  \item{pformula}{Right-hand sided formula for detection probability.}
\item{data}{unmarkedMultFrame object that supplies the data (see \code{\link{unmarkedMultFrame}}).}
\item{starts}{optionally, initial values for parameters in the optimization.}
\item{B}{number of bootstrap interations (the default 0 indicates no bootstrapping).}
\item{method}{Optimization method used by \code{\link{optim}}.}
\item{control}{Other arguments passed to \code{\link{optim}}.}
\item{se}{logical specifying whether or not to compute standard errors.}}
\examples{data(frogs)
umf <- formatMult(masspcru)
obsCovs(umf) <- scale(obsCovs(umf))
## constant transition rates
(fm <- colext(psiformula = ~ 1,
gammaformula = ~ 1,
epsilonformula = ~ 1,
pformula = ~ JulianDate + I(JulianDate^2), umf, control = list(trace=1, maxit=1e4)))

\dontrun{
## try yearly transition rates
yearlySiteCovs(umf) <- data.frame(year = factor(rep(1:7, numSites(umf))))
(fm.yearly <- colext(psiformula = ~ 1,
gammaformula = ~ year,
epsilonformula = ~ year,
pformula = ~ JulianDate + I(JulianDate^2), umf, 
	control = list(trace=1, maxit=1e4)))
}
}

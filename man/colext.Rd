\name{colext}
\alias{colext}
\title{Fit the colonization-extinction model.}
\usage{colext(formula=~1 ~ 1, data, starts, B=0, method=BFGS, control=list())}
\description{Estimate parameters of the colonization-extinction model, including covariate-dependent rates and detection process.}
\value{unmarkedFitColExt object describing model fit.}
\keyword{models}
\references{MacKenzie, D.I. et al., 2002. Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. Ecology, 83(8), 2248-2255.}
\arguments{\item{formula}{double right-hand side formula describing covariates of detection and occupancy in that order.}
\item{data}{unmarkedMultFrame object that supplies the data (see \link{unmarkedFrame}).}
\item{starts}{optionally, initial values for parameters in the optimization.}
\item{B}{number of bootstrap interations (the default 0 indicates no bootstrapping).}}
\examples{data(frogs)
umf <- formatMult(masspcru)
obsCovs(umf) <- scale(obsCovs(umf))
fm <- colext(~ JulianDate + I(JulianDate^2) ~ 1, umf, control = list(trace=1, maxit=1e4))
fm}

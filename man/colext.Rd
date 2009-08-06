\name{colext}
\alias{colext}
\title{Fit the general multistate multiseason occupancy model.}
\usage{colext(formula=~1 ~ 1, data, starts, B=0, method=BFGS, control=list())}
\description{Estimate parameters of the colonization-extinction model.}
\details{Each freely varying detection parameter can be modeled as a linear function of observation level covariates.
\code{detconstraint} is a matrix with number of rows equal to the detection matrix parameters and number of
columns equal to the number of covariates (plus 1 for the intercept).  Each column is a constraint vector beginning with 1 that indicates
which covariates are restricted to have the same effect on a given matrix parameter.}
\value{unmarkedFitColExt object describing model fit.}
\keyword{models}
\arguments{\item{formula}{double right-hand side formula describing covariates of detection and occupancy in that order.}
\item{data}{unmarkedMultFrame object that supplies the data (see \link{unmarkedFrame}).}
\item{starts}{optionally, initial values for parameters in the optimization.}
\item{B}{number of bootstrap interations (the default 0 indicates no bootstrapping).}}
\examples{data(frogs)
umf <- formatMult(masspcru)
obsCovs(umf) <- scale(obsCovs(umf))
fm <- colext(~ JulianDate + I(JulianDate^2) ~ 1, umf, control = list(trace=1, maxit=1e4))
fm}

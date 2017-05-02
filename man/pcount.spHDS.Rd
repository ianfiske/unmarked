\name{pcount.spHDS}
\alias{pcount.spHDS}
\encoding{UTF-8}

\title{
Fit spatial hierarchical distance sampling model
}
\description{
Function fits spatial hierarchical distance sampling model; \pkg{unmarked} model fitting function introduced in Section 9.8.4. See \code{\link[unmarked]{pcount}} for details.
}
\usage{
pcount.spHDS(formula, data, K, mixture = c("P", "NB", "ZIP"), starts,
  method = "BFGS", se = TRUE, ...)
}
\arguments{
  \item{formula}{
Double right-hand side formula describing covariates of detection and abundance, in that order

Detection model should be specified without an intercept, for example:
~ -1 + I(dist^2)
Internally this forces the intercept p(0) = 1, conventional for distance
sampling models (see section 9.8.4 for explanation). More general models
work but may not honor that constraint. e.g., 
~ 1,
~ dist,
~ I(dist^2),
~ dist + I(dist^2)
}
  \item{data}{
an \code{unmarkedFramePCount} object supplying data to the model.
}
  \item{K}{
Integer upper index of integration for N-mixture. This should be set high enough so that it does not affect the parameter estimates. Note that computation time will increase with K.
}
  \item{mixture}{
character specifying mixture: Poisson (P), Negative-Binomial (NB), or Zero Inflated Poisson (ZIP).
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
  \item{engine}{
Either "C" or "R" to use fast C++ code or native R code during the optimization.
}
  \item{\dots}{
Additional arguments to optim, such as lower and upper bounds
}
}

\value{
unmarkedFit object describing the model fit.
}
\references{
Kery & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.8.4
}
\author{
Kery & Royle
}


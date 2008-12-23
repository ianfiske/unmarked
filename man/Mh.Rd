\name{Mh}
\alias{Mh}
\title{Fit model Mh.}
\usage{Mh(nx)}
\description{Fits model Mh - closed population model with heterogeneity}
\details{This fits the logit-normal version of Model Mh in which the logit of
individual detection probability is assumed to have a normal distribution with
parameters  mu and sigma. These parameters and the unobserved number of
individuals, n0, are estimated by integrated likelihood.}
\value{optim output for now}
\references{Coull and Agresit (1999) Dorazio and Royle (2003)}
\author{Andy Royle \email{aroyle@usgs.gov}}
\keyword{models}
\note{Assumes balanced data right now.
this will be the engine for a poisson reg model for spec richness}
\arguments{\item{nx}{the detection frequency vector, of length J where J is the number of sample occasions}}
\examples{nx<-c(34, 16, 10, 4, 2, 2,0,0,0,0,0,0,0,0)
Mh(nx)}

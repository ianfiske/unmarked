\name{multinomPois}
\alias{multinomPois}
\title{Multinomial-Poisson Mixtures}
\usage{multinomPois(stateformula, detformula, piFun, umf)}
\description{Fit the multinomial-Poisson abundance mixture model.}
\details{This function takes advantage of the closed form of the integrated
likelihood when a latent Poisson distribution is assumed for abundance
at each site and a multinomial distribution is taken for the observation
state. Many common sampling methods can be framed in this context.  For
example, double-observer point counts, removal sampling, and distance
sampling can all be analyzed with this function by specifying the proper
multinomial cell probablilities.  This is done with by supplying the
appropriate function (piFun) argument.  \link{removalPiFun} and \link{doublePiFun}
are supplied as example cell probability functions.}
\value{unmarkedFit object describing the model fit.}
\author{Ian Fiske}
\keyword{models}
\references{Royle, J., Dawson, D., & Bates, S. (2004). Modeling abundance effects in distance sampling. Ecology, 85(6), 1591-1597.

Royle, J. A. (2004). Generalized estimators of avian abundance from count survey data. Animal Biodiversity and Conservation, 27(1), 375-386.

Royle, J. A., & Dorazio, R. M. (2006). Hierarchical Models of Animal Abundance and Occurrence. Journal Of Agricultural Biological And Environmental Statistics, 11(3), 249.}
\arguments{\item{stateformula}{Right-hand side formula describing covariates of abundance}
\item{detformula}{Right-hand side formula describing covariates of detection}
\item{piFun}{Function to define multinomial cell probabilities.}
\item{umf}{unmarkedFrame supplying data.}}
\examples{data(ovendata)
ovenFrame <- unmarkedFrame(ovendata.list$data,
siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])))
fm1 <- multinomPois(~ ufp + trba, ~ 1, removalPiFun, ovenFrame)
fm2 <- multinomPois(~ ufp, ~1, removalPiFun, ovenFrame)
fm3 <- multinomPois(~ trba, ~1, removalPiFun, ovenFrame)
fm4 <- multinomPois(~ 1, ~1, removalPiFun, ovenFrame)
fm4
fm1}

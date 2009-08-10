\name{pcount}
\alias{pcount}
\title{Fit the N-mixture point count model...}
\usage{pcount(formula, data, K, mixture=c("P", "NB"), starts, method=BFGS,
    control=list(), se=TRUE)}
\description{Fit the N-mixture point count model}
\details{This function fits binomial-Poisson mixture model for spatially replicated point count data.

See \code{\link{unmarkedFrame}} for a description of how to supply by creating
and unmarkedFrame.

This function fits the latent N-mixture model for point count data
(Royle 2004, Kery and Royle 2005).

The latent abundance distribution, \eqn{f(N | \mathbf{\theta})}{f(N |
theta)} can be set as either a Poisson or a negative binomial random
variable, depending on the setting of the \code{mixture} argument.
\code{mixture = "P"} or \code{mixture = "NB"} select the Poisson or
negative binomial distribution respectively.  The mean of \eqn{N_i} is
\eqn{\lambda_i}{lambda_i}.  If \eqn{N_i \sim NB}{N_i ~ NB}, then an
additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
\eqn{\alpha}{alpha} implies higher variance).

The detection process is modeled as binomial: \eqn{y_{ij} \sim
Binomial(N_i, p_{ij})}{y_ij ~ Binomial(N_i, p_ij)}.

Covariates of \eqn{\lambda_i}{lamdba_i} use the log link and
covariates of \eqn{p_{ij}}{p_ij} use the logit link.}
\value{unmarkedFit object describing the model fit.}
\author{Ian Fiske \email{ianfiske@gmail.com}}
\references{Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Kery, M. and Royle, J. A. (2005) Modeling Avaian Abundance from Replicated Counts Using Binomial Mixture Models. \emph{Ecological Applications} 15(4), pp. 1450--1461.}
\keyword{models}
\arguments{\item{stateformula}{Right-hand side formula describing covariates of abundance}
\item{detformula}{Right-hand side formula describing covariates of detection}
\item{umf}{an unmarkedFrame supplying data to the model.}
\item{K}{Integer upper index of integration for N-mixture.}
\item{mixture}{character specifying mixture: either "P" or "NB".}
\item{method}{Optimization method used by \code{\link{optim}}.}
\item{control}{Other arguments passed to \code{\link{optim}}.}
\item{se}{logical specifying whether or not to compute standard errors.}}
\examples{data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
obsCovs = mallard.obs)
(fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest, mallardUMF))
(fm.mallard.nb <- pcount(~ date + I(date^2) ~ length + elev, mixture = "NB", mallardUMF))}

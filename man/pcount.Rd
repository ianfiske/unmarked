\name{pcount}
\alias{pcount}
\title{ Fit N Mixture Point Count Models }
\description{
Fit the latent abundance N-mixture model (Royle 2004, Kery and Royle 2005). 
}
\usage{
pcount(stateformula, detformula, y, sitedata = NULL, obsdata = NULL,
K = NULL, mixture = "P")
}
\arguments{
  \item{stateformula}{ formula for covariates of occurrance. }
  \item{detformula}{ formula for covariates of detection. }
  \item{y}{M x J observation matrix of counts. }
  \item{covdata.site}{ optional data frame containing site-specific covariates. }
  \item{covdata.obs}{ optional array or list containing covariate data. }
  \item{K}{ the upper index of summation over abundance (defaults to
    \eqn{\max(y) + 20}{max(y) + 20)}
  \item{mixture}{ the latent abundance distribution, either "P" for
    Poisson (default), or "NB" for negative binomial.}
}
\details{
  
  See \link{unmarked} for detailed descriptions of passing data \code{y},
  \code{covdata.site}, and \code{covdata.obs}, and specifying covariates
  with \code{stateformula} and \code{detformula}.
  
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
  covariates of \eqn{p_{ij}}{p_ij} use the logit link.

}
\value{
  \item{estimates}{parameter estimates}
  \item{AIC}{AIC computed as \eqn{-2 \times log(like) + 2 \times nP}{-2ll + 2nP}}
}
\references{
Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from
Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Kery, M. and Royle, J. A. (2005) Modeling Avaian Abundance from
Replicated Counts Using Binomial Mixture Models. \emph{Ecological
  Applications} 15(4), pp. 1450--1461.
}
\section{Warning}{
Sometimes the MLE's can be unstable.  This is indicated by a sensitivity
to the upper summation argument, K.  We suggest trying several values of
K, ranging from 20 to 50.  If parameter estimates vary wildly, be wary
of drawing any conclusions and consider examining your data and/or
reducing the model.
}
\author{ Ian Fiske <ianfiske@gmail.com> }
\seealso{ \link{unmarked}  }
\examples{
data(mallard)
fm.nmx1 <- pcount(~ length + elev + forest,
  ~ ivel+ date + I(date^2), mallard.y, mallard.site, mallard.obs)
}
\keyword{ models }
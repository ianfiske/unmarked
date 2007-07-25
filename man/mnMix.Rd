\name{mnMix}
\alias{mnMix}
\title{ Fix Multinomial Mixture Abundance Models}
\description{
Fit latent abundance models to categorical data such as frog
calling data (Royle and Link 2005).  Currently, covariates can only be
supplied to the detection process (as in Royle and Link 2005).
}
\usage{
mnMix(stateformula = ~ 1, detformula = ~ 1, y, constraint = NULL, covdata.site = NULL, covdata.obs = NULL)
}
\arguments{
  \item{stateformula}{ formula for covariates of occurrance. }
  \item{detformula}{ formula for covariates of detection. }
  \item{y}{ M x J observation matrix of 0's and 1's. }
  \item{constraint}{ vector to describe which detection parameters are
    the same (see \emph{Details}}
  \item{sitedata}{ optional data frame containing site-specific covariates. }
  \item{obsdata}{ optional array or list containing covariate data (see \emph{Details}. }
}
\details{

  See \link{unmark} for detailed descriptions of passing data \code{y},
  \code{covdata.site}, and \code{covdata.obs}, and specifying covariates
  with \code{stateformula} and \code{detformula}.
  
  This function fits the model described in Royle and Link (2005).  This
function assumes that a site has a latent categorical state, taking on
values from 0 to K.  Thus, there are K + 1 multinomial states.  Here, we
estimate the probabilities \eqn{\mathbf{\psi}}{\bold{psi}}, of falling
into the K + 1 categories.

Further, the latent states are observed with error, and the detection
process also follows a multinomial distribution according to equation
(3) on page 2508 of Royle and Link (2005).  The user may specify
covariates of the \eqn{p}'s, which are modeled as in equation (4) on
p. 2509.

The \code{constraint} argument specifies which constraints are placed on
the multinomial detection parameters to reduce model complexity.  The
same constraint vector format is used as described on page 2508 of Royle
and Link (2005).  The default is to allow the \eqn{p_k}'s (\eqn{k = 1,
2, \ldots, K}{k = 1, 2, ..., K}) separate parameters, but constrain the
the \eqn{\beta}{beta}'s to be all equal.

}
\value{
  \item{estimates}{parameter estimates}
  \item{AIC}{AIC computed as \eqn{-2 \times log(like) + 2 \times nP}{-2ll + 2nP}}  
}
\references{
Royle, J. A. and Link W. A. (2005) A general class of multinomial
mixture models for anuran calling survey data. \emph{Ecology}, \bold{86}(9),
pp. 2505--2512.
}
\author{ Ian Fiske <ianfiske@gmail.com> }
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(gf)
fm.mmx1 <- mnMix(~ 1, ~ samp1 + samp3 + temp, con=c(1,2,2,3,3,3),
    y=gf.data, covdata.obs = gf.obs)
}
\keyword{ models }
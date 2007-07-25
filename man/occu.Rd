\name{occu}
\alias{occu}
\title{ Fit Occupancy Models }
\description{
  Fits occupancy models like in Chapter 3 of Royle and Dorazio 2008 and MacKenzie et al (2006).  The user can specify covariates of the probability of occurrence at the site level and the probability of detection at either the site or observation level.
}
\usage{
occu(stateformula, detformula, y, covdata.site = NULL,
     covdata.obs = NULL)
}
\arguments{
  \item{stateformula}{ formula for covariates of occurrance. }
  \item{detformula}{ formula for covariates of detection. }
  \item{y}{M x J observation matrix of 0's and 1's.  }
  \item{covdata.site}{ optional data frame containing site-specific covariates. }
  \item{covdata.obs}{ optional array or list containing covariate data. }
}
\details{

  See \link{unmarked} for detailed descriptions of passing data \code{y},
  \code{covdata.site}, and \code{covdata.obs}, and specifying covariates
  with \code{stateformula} and \code{detformula}.

  \command{occu} fits the traditional occupancy model based on the
  binomial mixture models (MacKenzie et al. 2006, Royle and Dorazio
  2008).  The occupancy state process (\eqn{z_i}) of site \eqn{i} is
  modeled as

  \deqn{z_i \sim Bernoulli(\psi_i)}{z_i ~ Bernoulli(psi_i)}

  The observation process is modeled as

  \deqn{y_{ij} \sim \text{Bernoulli}(p_{ij})}{y_ij ~ Bernoulli(p_ij)}

  Covariates of \eqn{\psi_i}{psi_i} and \eqn{p_{ij}}{p_ij} are modelled
  using the logit link.
  
}
\value{
  \item{estimates}{parameter estimates}
  \item{AIC}{AIC computed as \eqn{-2 \times log(like) + 2 \times nP}{-2ll + 2nP}}
}
\references{

  MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and
  Modeling}.  Amsterdam: Academic Press.

  Royle, J. A. and R. Dorazio. (2008) \emph{Book Name}.

}
\author{ Ian Fiske <ianfiske@gmail.com> }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(wt)
fm1 <- occu(~ forest + elev + I(elev^2) + length,
  ~ ivel + date + I(date^2), wt.y, wt.data, wt.obs)
}
\keyword{ models }
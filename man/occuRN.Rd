\name{occuRN}
\alias{occuRN}
\title{ Fit Occupancy Model with Abundance Linkage}
\description{
This function fits the the occupancy model described in Royle and Nichols (2003) and Royle and Dorazio (2008, Chapter 4).  The user can specify covariates for the latent abundance process and the detection process.
}
\usage{
occuRN(stateformula, detformula, data, y, covdata.site = NULL,
       covdata.obs = NULL)
}
\arguments{
  \item{stateformula}{ formula for covariates of occurrance. }
  \item{detformula}{ formula for covariates of detection. }
  \item{data}{ list containing y, covdata.site, and covdata.obs}
  \item{y}{ M x J observation matrix of 0's and 1's. }
  \item{covdata.site}{ optional data frame containing site-specific covariates. }
  \item{covdata.obs}{ optional array or list containing covariate data (see \emph{Details}. }
}
\details{

  See \link{unmarked} for detailed descriptions of passing data \code{y},
  \code{covdata.site}, and \code{covdata.obs}, and specifying covariates
  with \code{stateformula} and \code{detformula}.

  This function fits the latent abundance mixture model described in
  Royle and Nichols (2003).

  The latent abundance of site \eqn{i} is modelled as Poisson:

  \deqn{N_i \sim Poisson(\lambda_i)}{N_i ~ Poisson(lambda_i)}

  The detection of a single individual in site \eqn{i} during sample
  \eqn{j} is modelled as Bernoulli:

  \deqn{w_{ij} \sim Bernoulli(r_{ij})}{w_ij ~ Bernoulli(r_ij)}.

  Thus, the detection probability for a single site is linked to the
  detection probability for an individual by

  \deqn{p_{ij} = 1 - (1 - r_{ij}) ^ {N_i}}{p_ij = 1 - (1 - r_ij) ^ N_i}

  Covariates of \eqn{\lambda_i}{lambda_i} are modelled with the log link
  and covariates of \eqn{r_{ij}}{r_ij} are modelled with the logit link.

}
\value{
  \item{estimates}{parameter estimates}
  \item{AIC}{AIC computed as \eqn{-2 \times log(like) + 2 \times nP}{-2ll + 2nP}}
}
\references{
Royle, J. A. and Nichols, J. D. (2003) Estimating Abundance from
Repeated Presence-Absence Data or Point Counts. \emph{Ecology}, 84(3)
pp. 777--790.
}
\author{ Ian Fiske <ianfiske@gmail.com> }
\seealso{ \link{unmarked} }
\examples{
data(birds)
fm.wood.rn <- occuRN(~1, ~1, y = woodthrush.bin)
}
\keyword{ models }
\name{occuRN}
\alias{occuRN}
\title{Fit the Occupancy model of Royle and Nichols...}
\usage{occuRN(formula, data, K=25)}
\description{Fit the Occupancy model of Royle and Nichols}
\details{See \link{unmarked} for detailed descriptions of passing data \code{y},
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
and covariates of \eqn{r_{ij}}{r_ij} are modelled with the logit link.}
\value{unmarkedFit object describing the model fit.}
\author{Ian Fiske}
\references{Royle, J. A. and Nichols, J. D. (2003) Estimating Abundance from
Repeated Presence-Absence Data or Point Counts. \emph{Ecology}, 84(3)
pp. 777--790.}
\keyword{models}
\arguments{\item{stateformula}{Right-hand side formula describing covariates of abundance}
\item{detformula}{Right-hand side formula describing covariates of detection}
\item{umf}{unmarkedFrame supplying data to the model.}}
\examples{data(birds)
woodthrushUMF <- unmarkedFrame(woodthrush.bin)
fm.wood.rn <- occuRN(~ obs ~ 1, woodthrushUMF)
fm.wood.rn}

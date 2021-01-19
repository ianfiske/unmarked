\name{occuRN}

\alias{occuRN}

\title{Fit the occupancy model of Royle and Nichols (2003)}

\usage{occuRN(formula, data, K=25, starts, method="BFGS", se=TRUE, 
              engine=c("C","R"), threads=1, ...)}

\arguments{
    \item{formula}{double right-hand side formula describing covariates of
        detection and abundance, in that order.}
    \item{data}{Object of class \code{\link{unmarkedFrameOccu}} supplying data
        to the model.}
    \item{K}{the upper summation index used to numerically integrate out the
        latent abundance. This should be set high enough so that it does not
        affect the parameter estimates. Computation time will increase with K.}
    \item{starts}{initial values for the optimization.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{logical specifying whether or not to compute standard
      errors.}
    \item{engine}{Either "C" to use fast C++ code or "R" to use native R
      code during the optimization.}
    \item{threads}{Set the number of threads to use for optimization in C++, if
      OpenMP is available on your system. Increasing the number of threads
      may speed up optimization in some cases by running the likelihood 
      calculation in parallel. If \code{threads=1} (the default), OpenMP is disabled.} 
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{Fit the occupancy model of Royle and Nichols (2003), which relates probability
of detection of the species to the number of individuals available for detection at each site.
Probability of occupancy is a derived parameter: the probability that at least one
individual is available for detection at the site.}

\details{

This function fits the latent abundance mixture model described in
Royle and Nichols (2003).

The number of animals available for detection at site \eqn{i} is modelled as Poisson:

\deqn{N_i \sim Poisson(\lambda_i)}{N_i ~ Poisson(lambda_i)}

We assume that all individuals at site \eqn{i} during sample \eqn{j} have identical
detection probabilities, \eqn{r_{ij}}{r_ij}, and that detections are independent. The species
will be recorded if at least one individual is detected. Thus, the detection
probability for the species is linked to the
detection probability for an individual by

\deqn{p_{ij} = 1 - (1 - r_{ij}) ^ {N_i}}{p_ij = 1 - (1 - r_ij) ^ N_i}

Note that if \eqn{N_i = 0}, then \eqn{p_{ij} = 0}{p_ij = 0}, and increasing values of \eqn{N_i} lead to higher values of \eqn{p_{ij}}{p_ij} The equation for the detection history is then:

\deqn{y_{ij} \sim Bernoulli(p_{ij})}{y_ij ~ Bernoulli(p_ij)}

Covariates of \eqn{\lambda_i}{lambda_i} are modelled with the log link
and covariates of \eqn{r_{ij}}{r_ij} are modelled with the logit link.
}

\value{unmarkedFit object describing the model fit.}

\author{Ian Fiske}

\references{

Royle, J. A. and Nichols, J. D. (2003) Estimating Abundance from
Repeated Presence-Absence Data or Point Counts. \emph{Ecology}, 84(3)
pp. 777--790.

}

\examples{

\dontrun{

data(birds)
woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin)
# survey occasion-specific detection probabilities
(fm.wood.rn <- occuRN(~ obsNum ~ 1, woodthrushUMF))

# Empirical Bayes estimates of abundance at each site
re <- ranef(fm.wood.rn)
plot(re)


}

}

\keyword{models}

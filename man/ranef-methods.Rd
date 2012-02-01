\name{ranef-methods}
\docType{methods}
\alias{ranef}
\alias{ranef-methods}
\alias{ranef,unmarkedFitOccu-method}
\alias{ranef,unmarkedFitPCount-method}
\alias{ranef,unmarkedFitMPois-method}
\alias{ranef,unmarkedFitDS-method}
\title{ Methods for Function \code{ranef} in Package \pkg{unmarked} }
\description{
Compute Best Unbiased Predictor (BUP) of latent variables (abundance or
occurrence) using empirical Bayes methods.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "unmarkedFitOccu")}}{}
\item{\code{signature(object = "unmarkedFitPCount")}}{}
\item{\code{signature(object = "unmarkedFitMPois")}}{}
\item{\code{signature(object = "unmarkedFitDS")}}{}

}
}

\section{Warning}{
  Empirical Bayes methods typically under-estimate the variance of the
  posterior distribution because they do not account for uncertainty in
  the estimates of the fixed effects. See Laird and Louis (1987) for
  discussion and possible remedies.
}

\references{
  Laird, N.M. and T.A. Louis. 1987. Empirical Bayes Confidence Intervals
  Based on Bootstrap Samples. Journal of the American Statistical
  Association 82:739--750.
  }

\examples{
# Simulate data under N-mixture model
set.seed(4564)
R <- 20
J <- 5
N <- rpois(R, 3)
y <- matrix(NA, R, J)
y[] <- rbinom(R*J, N, 0.5)

# Fit model
umf <- unmarkedFramePCount(y=y)
fm <- pcount(~1 ~1, umf, K=15)

# Estimates of abundance at each site
(re <- ranef(fm))
coef(re)
confint(re, level=0.9)
plot(re, xlim=c(0,10))

# Compare estimates to truth
sum(N)
sum(coef(re))
colSums(confint(re))


}
\keyword{methods}


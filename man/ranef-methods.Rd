\name{ranef-methods}
\docType{methods}
\alias{ranef}
\alias{ranef-methods}
\alias{ranef,unmarkedFitOccu-method}
\alias{ranef,unmarkedFitOccuRN-method}
\alias{ranef,unmarkedFitPCount-method}
\alias{ranef,unmarkedFitMPois-method}
\alias{ranef,unmarkedFitDS-method}
\alias{ranef,unmarkedFitGMM-method}
\alias{ranef,unmarkedFitGDS-method}
\alias{ranef,unmarkedFitGMMorGDS-method}
\alias{ranef,unmarkedFitColExt-method}
\alias{ranef,unmarkedFitPCO-method}
\title{ Methods for Function \code{ranef} in Package \pkg{unmarked} }
\description{
Compute conditional posterior distributions of the random variables (latent abundance or
occurrence) using empirical Bayes methods. These methods return an
object storing the posterior distributions of the latent variables at
each site, and for each year (primary period) in the case of open
population models. See \link{unmarkedRanef-class} for methods used to
manipulate the returned object.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "unmarkedFitOccu")}}{Computes the
  conditional distribution of occurrence given the data and the estimates of the
  fixed effects, \eqn{Pr(z_i=1 | y_{ij}, \hat{\psi}_i, \hat{p}_{ij})}{Pr(z(i)=1 |
  y(i,j), psi(i), p(i,j))}}
\item{\code{signature(object = "unmarkedFitOccuRN")}}{Computes the
  conditional abundance distribution given the data and the estimates of the
  fixed effects, \eqn{Pr(N_i=k | y_{ij}, \hat{\psi}_i, \hat{r}_{ij}) for k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,j), psi(i), r(i,j)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitPCount")}}{\eqn{Pr(N_i=k |
  y_{ij}, \hat{\lambda}_i, \hat{p}_{ij})  for k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,j), lambda(i), p(i,j)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitMPois")}}{\eqn{Pr(N_i=k |
  y_{ij}, \hat{\lambda}_i, \hat{p}_{ij}) for k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,j), lambda(i), p(i,j)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitDS")}}{\eqn{Pr(N_i=k |
  y_{i,1:J}, \hat{\lambda}_i, \hat{\sigma}_{i}) for k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,1:J), lambda(i), sigma(i)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitGMM")}}{\eqn{Pr(M_i=k |
  y_{i,1:J,t}, \hat{\lambda}_i, \hat{\phi}_{it}, \hat{p}_{ijt})  for k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,1:J,t), lambda(i), phi(i,t), p(i,j,t)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitGDS")}}{\eqn{Pr(M_i=k |
  y_{i,1:J,t}, \hat{\lambda}_i, \hat{\phi}_{it}, \hat{\sigma}_{it}) for k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,1:J,t), lambda(i), phi(i,t), sigma(i,t)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitColExt")}}{\eqn{Pr(z_{it}=1 |
  y_{ijt}, \hat{\psi}_i, \hat{\gamma}_{it}, \hat{\epsilon}_{it}, \hat{p}_{ijt})}{Pr(z(i,t)=1 |
  y(i,j,t), psi(i), gamma(i,t), epsilon(i,t), p(i,j,t))} for k = 0,1,...,K}
\item{\code{signature(object = "unmarkedFitPCO")}}{\eqn{Pr(N_{it}=k |
  y_{ijt}, \hat{\lambda}_i, \hat{\gamma}_{it}, \hat{\omega}_{it},
  \hat{p}_{ijt}) for k = 0,1,\dots,K}{Pr(N(i,t)=k |
  y(i,j,t), lambda(i), gamma(i,t), omega(i,t), p(i,j,t)) for k = 0,1,...,K}}

}
}


\references{
  Royle, J.A and R.M. Dorazio. 2008. Hierarchical Modeling and Inference
  in Ecology. Academic Press.

  Laird, N.M. and T.A. Louis. 1987. Empirical Bayes confidence intervals
  based on bootstrap samples. Journal of the American Statistical
  Association 82:739--750.
  }

  \seealso{
  \link{unmarkedRanef-class}
  }

\examples{
# Simulate data under N-mixture model
set.seed(4564)
R <- 20
J <- 5
N <- rpois(R, 10)
y <- matrix(NA, R, J)
y[] <- rbinom(R*J, N, 0.5)

# Fit model
umf <- unmarkedFramePCount(y=y)
fm <- pcount(~1 ~1, umf, K=50)

# Estimates of conditional abundance distribution at each site
(re <- ranef(fm))
postMode(re)           # EBUP
confint(re, level=0.9) # 90\% CI

# Plots
plot(re, subset=site \%in\% c(1:10), layout=c(5, 2), xlim=c(-1,20))

# Compare estimates to truth
sum(N)
sum(postMode(re))
colSums(confint(re))

# Extract all values in convenient formats
post.df <- as(re, "data.frame")
head(post.df)
post.arr <- as(re, "array")


}
\keyword{methods}


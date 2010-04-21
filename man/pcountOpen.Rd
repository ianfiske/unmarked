\name{pcountOpen}
\alias{pcountOpen}
\title{Fit the open N-mixture model of Dail and Madsen (2010)}
\description{Fit the open N-mixture model of Dail and Madsen (2010)}
\usage{
pcountOpen(lambdaformula, gammaformula, omegaformula, pformula, data, 
	mixture = c("P", "NB"), K, starts, method = "BFGS", se = TRUE, ...)
}
\arguments{
  \item{lambdaformula}{Right-hand sided formula for initial abundance}
  \item{gammaformula}{Right-hand sided formula for recruitment rate}
  \item{omegaformula}{Right-hand sided formula for survival probability}
  \item{pformula}{Right-hand sided formula for detection probability}
  \item{data}{An object of class unmarkedFramePCountOpen. See details}
  \item{mixture}{character specifying mixture: either "P" or "NB".}
  \item{K}{Integer upper index of integration for N-mixture.}
  \item{starts}{vector of starting values}
  \item{method}{Optimization method used by \code{\link{optim}}.}
  \item{se}{logical specifying whether or not to compute standard errors.}
  \item{\dots}{additional arguments to be passed to \code{\link{optim}}.}
}
\details{
This model generalizes the Royle (2004) N-mixture model by relaxing the closure assumption.

The latent initial abundance distribution, 
\eqn{f(N | \mathbf{\theta})}{f(N | theta)} can be set as either a 
Poisson or a negative binomial random variable, depending on the setting of the 
\code{mixture} argument. \code{mixture = "P"} or \code{mixture = "NB"} select 
the Poisson or negative binomial distribution respectively.  The mean of 
\eqn{N_i} is \eqn{\lambda_i}{lambda_i}.  If \eqn{N_i \sim NB}{N_i ~ NB}, then an
additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
\eqn{\alpha}{alpha} implies higher variance).

The latent abundance state following the initial sampling period arises from a
Markovian process in which survivors are modeled as \eqn{S_it \sim 
Binomial(N_it-1, omega_it)}{G_it ~ Binomial(N_it-1, omega_it)}, and recruits
follow \eqn{G_it \sim Poisson(gamma_it)}{G_it ~ Poisson(gamma_it)}.

The detection process is modeled as binomial: \eqn{y_{it} \sim
Binomial(N_it, p_it)}{y_it ~ Binomial(N_it, p_it)}.

Covariates of \eqn{\lambda_i}{lamdba_i} use the log link and
covariates of \eqn{p_it}{p_it} use the logit link.
 
}
\value{
A unmarkedFitPCountOpen object.
}
\references{
Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.
}

\author{Richard Chandler \email{rchandler@nrc.umass.edu}}

\note{
When gamma or omega are modeled using observation-level covariates, the covariate data for the final survey occasion will be ignored; however, it must be supplied.

If time intervals vary among survey occasions (sampling periods), an M by T-1 matrix of time intervals should be supplied to unmarkedFramePCountOpen using the delta argument.
}

\section{Warning}{This function can be extremely slow, especially when there are covariates on gamma or omega. Consider testing the timing on a small subset of the data, perhaps with se=FALSE. Finding the lowest value of K that does not affect estimates (MLEs) will also help with speed. }

\seealso{
\code{\link{pcount}}
}

\examples{
# Simulate	
set.seed(3)
M <- 50
T <- 5
veght <- rnorm(M)
isolation <- matrix(rnorm(M*T), M, T)
time <- matrix(rnorm(M*T, 1), M, T)
lambda <- exp(-1)   # + 0.5*veght)
y <- p <- N <- gamma <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
gamma[] <- exp(-1 + -1*isolation)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], 0.8)
	G[,t] <- rpois(M, gamma[,t])
	N[,t+1] <- S[,t] + G[,t]
	}
p[] <- plogis(-1 + 1*time)
y[] <- rbinom(M*T, N, p)

                            
# Prepare data                               
umf <- unmarkedFramePCountOpen(y = y, siteCovs = data.frame(veght), 
	obsCovs = list(isolation=isolation, time=time))

head(umf)
summary(umf)
plot(umf, strip=FALSE)


# Fit some models
(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=15))   # K may be too small
backTransform(m1, "lambda")
backTransform(m1, "gamma")
backTransform(m1, "omega")
backTransform(m1, "det")

\dontrun{
(m2 <- pcountOpen(~1, ~isolation, ~1, ~time, umf, K=15, se=FALSE,
	control=list(maxit=30, trace=TRUE), starts=c(-1, -1, -1, 1.5, -1, 1)))

# Now simulate with uneven sampling period intervals
set.seed(3)
M <- 50
T <- 10
date <- matrix(NA, M, T)
datediff <- date[,-T]
date[,1] <- rpois(M, 2)
for(t in 2:T) {
    date[,t] <- date[,t-1] + rpois(M, 10)
    datediff[,t-1] <- date[,t] - date[,t-1]
    }
lambda <- exp(2)
y <- p <- N <- gamma <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], 0.95^datediff[,t])   # 0.95 is daily survival rate
	G[,t] <- rpois(M, 0.05*datediff[,t])           # 0.05 is daily recruitment
	N[,t+1] <- S[,t] + G[,t]
	}
p[] <- 0.7
y[] <- rbinom(M*T, N, p)


umfO <- unmarkedFramePCountOpen(y = y, delta = datediff)
umfO

(m3 <- pcountOpen(~1, ~1, ~1, ~1, umfO, se=FALSE, starts=c(2, -3, 3, 0.85),
    control=list(maxit=15, trace=T, REPORT=1)))
backTransform(m3, "lambda")
backTransform(m3, "gamma")
backTransform(m3, "omega")
backTransform(m3, "det")

       
}

}
\keyword{models}

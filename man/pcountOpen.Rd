\name{pcountOpen}
\alias{pcountOpen}
\title{Fit the open N-mixture model of Dail and Madsen (2010)}
\description{}
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
This model generalizes the Royle (2004) N-mixture model, by allowing N to vary among sites and survey occasions. By relying on a Markovian assumption, parameters are identifiable without secondary sampling periods. 
}
\value{
A unmarkedFitPCountOpen object.
}
\references{
%% ~put references to the literature/web site here ~
}

\author{Richard Chandler \email{rchandler@nrc.umass.edu}}

\note{
When gamma or omega are modeled using observation-level covariates, the covariate data for the final survey occasion will be ignored. However, it must be supplied because of the way that unmarkedFrame stores data.
}

\section{Warning }{This function can be extremely slow, especially when there are covariates on gamma or omega. Consider testing the timing on a small subset of the data.}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
# Simulate	
set.seed(1)
M <- 50
T <- 5
veght <- rnorm(M)
isolation <- matrix(rnorm(M*T), M, T)
date <- matrix(rnorm(M*T, 1), M, T)
lambda <- exp(-1 + 0.5*veght)
y <- p <- N <- gamma <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
gamma[] <- exp(-1 + -1*isolation)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], 0.8)
	G[,t] <- rpois(M, gamma[,t])
	N[,t+1] <- S[,t] + G[,t]
	}
p[] <- plogis(-1 + 1*date)
y[] <- rbinom(M*T, N, p)
#y[1, 1:2] <- NA
#isolation[1, 1:2] <- NA
#date[1, 1:2] <- NA
dat <- data.frame(veght)



# Prepare data                               
umf <- unmarkedFramePCountOpen(y = y, siteCovs = dat, 
	obsCovs = list(isolation=isolation, date=date))

head(umf)
summary(umf)
plot(umf)


# Fit some models
system.time(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10))
backTransform(m1, "lambda")
backTransform(m1, "gamma")
backTransform(m1, "omega")
backTransform(m1, "det")

\dontrun{
# Real slow
system.time(m2 <- pcountOpen(~veght, ~isolation, ~1, ~date, umf, K=10, se=F,
	control=list(maxit=20, trace=T), starts=c(-1, 0.5, -1, -1, 1.5, -1, 1)))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

\name{makePiFuns}
\alias{makePiFuns}
\alias{makeRemPiFun}
\alias{makeCrPiFun}
\alias{makeCrPiFunMb}
\alias{makeCrPiFunMh}

\title{Create functions to compute multinomial cell probabilities}

\description{
These are factory functions that generate piFuns with the required defaults, which are enclosed within the environment of the piFun. See the main entry for \link{piFuns}.
}

\usage{
makeRemPiFun(times)
makeCrPiFun(nOcc)
makeCrPiFunMb(nOcc)
makeCrPiFunMh(nOcc)
}

\arguments{
\item{times}{
a vector of times for each interval, \code{length(times)} is the number of survey occasions; can be all 1's if times are the same.
}
\item{nOcc}{
the number of survey occasions
}
}

\details{
\code{makeRemPiFun} produces a piFun for a removal model with the required number of occasions and potentially varying time intervals. The input to the piFun must be probabilities \emph{per unit time}. This is a generalisation of the piFun in the Examples section of \link{piFuns}.

\code{makeCrPiFun} produces a piFun for a standard capture-recapture model, M0, Mt or Mx. Probabilities of detection may vary across occasions. See Kery & Royle (2016) section 7.8.1.

\code{makeCrPiFunMb} produces a piFun for a capture-recapture model with a behavioral response after the first capture, Mb. Probabilities of detection are constant across occasions. The first column is the probability of detection for animals not caught before, column #2 is for animals after the first capture. The remaining columns are ignored. See Kery & Royle (2016) section 7.8.2.

\code{makeCrPiFunMh} produces a piFun for a capture-recapture model with individual heterogeneity in detection probability, Mh, using a logit-normal distribution. Probabilities of detection are constant across occasions. The first column is the mean of the logit-normal on the probability scale. Cell p[1, 2] is a value in [0, 1] which controls the spread of the distribution. The remaining cells are ignored. See Kery & Royle (2016) section 7.8.3.
}

\value{
A piFun with the appropriate defaults.
}

\references{
Kery, M., Royle, J. A. (2016) \emph{Applied Hierarchical Modeling in Ecology} Vol 1.
}

\examples{
# Generate piFuns and check their behaviour:

# makeRemPiFun
# ============
( pRem <- matrix(0.4, nrow=5, ncol=3) )
myPi <- makeRemPiFun(times=c(2,3,5))
myPi(pRem)
ls(environment(myPi))  # See what's in the environment
environment(myPi)$times

( pRem <- matrix(runif(15), 5, 3) )
myPi(pRem)

myPi <- makeRemPiFun(c(5,3,2))
environment(myPi)$times
myPi(pRem)

# More than 3 occasions
myPi <- makeRemPiFun(c(1,2,3,5))
try(myPi(pRem))  # Error
( pRem <- matrix(runif(20), 5, 4) )
myPi(pRem)
# Probability of escaping detection
1 - rowSums(myPi(pRem))

# makeCrPiFun
# ===========
p <- matrix(0.4, 2, 3)
myPi <- makeCrPiFun(3)
myPi(p)
myPi  # Look at the function
ls(environment(myPi))
environment(myPi)$histories

p <- matrix(runif(6, 0.1, 0.9), 2, 3)  # different p's everywhere
myPi(p)

p <- matrix(runif(4*5, 0.1, 0.9), 4, 5)  # > 3 occasions
try(myPi(p))  # Error
myPi <- makeCrPiFun(5)
( tmp <- myPi(p) )
1 - rowSums(tmp) # Probability of non-capture

# makeCrPiFunMb
# ==============
( pMb <- cbind(rep(0.7, 5), 0.3, NA) )
myPi <- makeCrPiFunMb(3)
myPi(pMb)

( pMb <- matrix(runif(15), 5, 3) )  # col #3 will be ignored
myPi(pMb)

# with > 3 occasions
( pMb <- matrix(runif(15), 3, 5) )
try(myPi(pMb))
myPi <- makeCrPiFunMb(5)
myPi(pMb)

# makeCrPiFunMh
# =============
pMh <- cbind(rep(0.4, 5), NA, NA)
pMh[1, 2] <- 0.3
pMh
myPi <- makeCrPiFunMh(3)
myPi(pMh)
pMh <- cbind(runif(5), NA, NA)
pMh[1, 2] <- 0.3
pMh
myPi(pMh)

# with > 3 occasions
pMh <- cbind(runif(5), NA, NA, NA, NA)
pMh[1, 2] <- 0.3
pMh
try(myPi(pMh))
myPi <- makeCrPiFunMh(5)
1 - rowSums(myPi(pMh))  # Probability of non-detection
}

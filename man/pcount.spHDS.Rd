\name{pcount.spHDS}
\alias{pcount.spHDS}
\encoding{UTF-8}

\title{
Fit spatial hierarchical distance sampling model.
}
\description{
Function fits an N-mixture model for a discrete state space with raster covariates, and a detection function which decreases with distance from the observer, assumed to be at the centre. See Kery & Royle (2016) Section 9.8.4 for details.
}
\usage{
pcount.spHDS(formula, data, K, mixture = c("P", "NB", "ZIP"), starts,
  method = "BFGS", se = TRUE, ...)
}
\arguments{
  \item{formula}{
Double right-hand side formula describing covariates of detection and abundance, in that order.

Detection model should be specified without an intercept, for example:
\code{~ -1 + I(dist^2)}, where \code{dist} is a covariate giving the distance of each cell of the raster from the observer.
Internally this forces the intercept \code{p(0) = 1}, conventional for distance
sampling models (see Kery & Royle (2016) for explanation). More general models
work but may not honor that constraint. e.g.,
\code{~ 1,
~ dist,
~ I(dist^2),
~ dist + I(dist^2)}
}
  \item{data}{
an \code{unmarkedFramePCount} object supplying data to the model.
}
  \item{K}{
Integer upper index of integration for N-mixture. This should be set high enough so that it does not affect the parameter estimates. Note that computation time will increase with K.
}
  \item{mixture}{
character specifying mixture: Poisson (P), Negative-Binomial (NB), or Zero Inflated Poisson (ZIP).
}
  \item{starts}{
vector of starting values
}
  \item{method}{
Optimization method used by \code{\link{optim}}.
}
  \item{se}{
logical specifying whether or not to compute standard errors.
}
  \item{\dots}{
Additional arguments to optim, such as lower and upper bounds
}
}

\value{
unmarkedFit object describing the model fit.
}
\references{
Kery & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.8.4
}
\author{
Kery & Royle
}
\examples{
## Simulate some data to analyse
# This is based on Kery and Royle (2016) section 9.8.3
# See AHMbook::sim.spatialDS for more simulation options.

# We will simulate distance data for a logit detection function with sigma = 1,
# for a 6x6 square, divided into a 30 x 30 grid of pixels (900 in all), with the
# observer in the centre.

set.seed(2017)

## 1. Create coordinates for 30 x 30 grid
grx <- seq(0.1, 5.9, 0.2)    # mid-point coordinates
gr <- expand.grid(grx, grx)  # data frame with coordinates of pixel centres

## 2a. Simulate spatially correlated Habitat covariate
# Get the pair-wise distances between pixel centres
tmp <- as.matrix(dist(gr))  # a 900 x 900 matrix
# Correlation is a negative exponential function of distance, with scale parameter = 1
V <- exp(-tmp/1)
Habitat <- crossprod(t(chol(V)), rnorm(900))

## 2b. Do a detection covariate: the distance of each pixel centre from the observer
dist <- sqrt((gr[,1]-3)^2 + (gr[,2]-3)^2)

## 3. Simulate the true population
# Probability that an animal is in a pixel depends on the Habitat covariate, with
#   coefficient beta:
beta <- 1
probs <- exp(beta*Habitat) / sum(exp(beta*Habitat))
# Allocate 600 animals to the 900 pixels, get the pixel ID for each animal
pixel.id <- sample(1:900, 600, replace=TRUE, prob=probs)

## 4. Simulate the detection process
# Get the distance of each animal from the observer
# (As an approximation, we'll treat animals as if they are at the pixel centre.)
d <- dist[pixel.id]
# Calculate probability of detection with logit detection function with
sigma <- 1
p <- 2*plogis(-d^2/(2*sigma^2))
# Simulate the 1/0 detection/nondetection vector
y <- rbinom(600, 1, p)
# Check the number of animals detected
sum(y)
# Select the pixel IDs for the animals detected and count the number in each pixel
detected.pixel.id <- pixel.id[y == 1]
pixel.count <- tabulate(detected.pixel.id, nbins=900)

## 5. Prepare the data for unmarked
# Centre the Habitat covariate
Habitat <- Habitat - mean(Habitat)
# Construct the unmarkedFramePCount object
umf <- unmarkedFramePCount(y=cbind(pixel.count),     # y needs to be a 1-column matrix
   siteCovs=data.frame(dist=dist, Habitat=Habitat))
summary(umf)

## 6. Fit some models
(fm0 <- pcount.spHDS(~ -1 + I(dist^2) ~ 1, umf, K = 20))
(fm1 <- pcount.spHDS(~ -1 + I(dist^2) ~ Habitat, umf, K = 20))
# The true Habitat coefficient (beta above) = 1
# fm1 has much lower AIC; look at the population estimate
sum(predict(fm1, type="state")[, 1])
}



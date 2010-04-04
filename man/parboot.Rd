\name{parboot}
\alias{parboot}
\alias{plot,parboot,missing-method}
\alias{show,parboot-method}
\title{Parametric bootstrap method for fitted models inheriting class.}
\description{Simulate datasets from a fitted model, refit the model, and 
  generate a sampling distribution for a user-specified fit-statistic.}
\arguments{
  \item{object}{a fitted model inheriting class "unmarkedFit"}
  \item{statistic}{a function returning a vector of fit-statistics. 
    First argument must be the fitted model. 
    Default is sum of squared residuals.}
  \item{nsim}{number of bootstrap replicates}
  \item{report}{print fit statistic every 'report' iterations during resampling}
  \item{...}{Additional arguments to be passed to statistic}}
\details{This function simulates datasets based upon a fitted model, 
  refits the model, and evaluates a user-specified fit-statistic for each 
  simulation. Comparing this sampling distribution to the observed statistic
  provides a means of evaluating goodness-of-fit or assessing uncertainty in
  a quantity of interest.}
\value{
  An object of class parboot with three slots: 
  \item{call}{parboot call}
  \item{t0}{Numeric vector of statistics for original fitted model.} 
  \item{t.star}{nsim by length(t0) matrix of statistics for each simulation fit.}}
\author{Richard Chandler \email{rchandler@nrc.umass.edu}}
\examples{

data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 
lengths <- linetran$Length

ltUMF <- with(linetran, {
	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
	siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
	tlength = lengths*1000, survey = "line", unitsIn = "m")
    })

# Fit a model
(fm <- distsamp(~ area ~habitat, ltUMF))

# Function returning two fit-stats: sum of squared errors and mean density at
# sampled sites.
fitStats <- function(fit) {
    sse <- SSE(fit)
    D <- mean(exp(coef(fit, type="state")), na.rm=TRUE)
    return(c(sse, D.hat=D))
    }

(pb <- parboot(fm, fitStats, nsim=25))
plot(pb, main="")

}

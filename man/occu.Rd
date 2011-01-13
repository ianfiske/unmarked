\name{occu}
\alias{occu}
\title{Fit the MacKenzie Occupancy Model}
\usage{occu(formula, data, knownOcc=numeric(0), starts, method="BFGS",
    control=list(), se=TRUE)}
\description{This function estimates the standard occupancy model of MacKenzie et al (2002).}
\details{See \link{unmarkedFrame} for a description of how to supply data to the \code{umf}
argument.
\code{occu} fits the standard occupancy model based on zero-inflated
binomial models (MacKenzie et al. 2006, Royle and Dorazio
2008).  The occupancy state process (\eqn{z_i}) of site \eqn{i} is
modeled as

\deqn{z_i \sim Bernoulli(\psi_i)}{z_i ~ Bernoulli(psi_i)}

The observation process is modeled as

\deqn{y_{ij}|z_i \sim Bernoulli(z_i p_{ij})}{y_ij | z_i ~ Bernoulli(z_i  * p_ij)}

Covariates of \eqn{\psi_i}{psi_i} and \eqn{p_{ij}}{p_ij} are modeled
using the logit link according to the \code{formula} argument.  The formula is a double right-hand sided formula
like \code{~ detform ~ occform} where \code{detform} is a formula for the detection process and \code{occform} is a 
formula for the partially observed occupancy state.  See \link{formula} for details on constructing model formulae
in \R.}
\value{unmarkedFitOccu object describing the model fit.}
\references{

MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege,
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy Rates
  When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.
  
MacKenzie, D. I. et al. 2006. \emph{Occupancy Estimation and Modeling}.  
  Amsterdam: Academic Press.  
  
Royle, J. A. and R. Dorazio. 2008. \emph{Hierarchical Modeling and Inference
  in Ecology}. Academic Press.

}
\author{Ian Fiske}
\keyword{models}
\arguments{\item{formula}{double right-hand side formula describing covariates of detection and occupancy in that order.}
\item{data}{an unmarkedFrameOccu object (see \link{unmarkedFrame})..}
\item{knownOcc}{vector of sites that are known to be occupied.}
\item{starts}{vector of parameter starting values.}
\item{method}{Optimization method used by \code{\link{optim}}.}
\item{control}{Other arguments passed to \code{\link{optim}}.}
\item{se}{logical specifying whether or not to compute standard errors.}}
\examples{

# Fake data
R <- 4 # number of sites
J <- 3 # number of visits
y <- matrix(c(
   1,1,0,
   0,0,0,
   1,1,1,
   1,0,1), nrow=R, ncol=J, byrow=TRUE)
y

site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
site.covs

obs.covs <- list(
   x3 = matrix(c(
      -1,0,1,
      -2,0,0,
      -3,1,0,
      0,0,0), nrow=R, ncol=J, byrow=TRUE),
   x4 = matrix(c(
      'a','b','c',
      'd','b','a',
      'a','a','c',
      'a','b','a'), nrow=R, ncol=J, byrow=TRUE))
obs.covs

umf <- unmarkedFrameOccu(y=y, siteCovs=site.covs, 
    obsCovs=obs.covs)   # organize data
umf                     # look at data
summary(umf)            # summarize      
fm <- occu(~1 ~1, umf)  # fit a model


# Real data
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
plot(pferUMF, panels=4)
# add some fake covariates for illustration
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))

# observation covariates are in site-major, observation-minor order
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))

(fm <- occu(~ obsvar1 ~ 1, pferUMF))

confint(fm, type='det', method = 'normal')
confint(fm, type='det', method = 'profile')

# estimate detection effect at obsvars=0.5
(lc <- linearComb(fm['det'],c(1,0.5)))

# transform this to probability (0 to 1) scale and get confidence limits
(btlc <- backTransform(lc))
confint(btlc, level = 0.9)}

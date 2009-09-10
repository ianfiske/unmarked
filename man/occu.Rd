\name{occu}
\alias{occu}
\title{Fit the MacKenzie Occupancy Model}
\usage{occu(formula, data, knownOcc=numeric(0), starts, method=BFGS,
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

\deqn{y_{ij} \sim Bernoulli(p_{ij})}{y_ij ~ Bernoulli(p_ij)}

Covariates of \eqn{\psi_i}{psi_i} and \eqn{p_{ij}}{p_ij} are modeled
using the logit link according to the \code{formula} argument.  The formula is a double right-hand sided formula
like \code{~ detform ~ occform} where \code{detform} is a formula for the detection process and \code{occform} is a 
formula for the partially observed occupancy state.  See \link{formula} for details on constructing model formulae
in \R.}
\value{unmarkedFitOccu object describing the model fit.}
\references{MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege, J. Andrew Royle, and C. A. Langtimm. Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. Ecology 83, no. 8 (2002): 2248-2255. \cr
MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.  Amsterdam: Academic Press.  Royle, J. A. and R. Dorazio. (2008).}
\author{Ian Fiske}
\keyword{models}
\arguments{\item{formula}{double right-hand side formula describing covariates of detection and occupancy in that order.}
\item{data}{an unmarkedFrameOccu object (see \link{unmarkedFrame})..}
\item{knownOcc}{vector of sites that are known to be occupied.}
\item{starts}{vector of parameter starting values.}
\item{method}{Optimization method used by \code{\link{optim}}.}
\item{control}{Other arguments passed to \code{\link{optim}}.}
\item{se}{logical specifying whether or not to compute standard errors.}}
\examples{data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)

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

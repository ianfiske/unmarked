\name{occu}
\alias{occu}
\title{Fit the MacKenzie Occupancy Model}
\usage{occu(stateformula, detformula, umf, knownOcc=numeric(0), profile=FALSE)}
\description{This function estimates the standard occupancy model.}
\details{See \link{unMarkedFrame} for a description of how to supply data to the \command{umf}
argument.

\command{occu} fits the traditional occupancy model based on the
binomial mixture models (MacKenzie et al. 2006, Royle and Dorazio
2008).  The occupancy state process (\eqn{z_i}) of site \eqn{i} is
modeled as

\deqn{z_i \sim Bernoulli(\psi_i)}{z_i ~ Bernoulli(psi_i)}

The observation process is modeled as

\deqn{y_{ij} \sim Bernoulli(p_{ij})}{y_ij ~ Bernoulli(p_ij)}

Covariates of \eqn{\psi_i}{psi_i} and \eqn{p_{ij}}{p_ij} are modelled
using the logit link.}
\value{unMarkedFit object describing the model fit.}
\references{MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege, J. Andrew Royle, and C. A. Langtimm. “Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. Ecology 83, no. 8 (2002): 2248-2255.
MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.  Amsterdam: Academic Press.  Royle, J. A. and R. Dorazio. (2008) \emph{Book Name}.}
\author{Ian Fiske}
\keyword{models}
\arguments{\item{stateformula}{right-hand side formula describing covariates of occurence.}
\item{detformula}{right-hand side formula describing covariates of detection.}
\item{umf}{unMarkedFrame object that supplies the data (see \link{unMarkedFrame})..}
\item{knownOcc}{vector of sites that are known to be occupied.}}
\examples{data(frogs)
pferUMF <- unMarkedFrame(pfer.bin)
fm <- occu(~1, ~1, pferUMF)
fm}

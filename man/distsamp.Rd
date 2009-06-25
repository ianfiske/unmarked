\name{distsamp}
\alias{distsamp}
\title{Fit the hierarchical distance sampling model...}
\usage{distsamp(formula, data, dist.breaks, tlength, keyfun=c("halfnorm",
    "exp", "hazard", "uniform"), survey, output=c("density",
    "abund"), unitsIn, unitsOut=c("ha", "kmsq"), starts,
    method=BFGS, control=list(), ...)}
\description{Fit the hierarchical distance sampling model}
\details{This functions fits the Royle et al. (2004) multinomial mixture model to line 
or point transect data recorded in discrete distance intervals.

Unlike conventional distance sampling, which uses the 'conditional on 
detection' likelihood formulation, this model is based upon the unconditional 
likelihood and thus allows for modeling both abundance and detection 
function parameters. 

The latent transect-level abundance distribution 
\eqn{f(N | \mathbf{\theta})}{f(N | theta)} is currently assumed to be Poisson 
with mean \eqn{\lambda}{lambda}.

The detection process is modeled as multinomial: 
\eqn{y_{ij} \sim Multinomial(N_i, pi_{ij})}{y_ij ~ Multinomial(N_i, pi_i1, pi_i2, ..., pi_iJ)}, 
where \eqn{pi_ij} is the multinomial cell probability for transect i in 
distance class j. These are computed based upon a detection function 
\eqn{g(x | \mathbf{\sigma})}{g(x | sigma)}, such as the half-normal, 
negative exponential, or hazard rate.  

Parameters \eqn{\lambda}{lambda} and \eqn{\sigma}{sigma} can be vectors 
affected by transect-specific covariates using the log link.}
\value{umDistsampFit object (child class of \link{unmarkedFit)} describing the 
model fit. Parameter estimates are displayed on the log-scale. 
Back-transformation can be achieved via the \link{predict} or 
\link{backTransform} methods.}
\note{The response matrix contains the counts of objects at each transect in each 
distance interval. These distance intervals must correspond to the distance 
break points vector dist.breaks. One must not change dist.breaks without also 
changing the response matrix (and vice versa). 

Currently, transect-level abundance is assumed to be Poisson distributed 
though other distributions such as the negative binomial may be added. 
Goodness-of-fit can be assessed using the \code{\link{parboot}} function.}
\author{Richard Chandler \email{rchandler@nrc.umass.edu}}
\references{Royle, J. A., D. K. Dawson, and S. Bates (2004) Modeling 
abundance effects in distance sampling. \emph{Ecology} 85, pp. 1591-1597.}
\seealso{\code{\link{umDistsampFit}} \code{\link{unmarkedFitList}}}
\keyword{models}
\arguments{\item{stateformula}{Right-hand or 2-sided side formula describing covariates 
of abundance or density. See examples for how to specify the response matrix 
in the stateformula when data is a data.frame.}
\item{detformula}{Right-hand side formula describing covariates of detection.}
\item{data}{data.frame or unmarkedFrame containing response and predictor 
variables. The dimensions of the response matrix must be M x J, where M is 
the number of transects and J is the number of distance intervals in which 
observation were recorded.}
\item{dist.breaks}{Numeric vector defining the distance intervals 
(in meters or km) in which observations were recorded.}
\item{tlength}{Numeric vector of transect lengths in meters or km.}
\item{keyfun}{One of the following detection functions: 
"halfnorm", "hazard", "exp", or "uniform." See details.}
\item{survey}{Either "line" or "point" transect.}
\item{output}{Model either "density" or "abund"}
\item{unitsIn}{Either "m" or "km" for BOTH dist.breaks and tlength}
\item{unitsOut}{Units of density. Either "ha" or "kmsq" for hectares and 
square kilometers, respectively.}
\item{starts}{Vector of starting values for parameters.}
\item{method}{Optimization method used by \code{\link{optim}}.}
\item{control}{Other arguments passed to \code{\link{optim}}.}}
\examples{## Line transect examples

# Distance cut points in meters
data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 

lengths <- linetran$Length

ltUMF <- with(linetran, {
unmarkedFrame(y = cbind(dc1, dc2, dc3, dc4), 
siteCovs = data.frame(Length, area, habitat))
})

# Half-normal detection function. Density output. No covariates. 
# lineDat$Length is transect lengths in km, so it has to be converted.
(fm1 <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
tlength=lengths*1000, survey="line", unitsIn="m"))

coef(fm1, type="det", altNames=TRUE)
backTransform(fm1, whichEstimate="det")
vcov(fm1, altNames=TRUE)
predict(fm1, type="state")

# Half-normal. Abundance output. No covariates. Note that transect length
# must be accounted for so abundance is animals per km of transect.
(fm2 <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
tlength=lengths*1000, output="abund", survey="line", unitsIn="m", 
unitsOut="kmsq"))

# Halfnormal. Covariates affecting both density and and detection.  
(fm3.1 <- distsamp(~ poly(area, 2) + habitat ~ habitat, 
ltUMF, dist.breaks=dbreaksLine, tlength=lengths*1000, 
survey="line", unitsIn="m", unitsOut="ha"))

# This won't run without starting values.
(fm3.2 <- distsamp(~ poly(area, 2) + habitat - 1 ~ habitat - 1, ltUMF, 
dist.breaks=dbreaksLine, tlength=lengths*1000, survey="line", 
unitsIn="m", unitsOut="ha", starts=c(1,0,0,1,2,2)))

# Negative exponential detection function. Density output in hectares. 
(fme <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
tlength=lengths*1000, key="exp", survey="line", unitsIn="m"))

# Hazard-rate detection function. Density output in hectares.
# This is real slow, especially for large datasets.
(fmhz <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
tlength=lengths*1000, keyfun="hazard", survey="line", unitsIn="m"))

plot(function(x) unmarked:::gxhaz(x, shape=8.38, scale=1.37), 0, 25, 
xlab="Distance(m)", ylab="Detection probability")

# Uniform detection function. Density output in hectars.
(fmu <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
tlength=lengths*1000, key="uniform", survey="line", unitsIn="m"))

## Point transect examples

# Radius cut points in meters
data(pointtran)
(dbreaksPt <- seq(0, 25, by=5))

ptUMF <- with(pointtran, {
unmarkedFrame(y = cbind(dc1, dc2, dc3, dc4, dc5), 
siteCovs = data.frame(area, habitat))
})

# Half-normal. Output is animals/point. No covariates.
(fmp1 <- distsamp(~ 1 ~ 1, ptUMF, dist.breaks=dbreaksPt, survey="point", 
unitsIn="m"))

# Negative exponential
(fmpe <- distsamp(~ 1 ~ 1, ptUMF, dist.breaks=dbreaksPt, key="exp", 
survey="point", output="density", unitsIn="m"))}

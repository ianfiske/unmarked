\name{distsamp}
\alias{distsamp}
\title{Fit the hierarchical distance sampling model of Royle et al.}
\usage{distsamp(formula, data, keyfun=c("halfnorm", "exp",
    "hazard", "uniform"), output=c("density", "abund"),
    unitsOut=c("ha", "kmsq"), starts, method=BFGS,
    control=list(), se=TRUE, ...)}
\description{Fit the hierarchical distance sampling model of Royle et al. (2004) to 
line or point transect data recorded in discrete distance intervals.}
\details{Unlike conventional distance sampling, which uses the 'conditional on 
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
\note{Currently, transect-level abundance is assumed to be Poisson distributed 
though other distributions such as the negative binomial may be added. 
Goodness-of-fit can be assessed using the \code{\link{parboot}} function.}
\value{unmarkedFitDS object (child class of \link{unmarkedFit}) describing the 
model fit. Parameter estimates are displayed on the log-scale. 
Back-transformation can be achieved via the \link{predict} or 
\link{backTransform} methods.}
\author{Richard Chandler \email{rchandler@nrc.umass.edu}}
\references{Royle, J. A., D. K. Dawson, and S. Bates (2004) Modeling 
abundance effects in distance sampling. \emph{Ecology} 85, pp. 1591-1597.}
\seealso{\code{\link{unmarkedFit-class}} \code{\link{fitList}} 
\code{\link{parboot}}}
\keyword{models}
\arguments{\item{formula}{'Two-level' right-hand formula describing state-level
covariates followed by detection covariates. ~1 ~1 would be a null model.}
\item{data}{object of class unmarkedFrameDS, containing response matrix, 
covariates, distance interval cut points, survey type ("line" or "point"), 
transect lengths (for survey = "line"), and units ("m" or "km") for cut 
points and transect lengths. See example for set up.}
\item{keyfun}{One of the following detection functions: 
"halfnorm", "hazard", "exp", or "uniform." See details.}
\item{output}{Model either "density" or "abund"}
\item{unitsOut}{Units of density. Either "ha" or "kmsq" for hectares and 
square kilometers, respectively.}
\item{starts}{Vector of starting values for parameters.}
\item{method}{Optimization method used by \code{\link{optim}}.}
\item{control}{Other arguments passed to \code{\link{optim}}.}
\item{se}{logical specifying whether or not to compute standard errors.}}
\examples{## Line transect examples

data(linetran)

ltUMF <- with(linetran, {
unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
siteCovs = data.frame(Length, area, habitat), 
dist.breaks = c(0, 5, 10, 15, 20),
tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
})

ltUMF
plot(ltUMF)
barplot(ltUMF)

# Half-normal detection function. Density output (log scale). No covariates. 
(fm1 <- distsamp(~ 1 ~ 1, ltUMF))

# Some methods to use on fitted model
coef(fm1, type="det", altNames=TRUE)
backTransform(fm1, whichEstimate="det")
vcov(fm1, altNames=TRUE)
confint(fm1, type = "state")
predict(fm1, type = "state")

# Half-normal. Abundance output. No covariates. Note that transect length
# must be accounted for so abundance is animals per km of transect.
(fm2 <- distsamp(~ 1 ~ 1, ltUMF, output="abund", unitsOut="kmsq"))

# Halfnormal. Covariates affecting both density and and detection.  
(fm3.1 <- distsamp(~ poly(area, 2) + habitat ~ habitat, ltUMF))

# This won't run without starting values.
(fm3.2 <- distsamp(~ poly(area, 2) + habitat - 1 ~ habitat - 1, ltUMF, 
dist.breaks=dbreaksLine, starts=c(1,0,0,1,2,2)))

# Negative exponential detection function. This also needs starting values.
(fm4 <- distsamp(~ 1 ~ 1, ltUMF, key="exp", starts=c(0,0)))

# Hazard-rate detection function. Density output in hectares.
# This is real slow, especially for large datasets.
(fmhz <- distsamp(~ 1 ~ 1, ltUMF, keyfun="hazard"))

plot(function(x) gxhaz(x, shape=8.38, scale=1.37), 0, 25, 
xlab="Distance(m)", ylab="Detection probability")

# Uniform detection function. Density output in hectars.
(fmu <- distsamp(~ 1 ~ 1, ltUMF, key="uniform"))

## Point transect example

data(pointtran)

ptUMF <- with(pointtran, {
unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4, dc5), 
siteCovs = data.frame(area, habitat), 
dist.breaks = seq(0, 25, by=5), survey = "point", unitsIn = "m")
})

# Half-normal. Output is animals / ha on log-scale. No covariates.
(fmp1 <- distsamp(~ 1 ~ 1, ptUMF))}

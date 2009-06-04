\name{distsamp}
\alias{distsamp}
\title{Fit the hierarchical distance sampling model...}
\usage{distsamp(stateformula, detformula=~1, data, dist.breaks, tlength,
    keyfun=c("halfnorm", "exp", "hazard", "uniform"), survey,
    output=c("density", "abund"), unitsIn, unitsOut=c("ha", "kmsq"),
    starts, method=BFGS, control=list(), ...)}
\description{Fit the hierarchical distance sampling model}
\details{This functions fits the multinomial-Poisson mixture model for line or point transect data where counts are recorded in discrete distance intervals.}
\keyword{models}
\arguments{\item{stateformula}{Right-hand side formula describing covariates of abundance or density.}
\item{stateformula}{Two-sided formula describing response matrix and covariates of abundance or density.}
\item{detformula}{Right-hand side formula describing covariates of detection.}
\item{data}{data.frame containing response and predictor variables. Each row is a transect.}
\item{dist.breaks}{Numeric vector defining the distance intervals (in meters or km) in which observations were recorded.}
\item{tlength}{Numeric vector of transect lengths in meters or km.}
\item{keyfun}{One of the following detection functions: "halfnorm", "hazard", "exp", or "uniform." See details}
\item{survey}{Either "line" or "point" transect.}
\item{output}{Either "density" or "abund."}
\item{unitsIn}{Either "m" or "km" for BOTH dist.breaks and tlength}
\item{unitsOut}{Units of density. Either "ha" or "kmsq."}
\item{starts}{Vector of starting values for parameters.}
\item{method}{Optimization method used by optim().}
\item{control}{Other arguments passed to optim().}}
\examples{### Line transect examples

## Distance cut points in meters
data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 

## Half-normal detection function. Density output. No covariates. 
## lineDat$Length is transect lengths in km, so it has to be converted.
(fm1 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
tlength=linetran$Length*1000, survey="line", unitsIn="m"))
# backtransform covariates to original scale and use specific names.
# lam(Intecept) is mean density in hectares.
# p(Intercept) is the standard deviation of half-normal function in meters.
exp(coef(fm1, altNames=TRUE))
# variance-covariance matrix
vcov(fm1, altNames=TRUE)


## Half-normal. Abundance output. No covariates. Note that transect length
## must be accounted for so abundance is animals per km of transect.
(fm2 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
tlength=linetran$Length*1000, output="abund", survey="line", unitsIn="m", 
unitsOut="kmsq"))
exp(coef(fm2, altNames=TRUE))

## Halfnormal. Covariates affecting both density and and detection.  
(fm3.1 <- distsamp(cbind(o1,o2,o3,o4) ~ poly(area, 2) + habitat, ~habitat, 
linetran, dist.breaks=dbreaksLine, tlength=linetran$Length*1000, 
survey="line", unitsIn="m", unitsOut="ha"))
exp(coef(fm3.1, altNames=TRUE))

## This won't run without starting values.
(fm3.2 <- distsamp(cbind(o1,o2,o3,o4) ~ poly(area, 2) + habitat - 1, 
~habitat - 1, linetran, dist.breaks=dbreaksLine, tlength=linetran$Length*1000, 
survey="line", unitsIn="m", unitsOut="ha", starts=c(1,0,0,1,2,2)))
exp(coef(fm3.2, altNames=TRUE))


## Negative exponential detection function. Density output in hectares. 
(fme <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
tlength=linetran$Length*1000, key="exp", survey="line", unitsIn="m", 
starts=c(0,-1)))
exp(coef(fme, altNames=TRUE))

## Hazard-rate detection function. Density output in hectares.
## This is real slow, especially for large datasets. Needs to be improved.
(fmhz <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
tlength=linetran$Length*1000, keyfun="hazard", survey="line", unitsIn="m"))
exp(coef(fmhz, altNames=TRUE))
plot(function(x) unmarked:::gxhaz(x, shape=8.38, scale=1.37), 0, 25)
plot(function(x) unmarked:::gxhaz(x, shape=8.38, scale=1.37), 0, 25, 
xlab="Distance(m)", ylab="Detection probability")

## Uniform detection function. Density output in hectars.
(fmu <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
tlength=linetran$Length*1000, key="uniform", survey="line", unitsIn="m"))
exp(coef(fmu, altNames=TRUE))



### Point transect examples

## Radius cut points in meters
data(pointtran)
(dbreaksPt <- seq(0, 25, by=5))

## Half-normal. Output is animals/point. No covariates.
(fmp1 <- distsamp(cbind(o1,o2,o3,o4,o5) ~ 1, ~1, pointtran, 
dist.breaks=dbreaksPt, survey="point", unitsIn="m"))
exp(coef(fmp1, altNames=TRUE))

## Negative exponential
(fmpe <- distsamp(cbind(o1,o2,o3,o4,o5) ~ 1, ~1, pointtran, 
dist.breaks=dbreaksPt, key="exp", survey="point", output="density", 
unitsIn="m"))
exp(coef(fmpe, altNames=TRUE))}

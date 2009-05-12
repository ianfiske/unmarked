\name{distsamp}
\alias{distsamp}
\title{ Distance sampling using unconditional likelihood }
\description{
Fit distance sampling model of Royle et al. (2004), which is a multinomial 
mixture model.
}
\usage{
distsamp(stateformula, detformula = ~1, data, dist.breaks, tlength = NULL, 
keyfun = c("halfnorm", "exp", "hazard", "uniform"), survey = c("line", "point"), 
output = c("density", "abund"), unitsIn = c("m", "km"), 
unitsOut = c("ha", "kmsq"), savedata = T, starts = NULL, hessian = T, 
method = "BFGS", control = list(maxit = 1000, reltol = 1e-20), ...)
}
\arguments{
  \item{stateformula}{ two-sided formula of density/abundance relationship. See details. }
  \item{detormula}{ right-hand sided formula of detection probability covariates. See details }
  \item{data}{ a data.frame containing response and predictor variables. Each row is a transect. }
  \item{dist.breaks}{ numeric vector of distance breaks in meters or km }
  \item{tlength}{ transect length in meters or km }
  \item{keyfun}{ detection function: "halfnorm", "hazard", "exp", or "uniform" }
  \item{survey}{ line transect or point }
  \item{output}{ density or abundance? }
  \item{unitsIn}{ are BOTH dist.breaks and tlength in meters or km }
  \item{unitsOut}{ for density, should it be in hectares or km }
  \item{savedata}{ Should the data.frame be saved in output? }
  \item{starts}{ vector of starting values for parameters }
  \item{hessian}{ return the hessian matrix and calculate SE? }
  \item{method}{ optimization method used by optim }
  \item{\dots}{ Other arguments passed to optim }
}
\details{
  The left-hand side is the columns (e.g. cbind(dist1, dist2, dist3) for the distance data. The right-hand side represents the predictor variables in the same format as used by standard R formulas. 
}
\value{
  \item{par}{Named vector of parameter estimates with attribute altNames}
  \item{value}{Final objective value (-logL) reached by optim}
  \item{counts}{Number of iterations for function and gradient}
  \item{convegence}{optim convergence code}
  \item{message}{optim message}
  \item{hessian}{named hessian matrix returned by optim}
  \item{SE}{standard errors}
  \item{y}{response matrix. lefthand side of stateformula}
  \item{Xlam}{design matrix returned by model.matrix applied to stateformula}
  \item{Xp}{design matrix from detformula}  
  \item{n}{nrow(y)}
  \item{K}{number of estimated parameters}
  \item{J}{ncol(y)}
  \item{keyfun}{keyfun}
  \item{dist.breaks}{dist.breaks in original units}
  \item{area}{area in output units}
  \item{tlength}{transect lengths in input units}
  \item{call}{returned by match.call()}
  \item{stateformula}{stateformula}      
  \item{detformula}{detformula}   
  }
\references{ Royle et al. 2004. Ecology 85:1591-1597}

\seealso{ \code{\link{ds.gof}, \link{lm}, \link{formula}} }

\examples{

### Line transect examples

## Distance cut points in meters
data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20)) 

## Half-normal detection function. Density output. No covariates. 
## lineDat$Length is transect lengths in km, so it has to be converted.
(fm1 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
	tlength=linetran$Length*1000))
# backtransform covariates to original scale and use specific names.
# lam(Intecept) is mean density in hectares.
# p(Intercept) is the standard deviation of half-normal function in meters.
exp(coef(fm1, altNames=T))
# variance-covariance matrix
vcov(fm1)


## Half-normal. Abundance output. No covariates. Note that transect length
## must be accounted for so abundance is animals per km of transect.
(fm2 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
	tlength=linetran$Length*1000, output="abund", unitsOut="kmsq"))
exp(coef(fm2, altNames=T))

## Halfnormal. Covariates affecting both density and and detection.  
(fm3.1 <- distsamp(cbind(o1,o2,o3,o4) ~ poly(area, 2) + habitat, ~habitat, 
   linetran, dist.breaks=dbreaksLine, tlength=linetran$Length*1000, 
   unitsOut="ha"))
exp(coef(fm3.1, altNames=T))

## This won't run without starting values.
(fm3.2 <- distsamp(cbind(o1,o2,o3,o4) ~ poly(area, 2) + habitat - 1, 
   ~habitat - 1, linetran, dist.breaks=dbreaksLine, tlength=linetran$Length*1000, 
   unitsOut="ha", starts=c(1,0,0,1,2,2)))
exp(coef(fm3, altNames=T))


## Negative exponential detection function. Density output in hectares. 
(fme <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
	tlength=linetran$Length*1000, key="exp"))
exp(coef(fme, altNames=T))

## Hazard-rate detection function. Density output in hectares.
## This is real slow, especially for large datasets. Needs to be improved.
(fmhz <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
	tlength=linetran$Length*1000, keyfun="hazard"))
exp(coef(fmhz, altNames=T))
plot(function(x) gxhaz(x, shape=8.38, scale=1.37), 0, 25)

## Uniform detection function. Density output in hectars.
(fmu <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
	tlength=linetran$Length*1000, key="uniform"))
exp(coef(fmu, altNames=T))



### Point transect examples

## Radius cut points in meters
data(pointtran)
(dbreaksPt <- seq(0, 25, by=5))

## Half-normal. Output is animals/point. No covariates.
(fmp1 <- distsamp(cbind(o1,o2,o3,o4,o5) ~ 1, ~1, pointtran, 
	dist.breaks=dbreaksPt, survey="point", output="abund"))
exp(coef(fmp1, altNames=T))	# 15.8 animals / point

## Half-normal. Output is animals per km^2. No covariates.
(fmp2 <- distsamp(cbind(o1,o2,o3,o4,o5) ~ 1, ~1, pointtran, 
	dist.breaks=dbreaksPt, survey="point", output="density", starts=c(8,2)))
exp(coef(fmp2, altNames=T))	#8050 / km^2

## Negative exponential
(fmpe <- distsamp(cbind(o1,o2,o3,o4,o5) ~ 1, ~1, pointtran, 
	dist.breaks=dbreaksPt, key="exp", survey="point", output="density"))
exp(coef(fmpe, altNames=T))








}
\keyword{ models }

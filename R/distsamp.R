# stateformula is 2-sided (eg cbind(catbird1,catbird2,catbird3) ~ log(elev)
# detformula is RHS only (eg ~ poly(area, 2)
# data is a data.frame
# dist.breaks is a vector of distance cut points in meters or km
# tlength is vector of transect lengths in meters or km
# keyfun is the detection function. 
# survey: line transects or points?
# output: standardize by area surveyed or not
# unitsIn: are BOTH dist.breaks and tlength in meters or km?
# unitsOut: hectares or km^2
# savedata? needed for prediction and parametric bootstrapping
# all the rest are passed to optim()


distsamp <- function(stateformula, detformula=~1, data, dist.breaks, 
	tlength=NULL, keyfun=c("halfnorm", "exp", "hazard", "uniform"), 
	survey=c("line", "point"), output=c("density", "abund"), 
	unitsIn=c("m", "km"), unitsOut=c("ha", "kmsq"), savedata=T, 
	starts=NULL, hessian=T, method="BFGS", control=list(maxit=1000, 
	reltol=1e-20), ...)
{
keyfun <- match.arg(keyfun)
survey <- match.arg(survey)
output <- match.arg(output)
unitsIn <- match.arg(unitsIn)
unitsOut <- match.arg(unitsOut)
mf <- model.frame(stateformula, data)
Y <- model.response(mf)
Xlam <- model.matrix(stateformula, data)
Xp <- model.matrix(detformula, data)
n <- nrow(Y)
J <- ncol(Y)
lampars <- colnames(Xlam)
ppars <- colnames(Xp)
nlampars <- length(lampars)
altlampars <- paste(rep("lam", nlampars), colnames(Xlam), sep="")
nppars <- length(ppars)
if(J != length(dist.breaks) - 1)
	stop("ncol(response matrix) must equal length(dist.breaks)-1")
if(dist.breaks[1] != 0)
	stop("dist.breaks[1] must be 0")
switch(unitsIn, 
	km = conv <- 1,
	m = conv <- 1000)
if(output=="density") {
   switch(survey, 
   line = {
	  stripwidths <- (((dist.breaks*2)[-1] - (dist.breaks*2)[-(J+1)])) / conv
	  tl <- tlength / conv
	  a <- rep(tl, each=J) * stripwidths	 # km^2
	  },
   point = {
	  W <- max(dist.breaks) / conv
	  a <- pi * W^2				             # km^2
	  })
   if(unitsOut=="ha")
	  a <- a * 100
   }
   else {
      switch(survey, 
   	  line = a <- tlength / conv, # transect length must be accounted for
	     point = a <- 1)
      }
switch(keyfun, 
halfnorm = { 
   altppars <- paste("sigma", colnames(Xp), sep="")
   if(is.null(starts)) {
	starts <- c(rep(0, length(lampars)), log(max(dist.breaks)), 
		rep(0, length(ppars)-1))
	names(starts) <- c(lampars, ppars)
	}
   else {
	if(is.null(names(starts))) 
		names(starts) <- c(lampars, ppars)
   }
   if(!all(names(starts) %in% c(lampars, ppars)))
	  stop("names(starts) does not agree with necessary model parameters")
   m <- optim(starts, ll.halfnorm, Y=Y, Xlam=Xlam, Xp=Xp, J=J, a=a, 
	  d=dist.breaks, lampars=lampars, survey=survey, method=method, 
	  hessian=hessian, control=control, ...)
   },
exp = { 
   altppars <- paste("rate", colnames(Xp), sep="")
   if(is.null(starts)) {
	  starts <- c(rep(0, length(lampars)), log(median(dist.breaks)), 
		 rep(0, length(ppars)-1))
	  names(starts) <- c(lampars, ppars)
	  }
   else {
	if(is.null(names(starts)))
		names(starts) <- c(lampars, ppars)
	  }
   if(!all(names(starts) %in% c(lampars, ppars)))
	  stop("names(starts) does not agree with necessary model parameters")
   m <- optim(starts, ll.exp, Y=Y, Xlam=Xlam, Xp=Xp, J=J, d=dist.breaks, 
	  a=a, lampars=lampars, survey=survey, method=method, hessian=hessian, 
	  control=control, ...)
	  },
hazard = {	
   ppars <- c(ppars, "scale")
   nppars <- length(ppars)
   altppars <- c(paste("shape", colnames(Xp), sep=""), "scale")
   if(is.null(starts)) {
      starts <- c(rep(0, length(lampars)), log(median(dist.breaks)), 
         rep(0, length(ppars)-2), 1)
      names(starts) <- c(lampars, ppars)
	}
   else {
	if(is.null(names(starts)))
		names(starts) <- c(lampars, ppars)
	}
   if(!all(names(starts) %in% c(lampars, ppars)))
      stop("names(starts) does not agree with necessary model parameters")
   m <- optim(starts, ll.hazard, Y=Y, Xlam=Xlam, Xp=Xp, J=J, d=dist.breaks, 
	  a=a, lampars=lampars, survey=survey, method=method, hessian=hessian, 
	  control=control, ...)
   }, 
uniform = {
   ppars <- character(0)
   nppars <- 0	
   altppars <- character(0)
   if(is.null(starts)) {
      starts <- rep(0, length(lampars))
      names(starts) <- lampars
      }
   else {
	if(is.null(names(starts)))
		names(starts) <- lampars
      }
   if(!all(names(starts) %in% lampars))
      stop("names(starts) does not agree with necessary model parameters")
   m <- optim(starts, ll.uniform, Y=Y, Xlam=Xlam, Xp=Xp, J=J, a=a, 
	  method=method, hessian=hessian, control=control, ...)
   })
attr(m$par, "altNames") <- c(altlampars, altppars)
attr(m$par, "type") <- c(rep("state", nlampars), rep("p", nppars))
if(!is.null(m$hessian)) {
	m$SE <- try(sqrt(diag(solve(m$hessian))))
	if(class(m$SE) != "try-error")
		names(m$SE) <- try(names(m$par))
	}
if(savedata) {
	m$y <- Y
	m$Xlam <- Xlam
	m$Xp <- Xp
	}
m$n <- n
m$K <- length(starts)
m$J <- J
m$keyfun <- keyfun
m$survey <- survey
m$output <- output
m$unitsIn <- unitsIn
m$unitsOut <- unitsOut
m$dist.breaks <- dist.breaks		# original units
m$area <- a		 					   # output units
m$tlength <- tlength					# original units
m$call <- match.call()
m$stateformula <- stateformula
m$detformula <- detformula
class(m) <- c("mixmod", "distsamp")
return(m)
}


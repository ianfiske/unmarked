#' include dsClasses.R
#' include dsUtils.R
roxygen()






setClass("umDistsampFit",
    representation(fitType = "character",
        call = "call",
        stateformula = "formula",
        detformula = "formula",
        data = "data.frame",
        keyfun = "character",
        dist.breaks = "numeric",
        tlength = "numeric",
        area = "numeric",
        survey = "character",
        unitsIn = "character",
        unitsOut = "character",
        estimates = "unMarkedEstimateList",
        AIC = "numeric",
        hessian = "matrix",
        negLogLike = "numeric")
        )





#' Fit the multinomial mixture distance sampling model
#'
#' This functions fits the multinomial-Poisson mixture model for line or point transect data where counts are recorded in discrete distance intervals.
#'
#' 




distsamp <- function(stateformula, detformula=~1, data, dist.breaks, 
	tlength=NULL, keyfun=c("halfnorm", "exp", "hazard", "uniform"), 
	survey=c("line", "point"), output=c("density", "abund"), 
	unitsIn=c("m", "km"), unitsOut=c("ha", "kmsq"), starts=NULL, method="BFGS", 
   control=list(maxit=1000, reltol=1e-20), ...)
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
lamParms <- colnames(Xlam)
detParms <- colnames(Xp)
altlamParms <- paste("lam", colnames(Xlam), sep="")
nAP <- length(lamParms)
nDP <- length(detParms)
nP <- nAP + nDP
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
   if(unitsOut=="ha") a <- a * 100
   }
   else {
      switch(survey, 
	     line = a <- tlength / conv, # transect length must be accounted for
	     point = a <- 1)
      }
switch(keyfun, 
halfnorm = { 
      altdetParms <- paste("sigma", colnames(Xp), sep="")
      if(is.null(starts)) {
      starts <- c(rep(0, nAP), log(max(dist.breaks)), rep(0, nDP-1))
      names(starts) <- c(lamParms, detParms)
      }
   else {
      if(is.null(names(starts))) names(starts) <- c(lamParms, detParms)
      }
   if(!all(names(starts) %in% c(lamParms, detParms)))
      stop("names(starts) does not agree with necessary model parameters")
   fm <- optim(starts, ll.halfnorm, Y=Y, Xlam=Xlam, Xp=Xp, J=J, a=a, 
      d=dist.breaks, nAP=nAP, nP=nP, survey=survey, method=method, hessian=T,
      control=control, ...)
   },
exp = { 
      altdetParms <- paste("rate", colnames(Xp), sep="")
      if(is.null(starts)) {
      starts <- c(rep(0, nAP), log(median(dist.breaks)), rep(0, nDP-1))
      names(starts) <- c(lamParms, detParms)
      }
   else {
	if(is.null(names(starts)))
		names(starts) <- c(lamParms, detParms)
      }
   if(!all(names(starts) %in% c(lamParms, detParms)))
      stop("names(starts) does not agree with necessary model parameters")
   fm <- optim(starts, ll.exp, Y=Y, Xlam=Xlam, Xp=Xp, J=J, d=dist.breaks, 
      a=a, nAP=nAP, nP=nP, survey=survey, method=method, hessian=T, 
      control=control, ...)
	  },
hazard = {	
   detParms <- c(detParms, "scale")
   nDP <- length(detParms)
   altdetParms <- c(paste("shape", colnames(Xp), sep=""), "scale")
   if(is.null(starts)) {
      starts <- c(rep(0, nAP), log(median(dist.breaks)), rep(0, nDP-2), 1)
      names(starts) <- c(lamParms, detParms)
	}
   else {
	if(is.null(names(starts)))
		names(starts) <- c(lamParms, detParms)
	}
   if(!all(names(starts) %in% c(lamParms, detParms)))
      stop("names(starts) does not agree with necessary model parameters")
   fm <- optim(starts, ll.hazard, Y=Y, Xlam=Xlam, Xp=Xp, J=J, d=dist.breaks, 
      a=a, nAP=nAP, nP=nP, survey=survey, method=method, hessian=T, 
      control=control, ...)
   }, 
uniform = {
   detParms <- character(0)
   altdetParms <- character(0)
   nDP <- 0	
   if(is.null(starts)) {
      starts <- rep(0, length(lamParms))
      names(starts) <- lamParms
      }
   else {
	if(is.null(names(starts)))
		names(starts) <- lamParms
      }
   if(!all(names(starts) %in% lamParms))
      stop("names(starts) does not agree with necessary model parameters")
   fm <- optim(starts, ll.uniform, Y=Y, Xlam=Xlam, Xp=Xp, J=J, a=a, 
      method=method, hessian=T, control=control, ...)
   })
ests <- fm$par
estsAP <- ests[1:nAP]
estsDP <- ests[(nAP+1):nP]
attr(estsAP, "altNames") <- altlamParms
attr(estsDP, "altNames") <- altdetParms
tryCatch(covMat <- solve(fm$hessian), 
   error=simpleError("Hessian is not invertible.  Try using fewer covariates."))
attr(fm$hessian, "altNames") <- list(c(altlamParms, altdetParms), 
   c(altlamParms, altdetParms))   
fmAIC <- 2 * fm$value + 2 * nP
stateEstimates <- unMarkedEstimate(name = "Abundance", estimates = estsAP,
   covMat = as.matrix(covMat[1:nAP,1:nAP]), invlink = "exp", 
   invlinkGrad = "exp")
detEstimates <- unMarkedEstimate(name = "Detection", estimates = estsDP, 
   covMat = as.matrix(covMat[(nAP + 1) : nP, (nAP + 1) : nP]), invlink = "exp",
   invlinkGrad = "exp")
estimateList <- unMarkedEstimateList(list(state=stateEstimates, 
   det=detEstimates))
dsfit <- new("umDistsampFit", fitType = "pcount", call = match.call(), 
   stateformula=stateformula, detformula=detformula, data = data, keyfun=keyfun, 
   dist.breaks=dist.breaks, tlength=tlength, area=a, survey=survey, 
   unitsIn=unitsIn, unitsOut=unitsOut, estimates = estimateList, AIC = fmAIC, 
   hessian = fm$hessian, negLogLike = fm$value)
return(dsfit)
}









setMethod("show", "umDistsampFit", function(object)
{
cat("\n", "Call: ", deparse(object@call), "\n", sep="", fill=T)
cat("Coefficients:\n")
show(object@estimates)
})


setGeneric("coef",
		def = function(object, ...) {
			standardGeneric("coef")		})
setGeneric("vcov",
		def = function(object, ...) {
			standardGeneric("vcov")		})


setMethod("coef", "umDistsampFit", function(object, type=NULL, altNames=F)
{
if(is.null(type)) {
   eap <- object@estimates["state"]@estimates
   edp <- object@estimates["det"]@estimates
   e <- c(eap, edp)
   if(altNames)
      names(e) <- c(attr(eap, "altNames"), attr(edp, "altNames"))
   } 
else {
   e <- object@estimates[type]@estimates
   if(altNames)
      names(e) <- attr(e, "altNames")
   }
attr(e, "altNames") <- NULL
return(e)
})



setMethod("vcov", "umDistsampFit", function(object, type=NULL, drop=F, 
   altNames=F)
{
he <- solve(object@hessian)
if(altNames)
   dimnames(he) <- attr(object@hessian, "altNames")
if(!is.null(type)) {
   nAP <- length(coef(object, type="state"))
   nP <- length(coef(object))
   switch(type, 
      state = he <- he[1:nAP, 1:nAP, drop=drop],
      det = he <- he[(nAP+1):nP, (nAP+1):nP, drop=drop])
   }
return(he)
})      


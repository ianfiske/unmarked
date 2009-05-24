#' @include dsUtils.R
#' @include unmarkedFit.R
#' @include unmarkedEstimate.R
#' @include utils.R
roxygen()




#' Model-based prediction 
#'
#' Predict expectations of a model with standard errors (and CI if alpha is specified).
#'
#' If newdata is not specified, original data is used.
#' Requires the deltamethod() function from package msm.
#'
#' @param object a fitted distance sampling model of class 'umDistsampFit'
#' @param type the component of the model to predict. Either "state" or "det"
#' @param link the link function to be used
#' @param newdata an optional data.frame including variables necessary to predict from fitted model
#' @param notconstant an optional column name of a variable in newdata which is not constant. Only useful for plotting expected relationships.
#' @param alpha alpha level for creating confidence intervals
#' 
#' @examples
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#' 
#' #Fit a model
#' (fmhnA.H <- distsamp(cbind(o1,o2,o3,o4) ~ area, ~area + habitat, linetran, 
#' 	dist.breaks=dbreaksLine, tlength=linetran$Length*1000, survey="line", 
#' 	unitsIn="m"))
#'
#' # Create new data.frame for prediction
#' newhabitat <- factor("A")
#' newarea <- seq(min(linetran$area), max(linetran$area), length=20)
#' (ndA <- data.frame(area=newarea, habitat=newhabitat))
#'
#' # Predict density based upon fitted model and 'new' data
#' (Elam.A <- predict(fmhnA.H, newdata=ndA, notconstant="area", type="state", 
#' 	link="log"))
#' with(Elam.A, { # Plot relationships between density and area
#' 	plot(Predictor, Est., ylim=c(0.5, 1.2), xlab="Area", ylab="Density")
#' 	segments(Predictor, Est.-SE, Predictor, Est.+SE)
#' 	})
#' 
#' # Same as above but for detection 
#' ndH <- data.frame(area=mean(linetran$area), habitat=factor(c("a", "b")))
#' ndH
#' 
#' (Ep.H <- predict(fmhnA.H, newdata=ndH, notconstant="habitat", type="det", 
#' 	link="log"))
#' with(Ep.H, {  # Plot relationship (lack of) difference in detectability between habitat types
#' 	bp <- barplot(Est., ylim=c(0, 15), names=Predictor, 
#' 		ylab="Sigma (half-normal shape parameter)") 
#' 	arrows(bp, Est., bp, Est.+SE, code=2, angle=90, length=0.2)
#' 	box()
#' 	})
#' # OR
#' with(Ep.H, {
#' 	plot(function(x) unmarked:::gxhn(x, Est.[1]), 0, 25, 
#' 		xlab="Distance (m)", ylab="Detection probability")
#' 	plot(function(x) unmarked:::gxhn(x, Est.[2]), 0, 25, add=T, lty=2, lwd=2)
#' 	legend(1, 20, c("Habitat a", "Habitat b"), lty=1:2)
#' 	})
#'
#' @exportClass umDistsampFit
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
				estimates = "unmarkedEstimateList",
				AIC = "numeric",
				hessian = "matrix",
				negLogLike = "numeric")
)





#' Fit the hierarchical distance sampling model
#'
#' This functions fits the multinomial-Poisson mixture model for line or point transect data where counts are recorded in discrete distance intervals.
#'
#' 
#' @param stateformula Right-hand side formula describing covariates of abundance or density.
#' @param detformula Right-hand side formula describing covariates of detection.
#' @param data data.frame containing response and predictor variables. Each row is a transect.
#' @param dist.breaks Numeric vector defining the distance intervals (in meters or km) in which observations were recorded.
#' @param tlength Numeric vector of transect lengths in meters or km.
#' @param keyfun One of the following detection functions: "halfnorm", "hazard", "exp", or "uniform." See details
#' @param survey Either "line" or "point" transect.
#' @param output Either "density" or "abund."
#' @param unitsIn Either "m" or "km" for BOTH dist.breaks and tlength
#' @param unitsOut Units of density. Either "ha" or "kmsq."
#' @param starts Vector of starting values for parameters.
#' @param method Optimization method used by optim().
#' @param control Other arguments passed to optim().
#'
#' @examples
#' ### Line transect examples
#' 
#' ## Distance cut points in meters
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#' 
#' ## Half-normal detection function. Density output. No covariates. 
#' ## lineDat$Length is transect lengths in km, so it has to be converted.
#' (fm1 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
#' 	tlength=linetran$Length*1000, survey="line", unitsIn="m"))
#' # backtransform covariates to original scale and use specific names.
#' # lam(Intecept) is mean density in hectares.
#' # p(Intercept) is the standard deviation of half-normal function in meters.
#' exp(coef(fm1, altNames=T))
#' # variance-covariance matrix
#' vcov(fm1)
#' 
#' 
#' ## Half-normal. Abundance output. No covariates. Note that transect length
#' ## must be accounted for so abundance is animals per km of transect.
#' (fm2 <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
#' 	tlength=linetran$Length*1000, output="abund", survey="line", unitsIn="m", 
#'    unitsOut="kmsq"))
#' exp(coef(fm2, altNames=T))
#' 
#' ## Halfnormal. Covariates affecting both density and and detection.  
#' (fm3.1 <- distsamp(cbind(o1,o2,o3,o4) ~ poly(area, 2) + habitat, ~habitat, 
#'    linetran, dist.breaks=dbreaksLine, tlength=linetran$Length*1000, 
#'    survey="line", unitsIn="m", unitsOut="ha"))
#' exp(coef(fm3.1, altNames=T))
#' 
#' ## This won't run without starting values.
#' (fm3.2 <- distsamp(cbind(o1,o2,o3,o4) ~ poly(area, 2) + habitat - 1, 
#'    ~habitat - 1, linetran, dist.breaks=dbreaksLine, tlength=linetran$Length*1000, 
#'    survey="line", unitsIn="m", unitsOut="ha", starts=c(1,0,0,1,2,2)))
#' exp(coef(fm3.2, altNames=T))
#' 
#' 
#' ## Negative exponential detection function. Density output in hectares. 
#' (fme <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
#' 	tlength=linetran$Length*1000, key="exp", survey="line", unitsIn="m", 
#'    starts=c(0,-1)))
#' exp(coef(fme, altNames=T))
#' 
#' ## Hazard-rate detection function. Density output in hectares.
#' ## This is real slow, especially for large datasets. Needs to be improved.
#' (fmhz <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
#' 	tlength=linetran$Length*1000, keyfun="hazard", survey="line", unitsIn="m"))
#' exp(coef(fmhz, altNames=T))
#' plot(function(x) unmarked:::gxhaz(x, shape=8.38, scale=1.37), 0, 25)
#' 
#' ## Uniform detection function. Density output in hectars.
#' (fmu <- distsamp(cbind(o1,o2,o3,o4) ~ 1, ~1, linetran, dist.breaks=dbreaksLine, 
#' 	tlength=linetran$Length*1000, key="uniform", survey="line", unitsIn="m"))
#' exp(coef(fmu, altNames=T))
#' 
#' 
#' 
#' ### Point transect examples
#' 
#' ## Radius cut points in meters
#' data(pointtran)
#' (dbreaksPt <- seq(0, 25, by=5))
#' 
#' ## Half-normal. Output is animals/point. No covariates.
#' (fmp1 <- distsamp(cbind(o1,o2,o3,o4,o5) ~ 1, ~1, pointtran, 
#' 	dist.breaks=dbreaksPt, survey="point", unitsIn="m"))
#' exp(coef(fmp1, altNames=T))
#' 
#' ## Negative exponential
#' (fmpe <- distsamp(cbind(o1,o2,o3,o4,o5) ~ 1, ~1, pointtran, 
#' 	dist.breaks=dbreaksPt, key="exp", survey="point", output="density", 
#'    unitsIn="m"))
#' exp(coef(fmpe, altNames=T))
#' 
#' 
#' 
#' 
#' @export
#' @keywords models
distsamp <- function(stateformula, detformula=~1, data, dist.breaks, 
		tlength=NULL, keyfun=c("halfnorm", "exp", "hazard", "uniform"), 
		survey, output=c("density", "abund"), unitsIn, unitsOut=c("ha", "kmsq"), 
		starts=NULL, method="BFGS", control=list(), ...)
{
	keyfun <- match.arg(keyfun)
	output <- match.arg(output)
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
	if(is.null(tlength)) tlength <- numeric(0)
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
				nP <- nAP + nDP
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
	try(covMat <- solve(fm$hessian))
	covMatAP <- as.matrix(covMat[1:nAP, 1:nAP])
	if(keyfun=="uniform")
		covMatDP <- matrix(numeric(0), 1)
	else
		covMatDP <- as.matrix(covMat[(nAP+1):nP, (nAP+1):nP])
	attr(covMatAP, "altNames") <- list(altlamParms, altlamParms)
	attr(covMatDP, "altNames") <- list(altdetParms, altdetParms)
	attr(fm$hessian, "altNames") <- list(c(altlamParms, altdetParms), 
			c(altlamParms, altdetParms))
	fmAIC <- 2 * fm$value + 2 * nP
	if (keyfun != "uniform") {
		stateEstimates <- unmarkedEstimate(name = "Abundance", estimates = estsAP,
				covMat = covMatAP, invlink = "exp", invlinkGrad = "exp")
		detEstimates <- unmarkedEstimate(name = "Detection", estimates = estsDP, 
				covMat = covMatDP, invlink = "exp", invlinkGrad = "exp")
		estimateList <- unmarkedEstimateList(list(state=stateEstimates, 
						det=detEstimates))
	} else {
		stateEstimates <- unmarkedEstimate(name = "Abundance", estimates = estsAP,
				covMat = covMatAP, invlink = "exp", invlinkGrad = "exp")
		estimateList <- unmarkedEstimateList(list(state=stateEstimates))
	}
	dsfit <- new("umDistsampFit", fitType = "distsamp", call = match.call(), 
			stateformula=stateformula, detformula=detformula, data = data, keyfun=keyfun, 
			dist.breaks=dist.breaks, tlength=tlength, area=a, survey=survey, 
			unitsIn=unitsIn, unitsOut=unitsOut, estimates = estimateList, AIC = fmAIC, 
			hessian = fm$hessian, negLogLike = fm$value)
	return(dsfit)
}



#' @importFrom graphics plot
#' @importFrom stats coef vcov predict update




#' @exportMethod show
setMethod("show", "umDistsampFit", function(object)
		{
			cat("\n", "Call: ", deparse(object@call), "\n", sep="", fill=T)
			cat("Coefficients:\n")
			show(object@estimates)
			cat("Key-function:", object@keyfun)
#cat("\nOutput:", object@output, "\n")
			cat("\nAIC:", object@AIC)
			cat("\nSample size: ", nrow(object@data), "\n")
		})

#
#setGeneric("coef",
#		def = function(object, ...) {
#			standardGeneric("coef")		})
#setGeneric("vcov",
#		def = function(object, ...) {
#			standardGeneric("vcov")		})


#' @exportMethod coef
setMethod("coef", "umDistsampFit", function(object, type=NULL, altNames=F)
		{
			if(is.null(type)) {
				eap <- object@estimates["state"]@estimates
				edp <- NULL
				if(object@keyfun != "uniform") edp <- object@estimates["det"]@estimates
				e <- c(eap, edp)
				if(altNames)
					names(e) <- c(attr(eap, "altNames"), attr(edp, "altNames"))
			}	else {
				e <- object@estimates[type]@estimates
				if(altNames)
					names(e) <- attr(e, "altNames")
			}
			attr(e, "altNames") <- NULL
			return(e)
		})


#' @exportMethod vcov
setMethod("vcov", "umDistsampFit", function(object, type=NULL, drop=F, 
				altNames=F)
		{
			if(is.null(type)) {
				vc <- solve(object@hessian)
				if(altNames)
					dimnames(vc) <- attr(object@hessian, "altNames")
			} 
			else {
				if(!is.null(type)) {
					vc <- object@estimates[type]@covMat
					if(altNames)
						dimnames(vc) <- attr(vc, "altNames")
				}
			}
			attr(vc, "altNames") <- NULL
			return(vc)
		})      


#' @include unmarkedFit.R
#' @include unmarkedEstimate.R
#' @include umDistsampFit.R
#' @include utils.R
{}


#' Fit the hierarchical distance sampling model
#'
#' This functions fits the Royle et al. (2004) multinomial mixture model to 
#' line or point transect data recorded in discrete distance intervals.
#'
#' Unlike conventional distance sampling, which uses the 'conditional on 
#' detection' likelihood formulation, this model is based upon the unconditional 
#' likelihood and thus allows for modeling both abundance and detection 
#' function parameters. 
#'
#' The latent transect-level abundance distribution 
#' \eqn{f(N | \mathbf{\theta})}{f(N | theta)} is currently assumed to be Poisson 
#' with mean \eqn{\lambda}{lambda}.
#'
#' The detection process is modeled as multinomial: 
#' \eqn{y_{ij} \sim Multinomial(N_i, pi_{ij})}{y_ij ~ Multinomial(N_i, pi_i1, pi_i2, ..., pi_iJ)}, 
#' where \eqn{pi_ij} is the multinomial cell probability for transect i in 
#' distance class j. These are computed based upon a detection function 
#' \eqn{g(x | \mathbf{\sigma})}{g(x | sigma)}, such as the half-normal, 
#' negative exponential, or hazard rate.  
#'
#' Parameters \eqn{\lambda}{lambda} and \eqn{\sigma}{sigma} can be vectors 
#' affected by transect-specific covariates using the log link.
#' @note 
#' The response matrix contains the counts of objects at each transect in each 
#' distance interval. These distance intervals must correspond to the distance 
#' break points vector dist.breaks. One must not change dist.breaks without also 
#' changing the response matrix (and vice versa). 
#'
#' Currently, transect-level abundance is assumed to be Poisson distributed 
#' though other distributions such as the negative binomial may be added. 
#' Goodness-of-fit can be assessed using the \code{\link{parboot}} function.  
#'
#' @param stateformula Right-hand or 2-sided side formula describing covariates 
#' of abundance or density. See examples for how to specify the response matrix 
#' in the stateformula when data is a data.frame.
#' @param detformula Right-hand side formula describing covariates of detection.
#' @param data data.frame or unmarkedFrame containing response and predictor 
#' variables. The dimensions of the response matrix must be M x J, where M is 
#' the number of transects and J is the number of distance intervals in which 
#' observation were recorded.
#' @param dist.breaks Numeric vector defining the distance intervals 
#' (in meters or km) in which observations were recorded.
#' @param tlength Numeric vector of transect lengths in meters or km.
#' @param keyfun One of the following detection functions: 
#' "halfnorm", "hazard", "exp", or "uniform." See details.
#' @param survey Either "line" or "point" transect.
#' @param output Model either "density" or "abund" 
#' @param unitsIn Either "m" or "km" for BOTH dist.breaks and tlength
#' @param unitsOut Units of density. Either "ha" or "kmsq" for hectares and 
#' square kilometers, respectively.
#' @param starts Vector of starting values for parameters.
#' @param method Optimization method used by \code{\link{optim}}.
#' @param control Other arguments passed to \code{\link{optim}}.
#'
#' @return umDistsampFit object (child class of \link{unmarkedFit}) describing the 
#' model fit. Parameter estimates are displayed on the log-scale. 
#' Back-transformation can be achieved via the \link{predict} or 
#' \link{backTransform} methods.
#'
#' @author Richard Chandler \email{rchandler@@nrc.umass.edu}
#'
#' @references Royle, J. A., D. K. Dawson, and S. Bates (2004) Modeling 
#' abundance effects in distance sampling. \emph{Ecology} 85, pp. 1591-1597.
#'
#' @seealso \code{\link{umDistsampFit}} \code{\link{unmarkedFitList}}
#'
#' @examples
#' ## Line transect examples
#' 
#' # Distance cut points in meters
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#'
#' lengths <- linetran$Length
#'
#' ltUMF <- with(linetran, {
#' 	unmarkedFrame(y = cbind(dc1, dc2, dc3, dc4), 
#' 		siteCovs = data.frame(Length, area, habitat))
#' 	})
#' 
#' # Half-normal detection function. Density output. No covariates. 
#' # lineDat$Length is transect lengths in km, so it has to be converted.
#'  (fm1 <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
#'  		tlength=lengths*1000, survey="line", unitsIn="m"))
#' 
#' coef(fm1, type="det", altNames=TRUE)
#' backTransform(fm1, whichEstimate="det")
#' vcov(fm1, altNames=TRUE)
#' predict(fm1, type="state")
#' 
#' # Half-normal. Abundance output. No covariates. Note that transect length
#' # must be accounted for so abundance is animals per km of transect.
#' (fm2 <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
#' 		tlength=lengths*1000, output="abund", survey="line", unitsIn="m", 
#' 		unitsOut="kmsq"))
#' 
#' # Halfnormal. Covariates affecting both density and and detection.  
#' (fm3.1 <- distsamp(~ poly(area, 2) + habitat ~ habitat, 
#' 		ltUMF, dist.breaks=dbreaksLine, tlength=lengths*1000, 
#' 		survey="line", unitsIn="m", unitsOut="ha"))
#' 
#' # This won't run without starting values.
#' (fm3.2 <- distsamp(~ poly(area, 2) + habitat - 1 ~ habitat - 1, ltUMF, 
#' 		dist.breaks=dbreaksLine, tlength=lengths*1000, survey="line", 
#' 		unitsIn="m", unitsOut="ha", starts=c(1,0,0,1,2,2)))
#' 
#' # Negative exponential detection function. This one also needs starting values.
#' (fme <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
#' 		tlength=lengths*1000, key="exp", survey="line", unitsIn="m", 
#'		starts=c(0,0)))
#' 
#' # Hazard-rate detection function. Density output in hectares.
#' # This is real slow, especially for large datasets.
#' (fmhz <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
#' 		tlength=lengths*1000, keyfun="hazard", survey="line", unitsIn="m"))
#' 
#' plot(function(x) unmarked:::gxhaz(x, shape=8.38, scale=1.37), 0, 25, 
#' 		xlab="Distance(m)", ylab="Detection probability")
#' 
#' # Uniform detection function. Density output in hectars.
#' (fmu <- distsamp(~ 1 ~ 1, ltUMF, dist.breaks=dbreaksLine, 
#' 		tlength=lengths*1000, key="uniform", survey="line", unitsIn="m"))
#'  
#' ## Point transect examples
#' 
#' # Radius cut points in meters
#' data(pointtran)
#' (dbreaksPt <- seq(0, 25, by=5))
#'
#' ptUMF <- with(pointtran, {
#' 	unmarkedFrame(y = cbind(dc1, dc2, dc3, dc4, dc5), 
#' 		siteCovs = data.frame(area, habitat))
#' 	})
#' 
#' # Half-normal. Output is animals/point. No covariates.
#' (fmp1 <- distsamp(~ 1 ~ 1, ptUMF, dist.breaks=dbreaksPt, survey="point", 
#' 		unitsIn="m"))
#' 
#' # Negative exponential
#' (fmpe <- distsamp(~ 1 ~ 1, ptUMF, dist.breaks=dbreaksPt, key="exp", 
#' 		survey="point", output="density", unitsIn="m"))
#' 
#' @export
#' @keywords models
distsamp <- function(formula, data, dist.breaks, tlength=NULL, 
	keyfun=c("halfnorm", "exp", "hazard", "uniform"), survey, 
	output=c("density", "abund"), unitsIn, unitsOut=c("ha", "kmsq"), 
	starts=NULL, method="BFGS", control=list(), ...)
{
	keyfun <- match.arg(keyfun)
	output <- match.arg(output)
	unitsOut <- match.arg(unitsOut)
	umf <- switch(class(data),
		data.frame = as(data, "unmarkedFrame"),
		unmarkedFrame = data,
		stop("Data is not a data frame or unmarkedFrame."))
	obsToY(umf) <- matrix(1, 1, ncol(y(umf)))
	designMats <- getDesign2(formula, umf)
	X <- designMats$X; V <- designMats$V; y <- designMats$y
	M <- nrow(y)
	J <- ncol(y)
	lamParms <- colnames(X)
	detParms <- colnames(V)
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
	a <- obsCovs(umf)$dcArea
	if(is.null(a))
		a <- calcAreas(dist.breaks, tlength, output, survey, unitsIn, unitsOut)	
	if(is.null(tlength)) 
		tlength <- numeric(0)
	switch(keyfun, 
		halfnorm = { 
			altdetParms <- paste("sigma", colnames(V), sep="")
			if(is.null(starts)) {
				starts <- c(rep(0, nAP), log(max(dist.breaks)), rep(0, nDP-1))
				names(starts) <- c(lamParms, detParms)
			} else {
				if(is.null(names(starts))) names(starts) <- c(lamParms, detParms)
			}
			if(!all(names(starts) %in% c(lamParms, detParms)))
				stop("names(starts) does not agree with necessary model parameters")
			fm <- optim(starts, ll.halfnorm, Y=y, X=X, V=V, J=J, a=a, 
				d=dist.breaks, nAP=nAP, nP=nP, survey=survey, method=method, hessian=T,
				control=control, ...)
		},
		exp = { 
			altdetParms <- paste("rate", colnames(V), sep="")
			if(is.null(starts)) {
				starts <- c(rep(0, nAP), log(median(dist.breaks)), rep(0, nDP-1))
				names(starts) <- c(lamParms, detParms)
			} else {
				if(is.null(names(starts)))
					names(starts) <- c(lamParms, detParms)
			}
			if(!all(names(starts) %in% c(lamParms, detParms)))
				stop("names(starts) does not agree with necessary model parameters")
			fm <- optim(starts, ll.exp, Y=y, X=X, V=V, J=J, d=dist.breaks, 
				a=a, nAP=nAP, nP=nP, survey=survey, method=method, hessian=T, 
				control=control, ...)
		},
		hazard = {	
			detParms <- c(detParms, "scale")
			nDP <- length(detParms)
			nP <- nAP + nDP
			altdetParms <- c(paste("shape", colnames(V), sep=""), "scale")
			if(is.null(starts)) {
				starts <- c(rep(0, nAP), log(median(dist.breaks)), rep(0, nDP-2), 1)
				names(starts) <- c(lamParms, detParms)
			} else {
				if(is.null(names(starts)))
					names(starts) <- c(lamParms, detParms)
			}
			if(!all(names(starts) %in% c(lamParms, detParms)))
				stop("names(starts) does not agree with necessary model parameters")
			fm <- optim(starts, ll.hazard, Y=y, X=X, V=V, J=J, d=dist.breaks, 
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
			} else {
				if(is.null(names(starts)))
					names(starts) <- lamParms
			}
			if(!all(names(starts) %in% lamParms))
				stop("names(starts) does not agree with necessary model parameters")
			fm <- optim(starts, ll.uniform, Y=y, X=X, V=V, J=J, a=a, 
				method=method, hessian=T, control=control, ...)
	})
	opt <- fm
	ests <- fm$par
	estsAP <- ests[1:nAP]
	estsDP <- ests[(nAP+1):nP]
  	covMat <- solve(fm$hessian)
	covMatAP <- covMat[1:nAP, 1:nAP, drop=F]
	if(keyfun=="uniform")
		covMatDP <- matrix(numeric(0), 1)
	else
		covMatDP <- covMat[(nAP+1):nP, (nAP+1):nP, drop=F]
	names(estsDP) <- altdetParms 
	fmAIC <- 2 * fm$value + 2 * nP
	stateEstimates <- unmarkedEstimate(name = "Abundance", 
		short.name = "lam", estimates = estsAP, covMat = covMatAP, 
		invlink = "exp", invlinkGrad = "exp")
	if (keyfun != "uniform") {
		detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p", 
			estimates = estsDP, covMat = covMatDP, invlink = "exp", 
			invlinkGrad = "exp")
		estimateList <- unmarkedEstimateList(list(state=stateEstimates, 
			det=detEstimates))
	} else {
		estimateList <- unmarkedEstimateList(list(state=stateEstimates))
	}

	dsfit <- new("umDistsampFit", fitType = "distsamp", call = match.call(), 
		opt = opt, formula = formula, data = umf, keyfun=keyfun, 
		dist.breaks=dist.breaks, tlength=tlength, area=a, survey=survey, 
		unitsIn=unitsIn, unitsOut=unitsOut, estimates = estimateList, 
		AIC = fmAIC, negLogLike = fm$value)
	return(dsfit)
}


	  		   




## Functions required by distsamp()





#' @title Detection functions used by distsamp()
#' @name detFuns
#' @aliases gxhn gxexp gxhaz grhn grexp grhaz
#' @usage gxhn(x, sigma) gxexp(x, rate) gxhaz(x, shape, scale) 
# @usage gxexp(x, rate)
# @usage gxhaz(x, shape, scale)
# @usage grhn(r, sigma)
# @usage grexp(r, rate)
# @usage grhaz(r, shape, scale)
#' @param x Perpendicular distance
#' @param r Radial distance
#' @param sigma Shape parameter of half-normal detection function
#' @param rate Shape parameter of negative-exponential detection function
#' @param shape Shape parameter of hazard-rate detection function
#' @param scale Scale parameter of hazard-rate detection function
#' @export
gxhn <- function(x, sigma) exp(-x^2/(2 * sigma^2))
#' @export
gxexp <- function(x, rate) exp(-x / rate) 
#' @export
gxhaz <- function(x, shape, scale)  1 - exp(-(x/shape)^-scale)
#' @export
grhn <- function(r, sigma) exp(-r^2/(2 * sigma^2)) * r
#' @export
grexp <- function(r, rate) exp(-r / rate) * r
#' @export
grhaz <- function(r, shape, scale)  (1 - exp(-(r/shape)^-scale)) * r


# Vectorized version of integrate()
vIntegrate <- Vectorize(integrate, c("lower", "upper"))


# Multinomial cell probabilities for line or point transects under half-normal model
cp.hn <- function(d, s, survey=c("line", "point")) 
{
	survey <- match.arg(survey)
	switch(survey, 
			line = {
				strip.widths <- diff(d)
				f.0 <- 2 * dnorm(0, 0, sd=s)
				int <- 2 * (pnorm(d[-1], 0, sd=s) - pnorm(d[-length(d)], 0, sd=s))
				cp <- int / f.0 / strip.widths 
			},
			point = {
				W <- max(d)
				int <- as.numeric(vIntegrate(grhn, d[-length(d)], d[-1], sigma=s)["value",])
				cp <- 2 / W^2 * int
			})
	return(cp)
}



cp.exp <- function(d, rate, survey=c("line", "point")) 
{
	survey <- match.arg(survey)
	switch(survey, 
			line = {
				strip.widths <- diff(d)
				f.0 <- dexp(0, rate=rate)
				int <- pexp(d[-1], rate=rate) - pexp(d[-length(d)], rate=rate)
				cp <- int / f.0 / strip.widths
			},
			point = {
				W <- max(d)
				int <- as.numeric(vIntegrate(grexp, d[-length(d)], d[-1], 
								rate=rate)["value",])
				cp <- 2 / W^2 * int
			})
	return(cp)
}



cp.haz <- function(d, shape, scale, survey=c("line", "point"))
{
	survey <- match.arg(survey)
	switch(survey, 
			line = {
				strip.widths <- diff(d)
				int <- as.numeric(vIntegrate(gxhaz, d[-length(d)], d[-1], shape=shape, 
								scale=scale)["value",])
				cp <- int / strip.widths
			},
			point = {
				W <- max(d)
				int <- as.numeric(vIntegrate(grhaz, d[-length(d)], d[-1], shape=shape, 
								scale=scale)["value",])
				cp <- 2 / W^2 *int
			})
	return(cp)
}






ll.halfnorm <- function(param, Y, X, V, J, d, a, nAP, nP, survey) 
{
	sigma <- as.numeric(exp(V %*% param[(nAP+1):nP]))
	lambda <- as.numeric(exp(X %*% param[1:nAP]))
	pvec <- c(sapply(sigma, function(x) cp.hn(d=d, s=x, survey=survey)))
	growlam <- rep(lambda, each=J)
	datavec <- c(t(Y))
	ll <- dpois(datavec, growlam * pvec * a, log=T)
	-sum(ll)
}




ll.exp <- function(param, Y, X, V, K, J, a, d, nAP, nP, survey=survey)
{
	rate <- as.numeric(exp(V %*% param[(nAP+1):nP]))
	lambda <- as.numeric(exp(X %*% param[1:nAP]))
	pvec <- c(sapply(rate, function(x) cp.exp(d=d, rate=x, survey=survey)))
	growlam <- rep(lambda, each=J)
	datavec <- c(t(Y))
	ll <- dpois(datavec, growlam * pvec * a, log=T)
	-sum(ll)
}


ll.hazard <- function(param, Y, X, V, J, a, d, nAP, nP, survey)
{
	shape <- as.numeric(exp(V %*% param[(nAP+1):(nP-1)]))
	scale <- as.numeric(exp(param[nP]))
	lambda <- as.numeric(exp(X %*% param[1:nAP]))
	pvec <- c(sapply(shape, function(x) cp.haz(d=d, shape=x, scale=scale, 
								survey=survey)))
	growlam <- rep(lambda, each=J)
	datavec <- c(t(Y))
	ll <- dpois(datavec, growlam * a * pvec, log=T)
	-sum(ll)
}




ll.uniform <- function(param, Y, X, V, J, a)
{
	lambda <- exp(X %*% param)
	growlam <- rep(lambda, each=J)
	datavec <- c(t(Y))
	ll <- dpois(datavec, growlam * a, log=T)
	-1*sum(ll)
}


calcAreas <- function(dist.breaks, tlength, output, survey, unitsIn, unitsOut)
{
	switch(unitsIn, 
		km = conv <- 1,
		m = conv <- 1000)
    J <- length(dist.breaks)-1
	switch(output, 
		density =  {
			switch(survey, 
				line = {
					stripwidths <- (((dist.breaks*2)[-1] - (dist.breaks*2)[-(J+1)])) / conv
					tl <- tlength / conv
					a <- rep(tl, each=J) * stripwidths		# km^2
					},
				point = {
					W <- max(dist.breaks) / conv
					a <- pi * W^2							# km^2
					})
			if(unitsOut=="ha") a <- a * 100
			},
		abund = {
			switch(survey, 
				line = a <- tlength / conv, # transect length must be accounted for
				point = a <- 1)
			})
	return(a)
}

#' @include unmarkedFit.R
#' @include unmarkedEstimate.R
#' @include utils.R
{}


#' Fit the hierarchical distance sampling model of Royle et al. (2004) to 
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
#'
#' @note 
#'
#' Currently, transect-level abundance is assumed to be Poisson distributed 
#' though other distributions such as the negative binomial may be added. 
#' Goodness-of-fit can be assessed using the \code{\link{parboot}} function.  
#'
#' @param formula 'Two-level' right-hand formula describing state-level
#' covariates followed by detection covariates. ~1 ~1 would be a null model.
#' @param data object of class unmarkedFrameDS, containing response matrix, 
#' covariates, distance interval cut points, survey type ("line" or "point"), 
#' transect lengths (for survey = "line"), and units ("m" or "km") for cut 
#' points and transect lengths. See example for set up.
#' @param keyfun One of the following detection functions: 
#' "halfnorm", "hazard", "exp", or "uniform." See details.
#' @param output Model either "density" or "abund" 
#' @param unitsOut Units of density. Either "ha" or "kmsq" for hectares and 
#' square kilometers, respectively.
#' @param starts Vector of starting values for parameters.
#' @param method Optimization method used by \code{\link{optim}}.
#' @param control Other arguments passed to \code{\link{optim}}.
#' @param se logical specifying whether or not to compute standard errors.
#' @return unmarkedFitDS object (child class of \link{unmarkedFit}) describing the 
#' model fit. Parameter estimates are displayed on the log-scale. 
#' Back-transformation can be achieved via the \link{predict} or 
#' \link{backTransform} methods.
#'
#' @author Richard Chandler \email{rchandler@@nrc.umass.edu}
#'
#' @references Royle, J. A., D. K. Dawson, and S. Bates (2004) Modeling 
#' abundance effects in distance sampling. \emph{Ecology} 85, pp. 1591-1597.
#'
#' @seealso \code{\link{unmarkedFit}} \code{\link{fitList}} 
#' \code{\link{parboot}}
#'
#' @examples
#' ## Line transect examples
#' 
#' data(linetran)
#'
#' ltUMF <- with(linetran, {
#' 		unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4), 
#' 		siteCovs = data.frame(Length, area, habitat), 
#'		dist.breaks = c(0, 5, 10, 15, 20),
#'		tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
#' 		})
#'
#' ltUMF
#' plot(ltUMF)
#' barplot(ltUMF)
#' 
#' # Half-normal detection function. Density output (log scale). No covariates. 
#'  (fm1 <- distsamp(~ 1 ~ 1, ltUMF))
#' 
#' # Some methods to use on fitted model
#' coef(fm1, type="det", altNames=TRUE)
#' backTransform(fm1, whichEstimate="det")
#' vcov(fm1, altNames=TRUE)
#' confint(fm1, type = "state")
#' predict(fm1, type = "state")
#' 
#' # Half-normal. Abundance output. No covariates. Note that transect length
#' # must be accounted for so abundance is animals per km of transect.
#' (fm2 <- distsamp(~ 1 ~ 1, ltUMF, output="abund", unitsOut="kmsq"))
#' 
#' # Halfnormal. Covariates affecting both density and and detection.  
#' (fm3.1 <- distsamp(~ poly(area, 2) + habitat ~ habitat, ltUMF))
#' 
#' # This won't run without starting values.
#' (fm3.2 <- distsamp(~ poly(area, 2) + habitat - 1 ~ habitat - 1, ltUMF, 
#' 		dist.breaks=dbreaksLine, starts=c(1,0,0,1,2,2)))
#' 
#' # Negative exponential detection function. This also needs starting values.
#' (fm4 <- distsamp(~ 1 ~ 1, ltUMF, key="exp", starts=c(0,0)))
#' 
#' # Hazard-rate detection function. Density output in hectares.
#' # This is real slow, especially for large datasets.
#' (fmhz <- distsamp(~ 1 ~ 1, ltUMF, keyfun="hazard"))
#' 
#' plot(function(x) gxhaz(x, shape=8.38, scale=1.37), 0, 25, 
#' 		xlab="Distance(m)", ylab="Detection probability")
#' 
#' # Uniform detection function. Density output in hectars.
#' (fmu <- distsamp(~ 1 ~ 1, ltUMF, key="uniform"))
#'  
#' ## Point transect example
#' 
#' data(pointtran)
#'
#' ptUMF <- with(pointtran, {
#' 		unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4, dc5), 
#' 		siteCovs = data.frame(area, habitat), 
#'		dist.breaks = seq(0, 25, by=5), survey = "point", unitsIn = "m")
#' 		})
#' 
#' # Half-normal. Output is animals / ha on log-scale. No covariates.
#' (fmp1 <- distsamp(~ 1 ~ 1, ptUMF))
#' 
#' @export
#' @keywords models
distsamp <- function(formula, data, 
	keyfun=c("halfnorm", "exp", "hazard", "uniform"), 
	output=c("density", "abund"), unitsOut=c("ha", "kmsq"), starts=NULL, 
	method="BFGS", control=list(), se = TRUE, ...)
{
	keyfun <- match.arg(keyfun)
	output <- match.arg(output)
	unitsOut <- match.arg(unitsOut)
	dist.breaks <- data@dist.breaks
	tlength <- data@tlength
	survey <- data@survey
	unitsIn <- data@unitsIn
	if(all(is.na(data@plotArea))) {
		a <- calcAreas(dist.breaks = dist.breaks, tlength = tlength, 
			survey = survey, output = output, M = numSites(data), 
			J = ncol(getY(data)), unitsIn = unitsIn, unitsOut = unitsOut)
		a <- c(t(a))
		data@plotArea <- a
		}
	designMats <- getDesign2(formula, data)
	X <- designMats$X; V <- designMats$V; y <- designMats$y
	a <- designMats$plotArea
	M <- nrow(y)
	J <- ncol(y)
	lamParms <- colnames(X)
	detParms <- colnames(V)
	nAP <- length(lamParms)
	nDP <- length(detParms)
	nP <- nAP + nDP
	switch(keyfun,
		halfnorm = { 
			altdetParms <- paste("sigma", colnames(V), sep="")
			if(is.null(starts)) {
				starts <- c(rep(0, nAP), log(max(dist.breaks)), rep(0, nDP-1))
				names(starts) <- c(lamParms, detParms)
			} else {
				if(is.null(names(starts))) names(starts) <- c(lamParms, 
					detParms)
				}
				nll <- function(params) {
					ll.halfnorm(params, Y=y, X=X, V=V, J=J, a=a, 
						d=dist.breaks, nAP=nAP, nP=nP, survey=survey)
				}
		},
		exp = { 
			altdetParms <- paste("rate", colnames(V), sep="")
			if(is.null(starts)) {
				starts <- c(rep(0, nAP), log(median(dist.breaks)), 
					rep(0, nDP-1))
				names(starts) <- c(lamParms, detParms)
			} else {
				if(is.null(names(starts)))
					names(starts) <- c(lamParms, detParms)
			}
			nll <- function(params) {
				ll.exp(params,  Y=y, X=X, V=V, J=J, d=dist.breaks, 
					a=a, nAP=nAP, nP=nP, survey=survey)
			}
		},
		hazard = {	
			detParms <- c(detParms, "scale")
			nDP <- length(detParms)
			nP <- nAP + nDP
			altdetParms <- c(paste("shape", colnames(V), sep=""), "scale")
			if(is.null(starts)) {
				starts <- c(rep(0, nAP), log(median(dist.breaks)), 
					rep(0, nDP-2), 1)
				names(starts) <- c(lamParms, detParms)
			} else {
				if(is.null(names(starts)))
					names(starts) <- c(lamParms, detParms)
			}
			nll <- function(params) {
				ll.hazard(params, Y=y, X=X, V=V, J=J, d=dist.breaks, 
						a=a, nAP=nAP, nP=nP, survey=survey)
			}
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
			nll <- function(params) {
				ll.uniform(params, Y=y, X=X, V=V, J=J, a=a)
			}
		})
	fm <- optim(starts, nll, method=method, hessian=se, control=control)
	opt <- fm
	ests <- fm$par
	estsAP <- ests[1:nAP]
	estsDP <- ests[(nAP+1):nP]
	if(se) {
		tryCatch(covMat <- solve(fm$hessian),
				error=function(x) simpleError("Hessian is not invertible.  Try using fewer covariates."))
	} else {
		covMat <- matrix(NA, nP, nP)
	}
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
	dsfit <- new("unmarkedFitDS", fitType = "distsamp", call = match.call(), 
		opt = opt, formula = formula, data = data, keyfun=keyfun, 
		sitesRemoved = designMats$removed.sites, unitsOut=unitsOut, 
		estimates = estimateList, AIC = fmAIC, negLogLike = fm$value, 
		nllFun = nll)
	return(dsfit)
}





## Functions required by distsamp() and a utility





#' @title Detection functions used by distsamp()
#' @name detFuns
#' @aliases detFuns
#' @aliases gxhn 
#' @aliases gxexp 
#' @aliases gxhaz 
#' @aliases grhn 
#' @aliases grexp 
#' @aliases grhaz
#' @param x Perpendicular distance
#' @param r Radial distance
#' @param sigma Shape parameter of half-normal detection function
#' @param rate Shape parameter of negative-exponential detection function
#' @param shape Shape parameter of hazard-rate detection function
#' @param scale Scale parameter of hazard-rate detection function
#' @export gxhn
gxhn <- function(x, sigma) exp(-x^2/(2 * sigma^2))
#' @export gxexp
gxexp <- function(x, rate) exp(-x / rate) 
#' @export gxhaz
gxhaz <- function(x, shape, scale)  1 - exp(-(x/shape)^-scale)
#' @export grhn
grhn <- function(r, sigma) exp(-r^2/(2 * sigma^2)) * r
#' @export grexp
grexp <- function(r, rate) exp(-r / rate) * r
#' @export grhaz
grhaz <- function(r, shape, scale)  (1 - exp(-(r/shape)^-scale)) * r


# Vectorized version of integrate()
vIntegrate <- Vectorize(integrate, c("lower", "upper"))


# Multinomial cell probabilities for line or point transects under half-normal model
cp.hn <- function(d, s, survey) 
{
	switch(survey, 
		line = {
			strip.widths <- diff(d)
			f.0 <- 2 * dnorm(0, 0, sd=s)
			int <- 2 * (pnorm(d[-1], 0, sd=s) - pnorm(d[-length(d)], 0, sd=s))
			cp <- int / f.0 / strip.widths 
		},
		point = {
			W <- max(d)
			int <- as.numeric(vIntegrate(grhn, d[-length(d)], d[-1], 
				sigma=s)["value",])
			cp <- 2 / W^2 * int
		})
	return(cp)
}



cp.exp <- function(d, rate, survey) 
{
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



cp.haz <- function(d, shape, scale, survey)
{
	switch(survey, 
		line = {
			strip.widths <- diff(d)
			int <- as.numeric(vIntegrate(gxhaz, d[-length(d)], d[-1], 
				shape=shape, scale=scale)["value",])
			cp <- int / strip.widths
		},
		point = {
			W <- max(d)
			int <- as.numeric(vIntegrate(grhaz, d[-length(d)], d[-1], 
				shape=shape, scale=scale)["value",])
			cp <- 2 / W^2 * int
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




ll.exp <- function(param, Y, X, V, K, J, a, d, nAP, nP, survey)
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


#' Prepare area argument for distsamp(). This is primarily for internal use
#' but see details. Caution should be used because the returned matrix has 
#' different interpretations for different survey and output types.
#'
#' @param dist.breaks numeric vector of distance class break poings
#' @param tlength numeric vector of transect lengths for line transects
#' @param M number of transects
#' @param J number of distance classes
#' @param survey either "line" or "point"
#' @param output either "abund" or "density"
#' @param unitsIn either "m" or "km" for units of both dist.breaks and tlength.
#' @param unitsOut either "ha" or "kmsq"
#'
#' @return An M x J numeric matrix. 
#' 
#' If output == "density" and survey == "line"
#' then the values are the areas of each distance class for each transect. If
#' output == "density" and survey == "point" then the values are the the radii
#' of each point transect. Currently, radii cannot vary.
#'
#' If survey == "point" and output == "abund" a matrix of 1s is returned. If
#' survey == "line" and output == "abund" a matrix of transect lengths is
#' returned because transect lengths must be taken into account even if density
#' is not of interest.
#'
#' @note This function might be useful if some distance classes for some 
#' transects were not surveyed. In such a case, missing values could be added
#' to the output of calcAreas() and the modified matrix could be supplied
#' to the plotArea argument of unmarkedFrameDS.
#' 
#' @examples
#'
#' data(linetran)
#' (dbreaksLine <- c(0, 5, 10, 15, 20)) 
#' lengths <- linetran$Length * 1000
#' 
#' calcAreas(dbreaksLine, lengths, "line", "density", M=nrow(linetran), 
#' 		J = length(dbreaksLine) - 1, "m", "ha")
#'
#' data(pointtran)
#' (dbreaksPt <- seq(0, 25, by=5))
#' 
#' calcAreas(dbreaksPt, survey="point", output="density", M=nrow(pointtran),
#' 		J = length(dbreaksPt) - 1, unitsIn = "m", unitsOut = "ha")
#' 
#' @export 
calcAreas <- function(dist.breaks, tlength, survey, output, M, J, unitsIn, 
	unitsOut)
{
if(J != length(dist.breaks) - 1)
	stop("J must equal length(dist.breaks) - 1")
if(!missing(tlength))
	if(length(tlength) > 0 & length(tlength) != M)
		stop("length(tlength) must equal M")
switch(unitsIn, 
	km = conv <- 1,
	m = conv <- 1000
	)
switch(output, 	
	density = {
		switch(survey, 
			line = {
				stripwidths <- (((dist.breaks*2)[-1] - 
					(dist.breaks*2)[-(J+1)])) / conv
				tl <- tlength / conv
				a <- rep(tl, each=J) * stripwidths						# km^2
				a <- matrix(a, nrow=M, ncol=J, byrow=TRUE)
				if(unitsOut == "ha") a <- a * 100
				},
			point = {
				W <- max(dist.breaks) / conv
				a <- matrix(rep(pi * W^2, each=J), M, J, byrow=TRUE) 	# km^2
				if(unitsOut == "ha") a <- a * 100			
				})
			},
	abund = { 
		switch(survey, 
			line = a <- matrix(rep(tlength, each=J), M, J, byrow=TRUE),
			point = a <- matrix(1, M, J))
		})
return(a)
}





#' Convert individual-level distance data to the transect-level format required 
#' by distsamp()
#' 
#' This function creates a site (M) by distance interval (J) response matrix
#' from a data.frame containing the detection distances for each individual and
#' the transect names.
#'
#' @param distData data.frame where each row is a detected individual. Must have
#' at least 2 columns. One for distances and the other for transect names.
#' @param distCol character, the column name containing distances
#' @param transectNameCol character, column name containing transect names
#' @param dist.breaks numeric vector of distance interval cutpoints. Length must
#' equal J+1.
#'
#' @return
#' An M x J data.frame containing the tabulated detections in each distance
#' interval for each transect. Transect names will become rownames and 
#' colnames will be y.1, y.2, ..., y.J. 
#'
#' @note
#' It is very important that the factor containing transect names contains
#' levels for all the transects surveyed. This includes those where birds were
#' not detected. See the example for how to add levels to a factor.
#'
#' @seealso unmarkedFrameDS distsamp
#'
#' @examples 
#' # Create a data.frame containing distances of animals detected
#' # along 4 transects.
#' dat <- data.frame(transect=gl(4,5, labels=letters[1:4]), 
#'		distance=rpois(20, 10))
#' dat
#'
#' # Look at your transect names.
#' levels(dat$transect)
#' 
#' # Suppose that you also surveyed a transect named "e" where no animals were
#' # detected. You must add it to the levels of dat$transect
#' levels(dat$transect) <- c(levels(dat$transect), "e")
#' levels(dat$transect)
#' 
#' # Distance cut points defining distance intervals
#' cp <- c(6, 8, 10, 12, 14, 18)
#' 
#' # Create formated response data.frame
#' yDat <- formatDistData(dat, "distance", "transect", cp) 
#' yDat
#'
#' # Now you could merge yDat with transect-level covariates and 
#' # then use unmarkedFrameDS to prepare data for distsamp
#' @export
formatDistData <- function(distData, distCol, transectNameCol, dist.breaks)
{
transects <- distData[,transectNameCol]
M <- nlevels(transects)
J <- length(dist.breaks) - 1
y <- matrix(NA, M, J, 
	dimnames = list(levels(transects), paste("y", 1:J, sep=".")))
for(i in 1:M) {
	sub <- subset(distData, transects==rownames(y)[i])
	y[i,] <- table(cut(sub[,distCol], dist.breaks, include.lowest=TRUE))
	}
return(data.frame(y))
}




 

# Royel et al. 2004 distance sampling

distsamp <- function(formula, data, 
	keyfun=c("halfnorm", "exp", "hazard", "uniform"), 
	output=c("density", "abund"), unitsOut=c("ha", "kmsq"), starts=NULL, 
	method="BFGS", control=list(), se = TRUE)
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
    if(output == "abund" & length(table(tlength)) > 1)
        warning("Response is individuals per unit transect length")		
	designMats <- getDesign2(formula, data)
	X <- designMats$X; V <- designMats$V; y <- designMats$y
	yvec <- as.numeric(t(y))
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
					ll.halfnorm(params, yvec=yvec, X=X, V=V, J=J, a=a, 
						d=dist.breaks, nAP=nAP, nP=nP, survey=survey)
				}
		},
		exp = { 
			altdetParms <- paste("rate", colnames(V), sep="")
			if(is.null(starts)) {
				starts <- c(rep(0, nAP), 0, rep(0, nDP-1))
				names(starts) <- c(lamParms, detParms)
			} else {
				if(is.null(names(starts)))
					names(starts) <- c(lamParms, detParms)
			}
			nll <- function(params) {
				ll.exp(params,  yvec=yvec, X=X, V=V, J=J, d=dist.breaks, 
					a=a, nAP=nAP, nP=nP, survey=survey)
			}
		},
		hazard = {	
			nDP <- length(detParms)
			nP <- nAP + nDP + 1
			altdetParms <- paste("shape", colnames(V), sep="")
			if(is.null(starts)) {
				starts <- c(rep(0, nAP), log(median(dist.breaks)), 
					rep(0, nDP-1), 1)
				names(starts) <- c(lamParms, detParms, "scale")
			} else {
				if(is.null(names(starts)))
					names(starts) <- c(lamParms, detParms, "scale")
			}
			nll <- function(params) {
				ll.hazard(params, yvec=yvec, X=X, V=V, J=J, d=dist.breaks, 
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
				ll.uniform(params, yvec=yvec, X=X, V=V, J=J, a=a)
			}
		})
	fm <- optim(starts, nll, method=method, hessian=se, control=control)
	opt <- fm
	ests <- fm$par
	if(se) {
		tryCatch(covMat <- solve(fm$hessian),
				error=function(x) simpleError("Hessian is not invertible.  Try using fewer covariates."))
	} else {
		covMat <- matrix(NA, nP, nP)
	}
	estsAP <- ests[1:nAP]
	if(keyfun == "hazard") {
		estsDP <- ests[(nAP+1):(nP-1)]
		estsScale <- ests[nP]
		}
	else 
		estsDP <- ests[(nAP+1):nP]
	covMatAP <- covMat[1:nAP, 1:nAP, drop=F]
	if(keyfun=="hazard") {
		covMatDP <- covMat[(nAP+1):(nP-1), (nAP+1):(nP-1), drop=F]
		covMatScale <- covMat[nP, nP, drop=F]
		}
	else if(keyfun!="uniform")
		covMatDP <- covMat[(nAP+1):nP, (nAP+1):nP, drop=F]
	names(estsDP) <- altdetParms 
	fmAIC <- 2 * fm$value + 2 * nP
	stateName <- switch(output, 
		abund = "Abundance",
		density = "Density")
	stateEstimates <- unmarkedEstimate(name = stateName, 
		short.name = "lam", estimates = estsAP, covMat = covMatAP, 
		invlink = "exp", invlinkGrad = "exp")
	if (keyfun != "uniform") {
		detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p", 
			estimates = estsDP, covMat = covMatDP, invlink = "exp", 
			invlinkGrad = "exp")
		if(keyfun != "hazard")
			estimateList <- unmarkedEstimateList(list(state=stateEstimates, 
				det=detEstimates))
		else {
			scaleEstimates <- unmarkedEstimate(name = "Detection(scale)", 
				short.name = "p", estimates = estsScale, 
				covMat = covMatScale, invlink = "exp", invlinkGrad = "exp")
			estimateList <- unmarkedEstimateList(list(state=stateEstimates, 
				det=detEstimates, scale=scaleEstimates))
			}			
	} else {
		estimateList <- unmarkedEstimateList(list(state=stateEstimates))
	}
	dsfit <- new("unmarkedFitDS", fitType = "distsamp", call = match.call(), 
		opt = opt, formula = formula, data = data, keyfun=keyfun, 
		sitesRemoved = designMats$removed.sites, unitsOut=unitsOut, 
		estimates = estimateList, AIC = fmAIC, negLogLike = fm$value, 
		nllFun = nll, output=output)
	return(dsfit)
}


# Detection functions

gxhn <- function(x, sigma) exp(-x^2/(2 * sigma^2))
gxexp <- function(x, rate) exp(-x / rate) 
gxhaz <- function(x, shape, scale)  1 - exp(-(x/shape)^-scale)
grhn <- function(r, sigma) exp(-r^2/(2 * sigma^2)) * r
grexp <- function(r, rate) exp(-r / rate) * r
grhaz <- function(r, shape, scale)  (1 - exp(-(r/shape)^-scale)) * r

dxhn <- function(x, sigma) 
	gxhn(x=x, sigma=sigma) / integrate(gxhn, 0, Inf, sigma=sigma)$value
drhn <- function(r, sigma) 
	grhn(r=r, sigma=sigma) / integrate(grhn, 0, Inf, sigma=sigma)$value
dxexp <- function(x, rate) 
	gxexp(x=x, rate=rate) / integrate(gxexp, 0, Inf, rate=rate)$value
drexp <- function(r, rate) 
	grexp(r=r, rate=rate) / integrate(grexp, 0, Inf, rate=rate)$value
dxhaz <- function(x, shape, scale)
	gxhaz(x=x, shape=shape, scale=scale) / integrate(gxhaz, 0, Inf, 
		shape=shape, scale=scale)$value
drhaz <- function(r, shape, scale)
	grhaz(r=r, shape=shape, scale=scale) / integrate(grhaz, 0, Inf, 
		shape=shape, scale=scale)$value
		



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
#			f.0 <- dexp(0, rate=rate)
#			int <- pexp(d[-1], rate=rate) - pexp(d[-length(d)], rate=rate)
			int <- as.numeric(vIntegrate(gxexp, d[-length(d)], d[-1],
				rate=rate)["value",])
			cp <- int / strip.widths
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




# Likelihood functions

ll.halfnorm <- function(param, yvec, X, V, J, d, a, nAP, nP, survey) 
{
	sigma <- as.numeric(exp(V %*% param[(nAP+1):nP]))
	lambda <- as.numeric(exp(X %*% param[1:nAP]))
	pvec <- c(sapply(sigma, function(x) cp.hn(d=d, s=x, survey=survey)))
	growlam <- rep(lambda, each=J)
	ll <- dpois(yvec, growlam * pvec * a, log=T)
	-sum(ll)
}




ll.exp <- function(param, yvec, X, V, K, J, a, d, nAP, nP, survey)
{
	rate <- as.numeric(exp(V %*% param[(nAP+1):nP]))
	lambda <- as.numeric(exp(X %*% param[1:nAP]))
	pvec <- c(sapply(rate, function(x) cp.exp(d=d, rate=x, survey=survey)))
	growlam <- rep(lambda, each=J)
	ll <- dpois(yvec, growlam * pvec * a, log=T)
	-sum(ll)
}




ll.hazard <- function(param, yvec, X, V, J, a, d, nAP, nP, survey)
{
	shape <- as.numeric(exp(V %*% param[(nAP+1):(nP-1)]))
	scale <- as.numeric(exp(param[nP]))
	lambda <- as.numeric(exp(X %*% param[1:nAP]))
	pvec <- c(sapply(shape, function(x) cp.haz(d=d, shape=x, scale=scale, 
		survey=survey)))
	growlam <- rep(lambda, each=J)
	ll <- dpois(yvec, growlam * a * pvec, log=T)
	-sum(ll)
}




ll.uniform <- function(param, yvec, X, V, J, a)
{
	lambda <- exp(X %*% param)
	growlam <- rep(lambda, each=J)
	ll <- dpois(yvec, growlam * a, log=T)
	-1*sum(ll)
}


# Prepare area argument for distsamp(). This is primarily for internal use

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





# Convert individual-level distance data to the transect-level format required 
# by distsamp()

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



## Sight distance to perpendicular distance

sight2perpdist <- function(sightdist, sightangle) 
{
if(any(0 > sightangle | sightangle > 180))
	stop("sightangle must be degrees in [0, 180]")
sightdist * sin(sightangle * pi / 180)
}
 
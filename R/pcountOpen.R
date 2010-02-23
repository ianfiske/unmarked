


# Class
setClass("unmarkedFramePCountOpen",
		representation(
			delta = "matrix"),
		contains = "unmarkedFrame")
		





# Constructor
unmarkedFramePCountOpen <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo,
 	plotArea = NULL, delta) 
{
	J <- ncol(y)
	if(missing(plotArea) || is.null(plotArea)) plotArea <- rep(1, nrow(y))		
	if(missing(delta))
		delta <- matrix(1, nrow(y), ncol(y) - 1)
	if(nrow(delta) != nrow(y) | ncol(delta) != ncol(y) - 1)
		stop("Dimensions of delta matrix should be nrow(y), ncol(y)-1")
	if(class(obsCovs) == "list") {
		obsVars <- names(obsCovs)
    for(i in seq(length(obsVars))) {
    	if(!(class(obsCovs[[i]]) %in% c("matrix", "data.frame")))
        	stop("At least one element of obsCovs is not a matrix or data frame.")
    	if(ncol(obsCovs[[i]]) != ncol(y) | nrow(obsCovs[[i]]) != nrow(y))
        	stop("At least one matrix in obsCovs has incorrect number of dimensions.")
    	}
	if(is.null(obsNum)) obsNum <- ncol(obsCovs[[1]])
	obsCovs <- data.frame(lapply(obsCovs, function(x) as.vector(t(x))))
  	}
	umf <- new("unmarkedFramePCountOpen", y = y, siteCovs = siteCovs, 
		obsCovs = obsCovs, obsToY = diag(J), 
		plotArea = plotArea, delta = delta)
	return(umf)
}






# Design matrices
getDesign4 <- function(formula, umf, na.rm = TRUE) 
{
	
lamformula <- formula$lambdaformula
gamformula <- formula$gammaformula
omformula <- formula$omegaformula
pformula <- formula$pformula

delta <- umf@delta
M <- numSites(umf)
R <- obsNum(umf)
	
if(is.null(siteCovs(umf))) {
	siteCovs <- data.frame(placeHolder = rep(1, M))
	} else {
	siteCovs <- siteCovs(umf)
	}
Xlam.mf <- model.frame(lamformula, siteCovs, na.action = NULL)
Xlam <- model.matrix(lamformula, Xlam.mf)

if(is.null(obsCovs(umf))) {
	obsCovs <- data.frame(placeHolder = rep(1, M*R))
	} else {
	obsCovs <- obsCovs(umf)
	}
	
colNames <- c(colnames(obsCovs), colnames(siteCovs))
	
obsCovs <- cbind(obsCovs, siteCovs[rep(1:M, each = R),])
colnames(obsCovs) <- colNames
	
if(!("obs" %in% names(obsCovs))) {
	obsCovs <- cbind(obsCovs, obs = as.factor(rep(1:R, M)))
	}
	
Xp.mf <- model.frame(pformula, obsCovs, na.action = NULL)
Xp <- model.matrix(pformula, Xp.mf)
Xgam.mf <- model.frame(gamformula, obsCovs, na.action = NULL)
Xgam <- model.matrix(gamformula, Xgam.mf)
Xom.mf <- model.frame(omformula, obsCovs, na.action = NULL)
Xom <- model.matrix(omformula, Xom.mf)
	
if(na.rm)
	out <- handleNA4(umf, Xlam, Xgam, Xom, Xp)
else
	out <- list(y=getY(umf), Xlam=Xlam, Xgam=Xgam, Xom=Xom, Xp=Xp, 
		plotArea=umf@plotArea, delta=umf@delta, removed.sites=integer(0))
	
return(list(y = out$y, Xlam = out$Xlam, Xgam = out$Xgam, Xom = Xom, Xp = Xp,
	plotArea = out$plotArea, delta = out$delta, 
	removed.sites = out$removed.sites))
}








# NAs
handleNA4 <- function(umf, Xlam, Xgam, Xom, Xp) {
	# Perhaps all rows with NAs need to be removed b/c of recursive algorithm?
	obsToY <- obsToY(umf)
	if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")
	
	J <- numY(umf)
	R <- obsNum(umf)
	M <- numSites(umf)

	plotArea <- umf@plotArea
	plotArea.na <- is.na(plotArea)
	delta <- umf@delta
	
	Xlam.long <- Xlam[rep(1:M, each = J),]
	Xlam.long.na <- is.na(Xlam.long)
	
	long.na <- function(x) {
		x.mat <- matrix(x, M, R, byrow = TRUE)
		x.mat <- is.na(x.mat)
		x.mat <- x.mat %*% obsToY
		x.long <- as.vector(t(x.mat))
		x.long == 1
		}
	
	Xp.long.na <- apply(Xp, 2, long.na)
	Xp.long.na <- apply(Xp.long.na, 1, any)
	Xgam.long.na <- apply(Xgam, 2, long.na)			# Xgam[rep(1:M, each = J),]
	Xgam.long.na <- apply(Xgam.long.na, 1, any)		# is.na(Xgam.long)
	Xom.long.na <- apply(Xom, 2, long.na)			# Xom[rep(1:M, each = J),]
	Xom.long.na <- apply(Xom.long.na, 1, any)		# is.na(Xom.long)
	
	y.long <- as.vector(t(getY(umf)))
	y.long.na <- is.na(y.long)
	
	covs.na <- apply(cbind(Xlam.long.na, Xgam.long.na, Xom.long.na, Xp.long.na), 
		1, any)
	
	## are any NA in covs not in y already?
	y.new.na <- covs.na & !y.long.na
	
	if(sum(y.new.na) > 0) {
		y.long[y.new.na] <- NA
		warning("Some observations have been discarded because correspoding covariates were missing.")
	}
	
	y <- matrix(y.long, M, J, byrow = TRUE)
	sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))
	#sites.to.remove <- sites.to.remove | plotArea.na
	
	num.to.remove <- sum(sites.to.remove)
	if(num.to.remove > 0) {
		y <- y[!sites.to.remove, ,drop = FALSE]
		Xlam <- Xlam[!sites.to.remove, ,drop = FALSE]
		Xgam <- Xgam[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
		Xom <- Xom[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
		Xp <- Xp[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
		plotArea <- plotArea[!sites.to.remove]
		delta <- delta[!sites.to.remove, ,drop =FALSE]
		warning(paste(num.to.remove,"sites have been discarded because of missing data."))
	}
	
	list(y = y, Xlam = Xlam, Xgam = Xgam, Xom = Xom, Xp = Xp, 
		plotArea = plotArea, delta = delta, 
		removed.sites = which(sites.to.remove))
}















pcountOpen <- function(lambdaformula, gammaformula, omegaformula, pformula, 
	data, mixture=c("P", "NB"), K, starts, method="BFGS", se=TRUE, ...)
{
mixture <- match.arg(mixture)
formlist <- list(lambdaformula=lambdaformula, gammaformula=gammaformula, 
	omegaformula=omegaformula, pformula=pformula)
D <- getDesign4(formlist, data)
y <- D$y; Xlam = D$Xlam; Xgam = D$Xgam; Xom = D$Xom; Xp = D$Xp
delta <- D$delta; plotArea <- D$plotArea
M <- nrow(y)
T <- ncol(y)
y <- matrix(y,  M, T)
if(missing(K)) K <- max(y, na.rm=T) + 20
if(K <= max(y, na.rm = TRUE))
	stop("specified K is too small. Try a value larger than any observation")
k <- 0:K
lk <- length(k)
y.kk <- apply(y, 2, rep, each=lk*lk)
delta.kk <- apply(delta, 2, rep, each=lk*lk)

# identical() returns FALSE b/c of environment differences
if(isTRUE(all.equal(gammaformula, ~1)) & isTRUE(all.equal(omegaformula, ~1)))
	are.null <- TRUE  	# Save time in likelihood evaluation
else are.null <- FALSE

lamParms <- colnames(Xlam)
gamParms <- colnames(Xgam)
omParms <- colnames(Xom)
detParms <- colnames(Xp)
nAP <- ncol(Xlam)
nGP <- ncol(Xgam)
nOP <- ncol(Xom)
nDP <- ncol(Xp)
nP <- nAP + nGP + nOP + nDP + ifelse(identical(mixture, "NB"), 1, 0)

cmin <- pmin(rep(k, times=lk), rep(k, each=lk))
transProbs <- function(x) #(N.itm1, N.it, om, gam, cmin)
	sum(dbinom(0:x[5], x[2], x[3]) * dpois(x[1]-(0:x[5]), x[4]))
mk.order <- matrix(1:(M*lk), M, lk)
mat.to.vec <- as.numeric(apply(mk.order, 1, rep, times=lk))

nll <- function(parms) { # No survey-specific NA handling.
	lambda <- exp(Xlam %*% parms[1 : nAP]) * plotArea
	p <- matrix(plogis(Xp %*% parms[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]), 
		M, T, byrow=TRUE)
	if(are.null) { 	# Save time by using scalars
		gamma <- exp(parms[(nAP+1) : (nAP+nGP)])
		omega <- plogis(parms[(nAP+nGP+1) : (nAP+nGP+nOP)])
		} 
	else { # could save time using vectors when only siteCovs are present
		gamma <- matrix(exp(Xgam %*% parms[(nAP+1) : (nAP+nGP)]), 
			M, T, byrow=TRUE)[,-T]
		omega <- matrix(plogis(Xom %*% parms[(nAP+nGP+1) : (nAP+nGP+nOP)]), 
			M, T, byrow=TRUE)[,-T]
		}
	L <- rep(NA, M)
	g1 <- sapply(k, function(x) dbinom(y[,1], x, p[,1])) 
	switch(mixture, 
		P = g2 <- sapply(k, function(x) dpois(x, lambda)),
		NB = g2 <- sapply(k, function(x) dnbinom(x, size=exp(parms[nP]), 
			mu=lambda)))
	g.star <- array(NA, c(M, lk, T-1))
	g3args <- cbind(rep(k, times=lk), rep(k, each=lk), 
		rep(omega, each=lk*lk), rep(gamma, each=lk*lk), cmin)	# recycle
	g3 <- apply(g3args, 1, transProbs) 	# Slow when covars on gamma/omega
										# Code in C?
	g3 <- array(g3, c(lk, lk, M, T-1))
	pT.kk <- rep(p[, T], each=lk*lk)
	g1.Tm1 <- dbinom(y.kk[,T], k, pT.kk) # recycle
	delta.Tkk <- delta.kk[,T-1]
	g3.Tm1 <- g3[,,, T-1]^delta.Tkk
	g.star[,, T-1] <- apply(g1.Tm1 * g3.Tm1, c(3,2), sum)	# recycle
	for(t in (T-1):2) {
		pt.kk <- rep(p[, t], each=lk*lk)
		g1.t <- dbinom(y.kk[,t], k, pt.kk)
		g.star.vec <- g.star[,, t][mat.to.vec]
		delta.tkk <- delta.kk[,t-1]
		g3.t <- g3[,,, t]^delta.tkk
		g.star[,, t-1] <- apply(g1.t * g3.t * g.star.vec, c(3,2), sum)
		}
	L <- rowSums(g1 * g2 * g.star[,, 1])
	-sum(log(L))
	} 
if(missing(starts))
	starts <- rep(0, nP)
fm <- optim(starts, nll, method=method, hessian=se, ...)
opt <- fm 
ests <- fm$par
if(identical(mixture,"NB")) nbParm <- "alpha"
	else nbParm <- character(0)
names(ests) <- c(lamParms, gamParms, omParms, detParms, nbParm)
if(se) {           
	tryCatch(covMat <- solve(fm$hessian),
		error=function(x) simpleError("Hessian is not invertible.  Try using fewer covariates."))
} else covMat <- matrix(NA, nP, nP)

fmAIC <- 2 * fm$value + 2 * nP

lamName <- ifelse(all(data@plotArea == 1), "Abundance", "Density")
lamEstimates <- unmarkedEstimate(name = lamName, short.name = "lam",
	estimates = ests[1:nAP], covMat = as.matrix(covMat[1:nAP,1:nAP]), 
	invlink = "exp", invlinkGrad = "exp")
gamEstimates <- unmarkedEstimate(name = "Recruitment", short.name = "gam",
	estimates = ests[(nAP+1) : (nAP+nGP)],
	covMat = as.matrix(covMat[(nAP+1) : (nAP+nGP), (nAP+1) : (nAP+nGP)]),
	invlink = "exp", invlinkGrad = "exp")
omEstimates <- unmarkedEstimate(name = "Survival", short.name = "omega",
	estimates = ests[(nAP+nGP+1) : (nAP+nGP+nOP)],
	covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP), 
		(nAP+nGP+1) : (nAP+nGP+nOP)]),
	invlink = "logistic", invlinkGrad = "logistic.grad")	
detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
	estimates = ests[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)], 
	covMat = as.matrix(covMat[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP), 
		(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]), 
		invlink = "logistic", invlinkGrad = "logistic.grad")
estimateList <- unmarked:::unmarkedEstimateList(list(lambda=lamEstimates,
	gamma = gamEstimates, omega = omEstimates, det=detEstimates))
if(identical(mixture, "NB")) {
	estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
		short.name = "alpha", estimates = ests[nP],
		covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
		invlinkGrad = "exp")
	}
umfit <- new("unmarkedFitPCount", fitType = "pcountOpen", call = match.call(),
	formula = as.formula(paste(unlist(formlist), collapse=" ")), 
	data = data, sitesRemoved=D$removed.sites,
	estimates = estimateList, AIC = fmAIC, opt = opt, negLogLike = fm$value,
	nllFun = nll, K = K, mixture = mixture)
return(umfit)
}


                                                       

## # Simulate	
## set.seed(1)
## M <- 50
## T <- 5
## veght <- rnorm(M)
## isolation <- matrix(rnorm(M*T), M, T)
## date <- matrix(rnorm(M*T, 1), M, T)
## lambda <- exp(-1 + 0.5*veght)
## y <- p <- N <- gamma <- matrix(NA, M, T)
## S <- G <- matrix(NA, M, T-1)
## gamma[] <- exp(-1 + -1*isolation)
## N[,1] <- rpois(M, lambda)
## for(t in 1:(T-1)) {
## 	S[,t] <- rbinom(M, N[,t], 0.8)
## 	G[,t] <- rpois(M, gamma[,t])
## 	N[,t+1] <- S[,t] + G[,t]
## 	}
## p[] <- plogis(-1 + 1*date)
## y[] <- rbinom(M*T, N, p)
## #y[1, 1:2] <- NA
## #isolation[1, 1:2] <- NA
## #date[1, 1:2] <- NA
## dat <- data.frame(veght)



## # Prepare data                               
## umf <- unmarkedFramePCountOpen(y = y, siteCovs = dat, 
## 	obsCovs = list(isolation=isolation, date=date))

## umf
## #plot(umf)
## summary(umf)


## # Fit some models
## system.time(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10))
## backTransform(m1, "lambda")
## backTransform(m1, "gamma")
## backTransform(m1, "omega")
## backTransform(m1, "det")

## # Real slow
## system.time(m2 <- pcountOpen(~veght, ~isolation, ~1, ~date, umf, K=10, se=F,
## 	control=list(maxit=20, trace=T), starts=c(-1, 0.5, -1, -1, 1.5, -1, 1)))

## #debugonce(pcountOpen)












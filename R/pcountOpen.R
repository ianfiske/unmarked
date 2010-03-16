


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










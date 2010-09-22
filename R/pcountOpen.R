


pcountOpen <- function(lambdaformula, gammaformula, omegaformula, pformula,
    data, mixture=c("P", "NB"), K, dynamics=c("constant", "autoreg", "notrend"), 
    fix=c("none", "gamma", "omega"), 
    starts, method="BFGS", se=TRUE, ...)
{
mixture <- match.arg(mixture)
dynamics <- match.arg(dynamics)
fix <- match.arg(fix)
formlist <- list(lambdaformula=lambdaformula, gammaformula=gammaformula,
    omegaformula=omegaformula, pformula=pformula)
formula <- as.formula(paste(unlist(formlist), collapse=" "))
D <- unmarked:::getDesign(data, formula)
y <- D$y; Xlam <- D$Xlam; Xgam <- D$Xgam; Xom <- D$Xom; Xp <- D$Xp
delta <- D$delta
M <- nrow(y)
T <- ncol(y)
y <- matrix(y, M, T)
if(missing(K)) K <- max(y, na.rm=T) + 20
if(K <= max(y, na.rm = TRUE))
    stop("specified K is too small. Try a value larger than any observation")
k <- 0:K
lk <- length(k)
k.times <- rep(k, times=lk)
k.each <- rep(k, each=lk)

lamParms <- colnames(Xlam)
gamParms <- colnames(Xgam)
omParms <- colnames(Xom)
detParms <- colnames(Xp)
nAP <- ncol(Xlam)
nGP <- ncol(Xgam)
nOP <- ncol(Xom)
nDP <- ncol(Xp)
if(identical(fix, "gamma")) {
    if(!identical(dynamics, "constant")) 
        stop("dynamics must be constant when fixing gamma or omega")
    if(nGP > 1) stop("gamma covariates not allowed when fix==gamma")
    else { nGP <- 0; gamParms <- character(0) }
    }
if(identical(fix, "omega")) {
    if(!identical(dynamics, "constant")) 
        stop("dynamics must be constant when fixing gamma or omega")    
    if(nOP > 1) stop("omega covariates not allowed when fix==omega")
    else { nOP <- 0; omParms <- character(0) }
    }
nP <- nAP + nGP + nOP + nDP + ifelse(identical(mixture, "NB"), 1, 0)

nll <- function(parms) {
    lambda <- drop(exp(Xlam %*% parms[1 : nAP]))
    p <- matrix(plogis(Xp %*% parms[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]),
                M, T, byrow=TRUE)
    if(identical(fix, "omega"))
        omega <- matrix(1, M, T-1)
    else
        omega <- matrix(plogis(Xom %*% parms[(nAP+nGP+1) : (nAP+nGP+nOP)]),
            M, T, byrow=TRUE)[,-T] ^ delta
    if(dynamics == "notrend")
        gamma <- (1-omega)*lambda
    else {
        if(identical(fix, "gamma")) 
            gamma <- matrix(0, M, T-1)
        else 
            gamma <- matrix(drop(exp(Xgam %*% parms[(nAP+1) : (nAP+nGP)])),
                M, T, byrow=TRUE)[,-T] * delta
        }
    L <- numeric(M)
    g.star <- matrix(NA, lk, T-1)
    for(i in 1:M) {
        g1 <- dbinom(y[i,1], k, p[i,1])
        switch(mixture, 
            P = g2 <- dpois(k, lambda[i]),
            NB = g2 <- dnbinom(k, size=exp(parms[nP]), mu=lambda[i]))    
        g3args <- cbind(k.times, k.each, 
            rep(omega[i,], each=lk*lk), 
            rep(gamma[i,], each=lk*lk)) # recycle
        if(dynamics == "autoreg" & fix != "gamma")
            g3args[,4] <- g3args[,4] * g3args[,2]
        convMat <- matrix(0, nrow(g3args), K+1)
        for(z in k) {
            # this could be pulled out of likelihood
            nonzero <- which(g3args[,2] >= z & g3args[,1] >= z) 
            convMat[nonzero, z+1] <- dbinom(z, g3args[nonzero,2], 
                g3args[nonzero,3]) * dpois(g3args[nonzero,1] - z, 
                g3args[nonzero,4])
            }
        g3 <- rowSums(convMat)
        g3 <- array(g3, c(lk, lk, T-1))
        p.Tk <- rep(p[i, T], each=lk)
        p.Tk <- rep(p[i, T], each=lk)
        y.Tk <- rep(y[i, T], each=lk)
        g1.T <- dbinom(y.Tk, k, p.Tk)
        g3.T <- g3[,, T-1]
        g.star[, T-1] <- colSums(g1.T * g3.T)	# or rowSums?
        # NA handling: this will properly determine last obs for each site
        g.star[,T-1][is.na(g.star[, T-1])] <- 1
        for(t in (T-1):2) {
            p.tk <- rep(p[i, t], each=lk)
            y.tk <- rep(y[i, t], each=lk)
            g1.t <- dbinom(y.tk, k, p.tk)
            g.star[, t-1] <- colSums(g1.t * g3[,,t-1] * g.star[,t])
            g.star.na <- is.na(g.star[, t-1])
            g.star[,t-1][g.star.na] <- g.star[,t][g.star.na]
            }
        L[i] <- sum(g1 * g2 * g.star[,1])
        }
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
	covMat <- tryCatch(solve(fm$hessian), error=function(x) 
		simpleError("Hessian is not invertible. Try using fewer covariates or providing starting values."))
	if(class(covMat)[1] == "simpleError") {
		print(covMat$message)
		covMat <- matrix(NA, nP, nP)
		} 
    } else covMat <- matrix(NA, nP, nP)

fmAIC <- 2*fm$value + 2*nP

lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lam",
    estimates = ests[1:nAP], covMat = as.matrix(covMat[1:nAP,1:nAP]),
    invlink = "exp", invlinkGrad = "exp")
detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
    estimates = ests[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)],
    covMat = as.matrix(covMat[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP),
        (nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]),
        invlink = "logistic", invlinkGrad = "logistic.grad")
estimateList <- unmarked:::unmarkedEstimateList(list(lambda=lamEstimates))
if(!identical(fix, "gamma")) 
    estimateList@estimates$gamma <- unmarkedEstimate(name = "Recruitment", 
        short.name = "gam", estimates = ests[(nAP+1) : (nAP+nGP)],
        covMat = as.matrix(covMat[(nAP+1) : (nAP+nGP), (nAP+1) : (nAP+nGP)]),
        invlink = "exp", invlinkGrad = "exp")
if(!identical(fix, "omega")) 
    estimateList@estimates$omega <- unmarkedEstimate(name = "Survival", 
        short.name = "omega", estimates = ests[(nAP+nGP+1) : (nAP+nGP+nOP)],
        covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
            (nAP+nGP+1) : (nAP+nGP+nOP)]),
        invlink = "logistic", invlinkGrad = "logistic.grad")   
estimateList@estimates$det <- detEstimates
if(identical(mixture, "NB")) {
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
        invlinkGrad = "exp")
    }
umfit <- new("unmarkedFitPCountOpen", fitType = "pcountOpen", call = match.call(),
    formula = formula,
    formlist = formlist, data = data, sitesRemoved=D$removed.sites,
    estimates = estimateList, AIC = fmAIC, opt = opt, negLogLike = fm$value,
    nllFun = nll, K = K, mixture = mixture)
return(umfit)
}




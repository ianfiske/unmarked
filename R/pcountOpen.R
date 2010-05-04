


pcountOpen <- function(lambdaformula, gammaformula, omegaformula, pformula,
    data, mixture=c("P", "NB"), K, fix=c("none", "gamma", "omega"), 
    starts, method="BFGS", se=TRUE, ...)
{
mixture <- match.arg(mixture)
fix <- match.arg(fix)
formlist <- list(lambdaformula=lambdaformula, gammaformula=gammaformula,
    omegaformula=omegaformula, pformula=pformula)
formula <- as.formula(paste(unlist(formlist), collapse=" "))
D <- unmarked:::getDesign(data, formula)
y <- D$y; Xlam <- D$Xlam; Xgam <- D$Xgam; Xom <- D$Xom; Xp <- D$Xp
delta <- D$delta
M <- nrow(y)
T <- ncol(y)
y <- matrix(y,  M, T)
yind <- matrix(1:(M*T), M, T)
first <- 1:M
for(i in 1:M) 
    if(any(is.na(y[i,])))
        first[i] <- yind[i, min(which(!is.na(y[i,])))]
if(missing(K)) K <- max(y, na.rm=T) + 20
if(K <= max(y, na.rm = TRUE))
    stop("specified K is too small. Try a value larger than any observation")
k <- 0:K
lk <- length(k)
y.kk <- apply(y, 2, rep, each=lk*lk)
if(all(delta==1)) delta <- 1
    # delta.kk <- 1 else delta.kk <- apply(delta, 2, rep, each=lk*lk)

lamParms <- colnames(Xlam)
gamParms <- colnames(Xgam)
omParms <- colnames(Xom)
detParms <- colnames(Xp)
nAP <- ncol(Xlam)
nGP <- ncol(Xgam)
nOP <- ncol(Xom)
nDP <- ncol(Xp)
if(identical(fix, "gamma")) {
    if(nGP > 1) stop("gamma covariates not allowed when fix==gamma")
    else { nGP <- 0; gamParms <- character(0) }
    }
if(identical(fix, "omega")) {
    if(nOP > 1) stop("omega covariates not allowed when fix==omega")
    else { nOP <- 0; omParms <- character(0) }
    }

nP <- nAP + nGP + nOP + nDP + ifelse(identical(mixture, "NB"), 1, 0)

# Save time in likelihood evaluation
# identical() returns FALSE b/c of environment differences
equal.ints <- identical(length(table(delta)), 1L)
if(isTRUE(all.equal(gammaformula, ~1)) & isTRUE(all.equal(omegaformula, ~1)) & 
    equal.ints)
    goDims <- "scalar"
    else {
        goParms <- unique(c(all.vars(gammaformula), all.vars(omegaformula)))
        if(!any(goParms %in% colnames(obsCovs(data))) & equal.ints)
            goDims <- "vector"
            else goDims <- "matrix"
        }
mk.order <- matrix(1:(M*lk), M, lk)
mat.to.vec <- as.numeric(apply(mk.order, 1, rep, times=lk))
g.star <- array(NA, c(M, lk, T-1))

nll <- function(parms) { # No survey-specific NA handling.
    lambda <- exp(Xlam %*% parms[1 : nAP])
    p <- matrix(plogis(Xp %*% parms[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]),
        M, T, byrow=TRUE)
    switch(goDims,
    scalar = {
        gamma <- drop(exp(parms[(nAP+1) : (nAP+nGP)]))
        omega <- drop(plogis(parms[(nAP+nGP+1) : (nAP+nGP+nOP)]))
        },
    vector = {
        gamma <- matrix(drop(exp(Xgam %*% parms[(nAP+1) : (nAP+nGP)])), 
            M, T, byrow=TRUE)[,1]
        omega <- matrix(drop(plogis(Xom %*% parms[(nAP+nGP+1) : (nAP+nGP+nOP)])), 
            M, T, byrow=TRUE)[,1]
        },
    matrix = {
        gamma <- matrix(drop(exp(Xgam %*% parms[(nAP+1) : (nAP+nGP)])),
            M, T, byrow=TRUE)[,-T] * delta
        omega <- matrix(drop(plogis(Xom %*% parms[(nAP+nGP+1) : (nAP+nGP+nOP)])),
            M, T, byrow=TRUE)[,-T] ^ delta
        })
    if(identical(fix, "gamma")) gamma[] <- 0
        else if(identical(fix, "omega")) omega[] <- 1    
    #g1 <- sapply(k, function(x) dbinom(y[first], x, p[first]))
    g1 <- sapply(k, function(x) dbinom(y[,1], x, p[,1]))    
    g1[is.na(g1)] <- 1
    switch(mixture,
        P = g2 <- sapply(k, function(x) dpois(x, lambda)),
        NB = g2 <- sapply(k, function(x) dnbinom(x, size=exp(parms[nP]),
           mu=lambda)))
    g2[is.na(g2)] <- 1
    g3args <- cbind(rep(k, times=lk), rep(k, each=lk), 
        rep(omega, each=lk*lk), #^delta.kk
        rep(gamma, each=lk*lk)) #*delta.kk)# recycle
    convMat <- matrix(NA, nrow(g3args), K+1)
    for(i in k)
        convMat[,i+1] <- dbinom(i, g3args[,2], g3args[,3]) * 
            dpois(g3args[,1] - i, g3args[,4])
    g3 <- rowSums(convMat)#^delta.kk
    g3 <- array(g3, c(lk, lk, M, T-1))
    pT.kk <- rep(p[, T], each=lk*lk)
    g1.T <- dbinom(y.kk[,T], k, pT.kk) # recycle
    #delta.Tkk <- delta.kk[,T-1]
    g3.T <- g3[,,, T-1]#^delta.Tkk
    g.star[,, T-1] <- apply(g1.T * g3.T, 2, colSums)	# recycle
    # NA handling: this will properly determine last obs for each site
    g.star[,,T-1][is.na(g.star[,, T-1])] <- 1
    for(t in (T-1):2) {
        pt.kk <- rep(p[, t], each=lk*lk)
        g1.t <- dbinom(y.kk[,t], k, pt.kk)
        g.star.vec <- g.star[,, t][mat.to.vec]
        g3.t <- g3[,,, t-1]#^delta.tkk
        #delta.tkk <- delta.kk[,t-1]
        g.star[,, t-1] <- apply(g1.t * g3.t * g.star.vec, 2, colSums)
        g.star[,,t-1][is.na(g.star[,, t-1])] <- g.star[,,t][is.na(g.star[,, t-1])]
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







pcountOpen <- function(lambdaformula, gammaformula, omegaformula, pformula,
    data, mixture=c("P", "NB"), K, 
    dynamics=c("constant", "autoreg", "notrend"), 
    fix=c("none", "gamma", "omega"), 
    starts, method="BFGS", se=TRUE, ...)
{
mixture <- match.arg(mixture)
dynamics <- match.arg(dynamics)
fix <- match.arg(fix)
if(identical(dynamics, "notrend") & !identical(lambdaformula, omegaformula))
    stop("lambdaformula and omegaformula must be identical for notrend model") 
formlist <- list(lambdaformula=lambdaformula, gammaformula=gammaformula,
    omegaformula=omegaformula, pformula=pformula)
formula <- as.formula(paste(unlist(formlist), collapse=" "))
D <- unmarked:::getDesign(data, formula)
y <- D$y; Xlam <- D$Xlam; Xgam <- D$Xgam; Xom <- D$Xom; Xp <- D$Xp
delta <- D$delta; go.dims <- D$go.dims
deltamax <- max(delta, na.rm=TRUE)
M <- nrow(y)
T <- ncol(y)
y <- matrix(y, M, T)
if(missing(K)) K <- max(y, na.rm=T) + 20
if(K <= max(y, na.rm = TRUE))
    stop("specified K is too small. Try a value larger than any observation")
k <- 0:K
lk <- length(k)

first <- apply(y, 1, function(x) min(which(!is.na(x))))
last  <- apply(y, 1, function(x) max(which(!is.na(x)))) 

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
else 
    if(identical(dynamics, "notrend")) {
        if(nGP > 1) stop("gamma covariates not allowed when dyamics==notrend")
        else { nGP <- 0; gamParms <- character(0) }
        }
if(identical(fix, "omega")) {
    if(!identical(dynamics, "constant")) 
        stop("dynamics must be constant when fixing gamma or omega")    
    if(nOP > 1) stop("omega covariates not allowed when fix==omega")
    else { nOP <- 0; omParms <- character(0) }
    }
nP <- nAP + nGP + nOP + nDP + ifelse(identical(mixture, "NB"), 1, 0)
if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP)) 

nll <- function(parms) {
    lambda <- drop(exp(Xlam %*% parms[1 : nAP]))
    p <- matrix(plogis(Xp %*% parms[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]),
        M, T, byrow=TRUE)
    if(identical(fix, "omega"))
        omega <- matrix(1, M, T-1)
    else
        omega <- matrix(plogis(Xom %*% parms[(nAP+nGP+1) : (nAP+nGP+nOP)]),
            M, T-1, byrow=TRUE)
    if(dynamics == "notrend")
        gamma <- (1-omega)*lambda
    else {
        if(identical(fix, "gamma")) 
            gamma <- matrix(0, M, T-1)
        else 
            gamma <- matrix(exp(Xgam %*% parms[(nAP+1) : (nAP+nGP)]),
                M, T-1, byrow=TRUE)
        }
    if(identical(go.dims, "scalar")) {
        g3 <- array(NA, c(lk, lk, deltamax))
        g3[,,1] <- tranProbs(k, omega[1,1], gamma[1,1], 1, dynamics)
        if(deltamax > 1) {
            for(d in 2:deltamax)
                g3[,,d] <- g3[,,d-1] %*% g3[,,d-1]
            }
        }
    L <- numeric(M)
    for(i in 1:M) {
        first.i <- first[i]
        last.i <- last[i]
        Nsub1 <- k >= y[i, first.i]
        N1 <- k[Nsub1]
        g1 <- dbinom(y[i, first.i], N1, p[i, first.i])
        switch(mixture, 
            P = g2 <- dpois(N1, lambda[i]),
            NB = g2 <- dnbinom(N1, size=exp(parms[nP]), mu=lambda[i]))    
        if(first.i == last.i & first.i == 1) {
            L[i] <- sum(g1 * g2)
            next
            }
        NsubT <- k >= y[i, last.i]
        NT <- k[NsubT]
        g.star <- matrix(NA, lk, last.i-1) # must have lk rows
        if(identical(go.dims, "scalar"))
            g3.T <- g3[NsubT,, delta[i, last.i]]
        else
            g3.T <- tranProbs(k, omega[i, last.i-1], gamma[i, last.i-1], 
                delta[i, last.i], dynamics)[NsubT,] # not delta[i, last.i-1]
        g1.T <- dbinom(y[i, last.i], NT, p[i, last.i])
        g.star[, last.i-1] <- colSums(g1.T * g3.T)
        if(first.i == last.i & first.i > 1) {
            L[i] <- sum(g2 * colSums(g1 * g3.T * g.star[NsubT,last.i-1])[NsubT])
            next
            }
        if((last.i - first.i) > 1) { 
            for(t in (last.i-1):(first.i+1)) {
                if(is.na(y[i, t])) # time gap dealt with by delta
                    g.star[,t-1] <- g.star[,t]
                else {
                    Nsub <- k >= y[i, t]
                    N <- k[Nsub]
                    g1.t <- dbinom(y[i, t], N, p[i, t])
                    if(identical(go.dims, "scalar"))
                        g3.t <- g3[,,delta[i, t]]
                    else
                        g3.t <- tranProbs(k, omega[i, t-1], gamma[i, t-1], 
                            delta[i, t], dynamics)
                    g.star[,t-1] <- colSums(g1.t * g3.t[Nsub,] * g.star[Nsub,t])
                    }
                }
            }
        if(first.i == 1)
            L[i] <- sum(g1 * g2 * g.star[Nsub1, first.i])
        else {
            # g3.t should use delta[i, first.i]
            if(identical(go.dims, "scalar"))
                g3.1 <- g3[,,delta[i, first.i]]
            else
                g3.1 <- tranProbs(k, omega[i, first.i], gamma[i, first.i], 
                    delta[i, first.i], dynamics)
            L[i] <- sum(g2 * colSums(g1 * g3.1[Nsub1,] * 
                g.star[Nsub1,first.i])[Nsub1])
            }
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
if(!(identical(fix, "gamma") | identical(dynamics, "notrend"))) 
    estimateList@estimates$gamma <- unmarkedEstimate(name = "Recruitment", 
        short.name = "gam", estimates = ests[(nAP+1) : (nAP+nGP)],
        covMat = as.matrix(covMat[(nAP+1) : (nAP+nGP), (nAP+1) : (nAP+nGP)]),
        invlink = "exp", invlinkGrad = "exp")
if(!identical(fix, "omega")) 
    estimateList@estimates$omega <- unmarkedEstimate(name="Apparent Survival", 
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
umfit <- new("unmarkedFitPCountOpen", fitType = "pcountOpen", 
    call = match.call(), formula = formula, formlist = formlist, data = data, 
    sitesRemoved=D$removed.sites, estimates = estimateList, AIC = fmAIC, 
    opt = opt, negLogLike = fm$value, nllFun = nll, K = K, mixture = mixture, 
    dynamics = dynamics)
return(umfit)
}




   
    




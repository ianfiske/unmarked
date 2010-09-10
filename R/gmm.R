
# data will need to be an unmarkedMultFrame
gmm <- function(lambdaformula, phiformula, pformula, data, mixture=c('P', 'NB'),
    K, starts, method = "BFGS", control = list(), se = TRUE)
{
if(!is(data, "unmarkedFrameGMM"))
    stop("Data is not of class unmarkedFrameGMM.")

mixture <- match.arg(mixture)

formlist <- list(lambdaformula = lambdaformula, phiformula = phiformula, 
    pformula = pformula)
form <- as.formula(paste(unlist(formlist), collapse=" "))
D <- getDesign(data, formula = form)

Xlam <- D$Xlam
Xphi <- D$Xphi 
Xdet <- D$Xdet
y <- D$y  # MxJT 


Xlam.offset <- D$X.offset
Xphi.offset <- D$Xphi.offset
Xdet.offset <- D$Xdet.offset
if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

if(missing(K) || is.null(K)) K <- max(y, na.rm=TRUE) + 20
k <- 0:K
lk <- length(k)
M <- nrow(y)  
T <- data@numPrimary
R <- ncol(y)
J <- R / T

y <- array(y, c(M, T, J))
yt <- apply(y, 1:2, sum)

piFun <- data@piFun

lamPars <- colnames(Xlam)
phiPars <- colnames(Xphi)
detPars <- colnames(Xdet)
nLP <- ncol(Xlam)
nPP <- ncol(Xphi)
nDP <- ncol(Xdet)
nP <- nLP + nPP + nDP + ifelse(mixture=='NB', 1, 0)

p <- array(as.numeric(NA), c(M, T, J))
cp <- array(as.numeric(NA), c(M, T, J+1))
A <- matrix(0, lk, T)
g <- matrix(as.numeric(NA), M, lk)

lfac.k <- lgamma(k+1)
kmn <- array(NA, c(M, T, lk))
lfac.kmn <- array(0, c(M, T, lk))
for(i in 1:M) {
    for(t in 1:T) {
        kmn[i,t,] <- k - yt[i,t]
        zp <- kmn[i,t,] >= 0
        lfac.kmn[i, t, zp] <- lgamma(kmn[i,t,zp] + 1)
        }
    }

nll <- function(pars) {
    lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset) 
    phi <- drop(plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset))
    p[] <- plogis(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)] + Xdet.offset)
    for(i in 1:T) cp[,i,1:J] <- do.call(piFun, list(p[,i,]))
    cp[,,1:J] <- cp[,,1:J] * phi
    cp[,,J+1] <- apply(cp[,,1:J], 1:2, sum) 
    switch(mixture, 
        P = f <- sapply(k, function(x) dpois(x, lambda)),
        NB = f <- sapply(k, dnbinom(x, mu=lambda, size=pars[nP]))
        )
    for(i in 1:M) {
        for(t in 1:T)
            A[,t] <- lfac.k - lfac.kmn[i,t,] + sum(y[i,t,]*log(cp[i,t,1:J])) + 
                kmn[i,t,]*log(cp[i,t,J+1])
        g[i,] <- exp(rowSums(A))
        }
    ll <- rowSums(f*g)
    -sum(log(ll))
    }
    
if(missing(starts)) starts <- rep(0, nP)
fm <- optim(starts, nll, method = method, hessian = se, control = control)
opt <- fm
if(se) {
		covMat <- tryCatch(solve(fm$hessian), error=function(x) 
        simpleError("Hessian is singular. Try using fewer covariates."))
    if(identical(class(covMat)[1], "simpleError")) {
        warning(covMat$message)
        covMat <- matrix(NA, nP, nP)
        }
    } else covMat <- matrix(NA, nP, nP)
ests <- fm$par
fmAIC <- 2 * fm$value + 2 * nP
names(ests) <- c(lamPars, phiPars, detPars)

lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
    estimates = ests[1:nLP],
    covMat = as.matrix(covMat[1:nLP, 1:nLP]), invlink = "exp",
    invlinkGrad = "exp")
phiEstimates <- unmarkedEstimate(name = "Availability", short.name = "phi",
    estimates = ests[(nLP+1):(nLP+nPP)],
    covMat = as.matrix(covMat[(nLP+1):(nLP+nPP)]), invlink = "logistic",
    invlinkGrad = "logistic.grad")
detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
    estimates = ests[(nLP+nPP+1):(nLP+nPP+nDP)],
    covMat = as.matrix(
        covMat[(nLP+nPP+1):(nLP+nPP+nDP), (nLP+nPP+1):(nLP+nPP+nDP)]), 
    invlink = "logistic", invlinkGrad = "logistic.grad")
estimateList <- unmarkedEstimateList(list(lambda=lamEstimates,
    phi=phiEstimates, det=detEstimates))

umfit <- new("unmarkedFitGMM", fitType = "gmn", 
    call = match.call(), formula = form, formlist = formlist,    
    data = data, estimates = estimateList, sitesRemoved = D$removed.sites, 
    AIC = fmAIC, opt = opt, negLogLike = fm$value, nllFun = nll)

return(umfit)
}





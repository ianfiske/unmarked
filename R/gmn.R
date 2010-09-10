
# data will need to be an unmarkedMultFrame
gmn <- function(lambdaformula, phiformula, pformula, data, mixture=c('P', 'NB'),
    K, starts, method = "BFGS", control = list(), se = TRUE)
{
if(!is(data, "unmarkedMultFrame"))
    stop("Data is not an unmarkedMultFrame.")

mixture <- match.arg(mixture)
k <- 0:K
lk <- length(k)

formula <- list(lambdaformula = lambdaformula, phiformula = phiformula, 
    pformula = pformula)
D <- getDesign(data, formula = as.formula(paste(unlist(formula), collapse=" ")))
Xlam <- D$Xlam
Xphi <- D$Xphi 
Xdet <- D$Xdet
y.mto <- D$y  # MxJxL
n.mt <- apply(y, 1:2, sum)

Xlam.offset <- D$X.offset
Xphi.offset <- D$Xphi.offset
Xdet.offset <- D$Xdet.offset
if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

M <- nrow(y)  
T <- ncol(y)
O <- dim(y)[3]
R <- obsNum(data)

piFun <- data@piFun

lamPars <- colnames(Xlam)
phiPars <- colnames(Xphi)
detPars <- colnames(Xdet)
nLP <- ncol(Xlam)
nPP <- ncol(Xphi)
nDP <- ncol(Xdet)
nP <- nLP + nPP + nDP + ifelse(mixture=='NB', 1, 0)

p <- array(as.numeric(NA), c(M, J, L))
cp <- array(as.numeric(NA), c(M, J, R+1))   # L does not necessarily equal R
A <- matrix(0, M, lk)
g <- matrix(as.numeric(NA), M, lk)

lfac.k <- lgamma(k+1)
kmn <- array(NA, c(M, T, lk))
lfac.kmn <- array(0, c(M, T, lk))
for(i in 1:M) {
    for(t in 1:T) {
        kmn[i,t,] <- k-n[i,t]
        zp <- kmn >= 0
        lfac.kmn[i,t,zp] <- lgamma(kmn[zp]+1)
        }
    }

nll <- function(pars) {
		lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset) 
    phi <- plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset)
    p[] <- plogis(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)] + Xdet.offset)
    cp.mat <- apply(p, 2, function(x) do.call(piFun, list(x))) * phi
    cp[] <- cbind(cp.mat, 1-rowSums(pi.mat))
    switch(mixture, 
        P = f <- sapply(k, function(x) dpois(x, lambda)),
        NB = f <- sapply(k, dnbinom(x, mu=lambda, size=pars[nP]))
        )
    for(i in 1:M) {
        for(t in 1:T)
            A[,T] <- lfac.k - lfac.kmn[i,t,] + sum(y[i,j,]*log(cp[i,j,1:R])) + 
                kmn[i,t,]*log(cp[i,j,R+1])
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
    estimates = ests[1:nAP],
		covMat = as.matrix(covMat[1:nAP, 1:nAP]), invlink = "exp",
		invlinkGrad = "exp")
phiEstimates <- unmarkedEstimate(name = "Availability", short.name = "phi",
    estimates = ests[(nLP+1):(nLP+nPP)],
		covMat = as.matrix(covMat[(nLP+1):(nLP+nPP)]), invlink = "logistic",
		invlinkGrad = "logistic.grad")
detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
		estimates = ests[(nPP+1):(nLP+nPP+nDP)],
		covMat = as.matrix(covMat[(nPP+1):(nLP+nPP+nDP), (nPP+1):(nLP+nPP+nDP)]), 
		invlink = "logistic", invlinkGrad = "logistic.grad")
estimateList <- unmarkedEstimateList(list(lambda=stateEstimates,
		phi=phiEstimates, det=detEstimates))

umfit <- new("unmarkedFitGMN", fitType = "gmn", 
		call = match.call(), formula = formula, data = data, 
		estimates = estimateList, sitesRemoved = D$removed.sites, 
		AIC = fmAIC, opt = opt, negLogLike = fm$value, nllFun = nll)

return(umfit)
}





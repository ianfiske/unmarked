


distsampOpen <- function(lambdaformula, gammaformula, omegaformula, sigmaformula,    data,
    keyfun=c("halfnorm", "exp", "hazard", "uniform"),
    output=c("abund", "density"), unitsOut=c("ha", "kmsq"),
                         mixture=c("P", "NB", "ZIP"), K,
    dynamics=c("constant", "autoreg", "notrend", "trend"),
    fix=c("none", "gamma", "omega"),
    iotaformula = ~1,
    starts, method="BFGS", se=TRUE, immigration=FALSE, nintervals=10, ...)
{

### This block below andy put here from gdistsamp
if(!is(data, "unmarkedFrameDSO"))
    stop("Data is not of class unmarkedFrameDSO.")

keyfun <- match.arg(keyfun)
if(!keyfun %in% c("halfnorm", "exp", "hazard", "uniform"))
    stop("keyfun must be 'halfnorm', 'exp', 'hazard', or 'uniform'")
output <- match.arg(output)
unitsOut <- match.arg(unitsOut)
db <- data@dist.breaks
w <- diff(db)
tlength <- data@tlength
survey <- data@survey
unitsIn <- data@unitsIn
#####


mixture <- match.arg(mixture)
dynamics <- match.arg(dynamics)
fix <- match.arg(fix)
if(identical(dynamics, "notrend") &
   !identical(lambdaformula, omegaformula))
    stop("lambdaformula and omegaformula must be identical for notrend model")
formlist <- list(lambdaformula=lambdaformula, gammaformula=gammaformula,
    omegaformula=omegaformula, sigmaformula=sigmaformula, iotaformula=iotaformula)

formula <- as.formula(paste(unlist(formlist), collapse=" "))
D <- unmarked:::getDesign(data, formula)
y <- D$y

Xlam <- D$Xlam
Xgam <- D$Xgam
Xom <- D$Xom
Xsig <- D$Xp # Wrong length. Need new getDesign Function
Xiota<- D$Xiota

delta <- D$delta; go.dims <- D$go.dims
deltamax <- max(delta, na.rm=TRUE)
M <- nrow(y)
T <- data@numPrimary
J <- ncol(getY(data)) / T

# FIXME: temporary fix until new getDesign is ready
Xsig <- Xsig[1:(M*T),,drop=FALSE]

Xlam.offset <- D$Xlam.offset
Xgam.offset <- D$Xgam.offset
Xom.offset <- D$Xom.offset
Xsig.offset <- D$Xp.offset   # Wrong length
Xiota.offset<- D$Xiota.offset

if(is.null(Xlam.offset)) Xlam.offset <- rep(0, M)
if(is.null(Xgam.offset)) Xgam.offset <- rep(0, M*(T-1))
if(is.null(Xom.offset)) Xom.offset <- rep(0, M*(T-1))
if(is.null(Xsig.offset)) Xsig.offset <- rep(0, M*T)
if(is.null(Xiota.offset)) Xiota.offset<- rep(0, M*(T-1))

yna <- is.na(y)
yna[] <- as.integer(yna)
y <- array(y, c(M, J, T))
yt <- apply(y, c(1,3), function(x) {
    if(all(is.na(x)))
        return(NA)
    else return(sum(x, na.rm=TRUE))
    })
ytna <- apply(is.na(y), c(1,3), all)
ytna <- matrix(ytna, nrow=M)
ytna[] <- as.integer(ytna)

first <- apply(!ytna, 1, function(x) min(which(x)))
last  <- apply(!ytna, 1, function(x) max(which(x)))

if(missing(K)) {
    K <- max(y, na.rm=T) + 20
    warning("K was not specified and was set to ", K, ".")
}
if(K <= max(y, na.rm = TRUE))
    stop("specified K is too small. Try a value larger than any observation")
k <- 0:K
lk <- length(k)

# FIXME: These should be arguments
#db <- c(0, 25, 50, 75, 100)
#tlength <- 1

if(length(db)-1 != J)
    stop("duh")
####a <- numeric(J)
a<- u<- matrix(NA, M, J)

if(survey=="point"){
                for(i in 1:M) {
                a[i, 1] <- pi*db[2]^2
                for(j in 2:J)
                    a[i, j] <- pi*db[j+1]^2 - sum(a[i, 1:(j-1)])
                u[i,] <- a[i,] / sum(a[i,])
            }
            }
if(survey=="line"){
        for(i in 1:M) {
            a[i,] <- tlength[i] * w
            u[i,] <- a[i,] / sum(a[i,])
            }
    }
#   a[1] <- tlength*db[2]   # Note: should be vector valued.... a and u are passed to cpp routine
  # They both seem to be (i,j) matrices within the .cpp function
#for(j in 2:J)
#    a[j] <- tlength*db[j+1] - sum(a[1:(j-1)])
#}
#u <- a / sum(a)
#a <- t(matrix(a, J, M))     # re: above note --- right just like this
#u <- t(matrix(u, J, M))


switch(survey,
    line = A <- rowSums(a) * 2,
    point = A <- rowSums(a))
switch(unitsIn,
    m = A <- A / 1e6,
    km = A <- A)
switch(unitsOut,
    ha = A <- A * 100,
    kmsq = A <- A)



lamParms <- colnames(Xlam)
gamParms <- colnames(Xgam)
omParms <- colnames(Xom)
detParms <- colnames(Xsig)
nAP <- ncol(Xlam)
nGP <- ncol(Xgam)
nOP <- ncol(Xom)
nDP <- ncol(Xsig)

if(immigration){
    iotaParms <- colnames(Xiota)
    nIP <- ncol(Xiota)
}else{
    nIP <- 0
    iotaParms<- character(0)
    }


if(identical(fix, "gamma")) {
    if(!identical(dynamics, "constant"))
        stop("dynamics must be constant when fixing gamma or omega")
    if(nGP > 1)
        stop("gamma covariates not allowed when fix==gamma")
    else {
        nGP <- 0
        gamParms <- character(0)
    }
}
else if(identical(dynamics, "notrend")) {
    if(nGP > 1)
        stop("gamma covariates not allowed when dyamics==notrend")
    else {
        nGP <- 0
        gamParms <- character(0)
    }
}

if(identical(fix, "omega")) {
    if(!identical(dynamics, "constant"))
        stop("dynamics must be constant when fixing gamma or omega")
    if(nOP > 1)
        stop("omega covariates not allowed when fix==omega")
    else {
        nOP <- 0
        omParms <- character(0)
    }
} else if(identical(dynamics, "trend")) {
    if(nOP > 1)
        stop("omega covariates not allowed when dynamics='trend'")
    else {
        nOP <- 0
        omParms <- character(0)
    }
}


nP <- nAP + nGP + nOP + nDP + (mixture!="P") + (keyfun == "hazard")
if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))

ym <- matrix(y, nrow=M)


# Create indices (should be written in C++), all possible combinatinos of survivors and recruits, finding all unique likelihood transitions
I <- cbind(rep(k, times=lk),
           rep(k, each=lk))
I1 <- I[I[,1] <= I[,2],]
Ib <- Ip <- list()
for(i in 1:nrow(I)) {
    Z <- 0:min(I[i,])
    Ib[[i]] <- which((I1[,1] %in% Z) & (I1[,2]==I[i,1])) - 1
    Ip[[i]] <- as.integer(I[i,2]-Z)
}

#Some k-related indices to avoid repeated calculations in likelihood
lfac.k <- lgamma(k+1)
kmyt <- array(0, c(M, T, lk))
lfac.kmyt <- array(0, c(M, T, lk))
fin <- array(NA, c(M, T, lk))
naflag <- array(NA, c(M, T, J))
for(i in 1:M) {
  for(t in 1:T) {
    fin[i,t,] <- k - yt[i,t] >= 0
    naflag[i,t,] <- is.na(y[i,,t])
    if(!all(naflag[i,t,])) {
      kmyt[i,t,] <- k - yt[i,t]
      lfac.kmyt[i, t, ] <- lgamma(kmyt[i, t, ] + 1)
      #lfac.kmyt[i, t, fin[i,t,]] <- lgamma(kmyt[i, t, fin[i,t,]] + 1)
    }
  }
}
fin <- fin*1 #convert to numeric

nll <- function(parms) {
    beta.lam <- parms[1:nAP]
    beta.gam <- parms[(nAP+1):(nAP+nGP)]
    beta.om <- parms[(nAP+nGP+1):(nAP+nGP+nOP)]
    beta.sig <- parms[(nAP+nGP+nOP+1):(nAP+nGP+nOP+nDP)]
    beta.iota<- parms[(nAP + nGP + nOP + 1):(nAP + nGP + nOP + nDP + nIP)]
# 12/31/2015 Andy added arguments "survey" and also "nintervals" -- number of integration slices
    log.alpha <- 1
    scale<- -99.0
    if(mixture %in% c("NB", "ZIP") & keyfun=="hazard"){
        log.alpha <- parms[nP]
        scale<- -1*exp(parms[nP+1])
    }
    if(mixture %in% c("P") & keyfun=="hazard"){
        log.alpha <- 1
        scale<- -1*exp(parms[nP])
    }
    yperm <- aperm(y, c(1,3,2)) #easier orientation to use in c++
    lgy1 <- lgamma(yperm + 1)


    .Call("nll_distsampOpen",
          yperm, yt,
          Xlam, Xgam, Xom, Xsig, Xiota,
          beta.lam, beta.gam, beta.om, beta.sig, beta.iota, log.alpha,
          Xlam.offset, Xgam.offset, Xom.offset, Xsig.offset, Xiota.offset,
          ytna, yna, #yna not needed
          lk, mixture, first, last, M, J, T,
          delta, dynamics, survey, fix, go.dims, immigration,
          I, I1, Ib, Ip,
          scale,
          a, 
          t(u), #easier to use transpose
          w, db, keyfun, lfac.k, lfac.kmyt, kmyt, lgy1, fin,

          PACKAGE = "unmarked")
}

if(missing(starts))
    starts <- rep(0, nP)
fm <- optim(starts, nll, method=method, hessian=se, ...)
opt <- fm
ests <- fm$par
if(identical(mixture,"NB"))
    nbParm <- "alpha"
else if(identical(mixture, "ZIP"))
    nbParm <- "psi"
else
    nbParm <- character(0)
names(ests) <- c(lamParms, gamParms, omParms, detParms, iotaParms, nbParm)
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
detEstimates <- unmarkedEstimate(name = "Detection", short.name = "sigma",
    estimates = ests[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)],
    covMat = as.matrix(covMat[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP),
        (nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]),
        invlink = "exo", invlinkGrad = "exp")
estimateList <- unmarked:::unmarkedEstimateList(list(lambda=lamEstimates))
gamName <- switch(dynamics,
                  constant = "gamConst",
                  autoreg = "gamAR",
                  notrend = "",
                  trend = "gamTrend")
if(!(identical(fix, "gamma") | identical(dynamics, "notrend")))
    estimateList@estimates$gamma <- unmarkedEstimate(name = "Recruitment",
        short.name = gamName, estimates = ests[(nAP+1) : (nAP+nGP)],
        covMat = as.matrix(covMat[(nAP+1) :
                           (nAP+nGP), (nAP+1) : (nAP+nGP)]),
        invlink = "exp", invlinkGrad = "exp")
if(!(identical(fix, "omega") | identical(dynamics, "trend")))
    estimateList@estimates$omega <- unmarkedEstimate(
        name="Apparent Survival",
        short.name = "omega", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
        covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
            (nAP+nGP+1) : (nAP+nGP+nOP)]),
        invlink = "logistic", invlinkGrad = "logistic.grad")
estimateList@estimates$det <- detEstimates
if(immigration) {
    estimateList@estimates$iota <- unmarkedEstimate(
        name="Immigration",
        short.name = "iota", estimates = ests[(nAP+nGP+nOP+nDP+1) :(nAP+nGP+nOP+nDP+nIP)],
        covMat = as.matrix(covMat[(nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP),
            (nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP)]),
        invlink = "exp", invlinkGrad = "exp")
}
if(identical(mixture, "NB")) {
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
        invlinkGrad = "exp")
    }
if(identical(mixture, "ZIP")) {
    estimateList@estimates$psi <- unmarkedEstimate(name = "Zero-inflation",
        short.name = "psi", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "logistic",
        invlinkGrad = "logistic.grad")
    }
umfit <- new("unmarkedFitDSO", fitType = "distsampOpen",
    call = match.call(), formula = formula, formlist = formlist, data = data,
    sitesRemoved=D$removed.sites, estimates = estimateList, AIC = fmAIC,
    opt = opt, negLogLike = fm$value, nllFun = nll, K = K, mixture = mixture,
    dynamics = dynamics, immigration=immigration, keyfun=keyfun, unitsOut=unitsOut)
return(umfit)
}










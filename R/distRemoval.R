setClass("unmarkedFrameDistRemoval",
  representation(
    yDistance = "matrix",
    yRemoval = "matrix",
    survey = "character",
    dist.breaks = "numeric",
    unitsIn = "character"
  ),
  contains="unmarkedMultFrame"
)

unmarkedFrameDistRemoval <- function(yDistance, yRemoval, numPrimary=1,
                                     siteCovs=NULL, obsCovs=NULL,
                                     yearlySiteCovs=NULL, dist.breaks,
                                     unitsIn){

  # input checking here eventually
  umf <- new("unmarkedFrameDistRemoval", y=yRemoval, yDistance=yDistance,
             yRemoval=yRemoval, numPrimary=numPrimary, siteCovs=siteCovs,
             obsCovs=obsCovs, yearlySiteCovs=yearlySiteCovs, survey="point",
             dist.breaks=dist.breaks, unitsIn=unitsIn, obsToY=diag(ncol(yRemoval)))
  umf <- umf_to_factor(umf)
  umf
}

setAs("unmarkedFrameDistRemoval", "data.frame", function(from){


  out <- callNextMethod(from, "data.frame")
  J <- obsNum(from)
  out <- out[,5:ncol(out), drop=FALSE]

  yDistance <- from@yDistance
  colnames(yDistance) <- paste0("yDist.",1:ncol(yDistance))

  yRemoval <- from@yRemoval
  colnames(yRemoval) <- paste0("yRem.",1:ncol(yRemoval))

  data.frame(yDistance, yRemoval, out)
})

 # bracketing doesn't work yet


setMethod("getDesign", "unmarkedFrameDistRemoval",
  function(umf, formula, na.rm=TRUE){

  M <- numSites(umf)
  T <- umf@numPrimary
  Rdist <- ncol(umf@yDistance)
  Jdist <- Rdist/T
  Rrem <- ncol(umf@yRemoval)
  Jrem <- Rrem/T
  yRem <- as.vector(t(umf@yRemoval))
  yDist <- as.vector(t(umf@yDistance))

  sc <- siteCovs(umf)
  oc <- obsCovs(umf)
  ysc <- yearlySiteCovs(umf)

  if(is.null(sc)) sc <- data.frame(.dummy=rep(0, M))
  if(is.null(ysc)) ysc <- data.frame(.dummy=rep(0, M*T))
  if(is.null(oc)) oc <- data.frame(.dummy=rep(0, M*Rrem))

  ysc <- cbind(ysc, sc[rep(1:M, each=T),,drop=FALSE])
  oc <- cbind(oc, ysc[rep(1:nrow(ysc), each=Jrem),,drop=FALSE])

  Xlam <- model.matrix(formula$lambdaformula,
            model.frame(formula$lambdaformula, sc, na.action=NULL))

  Xphi <- model.matrix(formula$phiformula,
            model.frame(formula$phiformula, ysc, na.action=NULL))

  Xdist <- model.matrix(formula$distanceformula,
            model.frame(formula$distanceformula, ysc, na.action=NULL))

  Xrem <- model.matrix(formula$removalformula,
            model.frame(formula$removalformula, oc, na.action=NULL))

  list(yDist=yDist, yRem=yRem, Xlam=Xlam, Xphi=Xphi, Xdist=Xdist, Xrem=Xrem)
})

setClass("unmarkedFitDistRemoval", contains = "unmarkedFitGDS")

distRemoval <- function(lambdaformula=~1, phiformula=~1, removalformula=~1,
  distanceformula=~1, data, keyfun=c("halfnorm", "exp", "hazard", "uniform"),
  output=c("abund", "density"), unitsOut=c("ha", "kmsq"),
  mixture=c('P', 'NB'), K, starts, method = "BFGS", se = TRUE, engine=c("C","R"),
  rel.tol=1e-4, threads=1, ...){

  keyfun <- match.arg(keyfun)
  output <- match.arg(output)
  unitsOut <- match.arg(unitsOut)
  mixture <- match.arg(mixture)

  formlist <- mget(c("lambdaformula", "phiformula", "distanceformula", "removalformula"))

  M <- numSites(data)
  T <- data@numPrimary
  Rdist <- ncol(data@yDistance)
  Rrem <- ncol(data@yRemoval)
  mixture_code <- switch(mixture, P={1}, NB={2})

  gd <- getDesign(data, formlist)

  # Parameters-----------------------------------------------------------------
  n_param <- c(ncol(gd$Xlam), ifelse(mixture=="P",0,1),
              ifelse(T>1,ncol(gd$Xphi),0),
              ifelse(keyfun=="uniform", 0, ncol(gd$Xdist)),
              ifelse(keyfun=="hazard",1,0),
              ncol(gd$Xrem))
  nP <- sum(n_param)

  pnames <- colnames(gd$Xlam)
  if(mixture!="P") pnames <- c(pnames, "alpha")
  if(data@numPrimary > 1) pnames <- c(pnames, colnames(gd$Xphi))
  if(keyfun!="uniform") pnames <- c(pnames, colnames(gd$Xdist))
  if(keyfun=="hazard") pnames <- c(pnames, "scale")
  pnames <- c(pnames, colnames(gd$Xrem))

  # Distance info--------------------------------------------------------------
  db <- data@dist.breaks
  w <- diff(db)
  ua <- getUA(data) #in utils.R
  u <- ua$u; a <- ua$a
  A <- rowSums(a)
  switch(data@unitsIn, m = A <- A / 1e6, km = A <- A)
  switch(unitsOut,ha = A <- A * 100, kmsq = A <- A)
  if(output=='abund') A <- rep(1, numSites(data))

  # Get K----------------------------------------------------------------------
  if(missing(K) || is.null(K)) K <- max(c(gd$yDist), na.rm=TRUE) + 20
  k <- 0:K
  lk <- length(k)
  lfac.k <- lgamma(k+1)
  #Kmin <- rep(0, M) # Needs to be all 0 for this to work

  # Need separate k info for each part of detection model
  # Distance
  ytDist <- array(gd$yDist, c(Rdist/T, T, numSites(data)))
  ytDist <- aperm(ytDist, c(3,1,2))
  ytDist <- apply(ytDist, c(1,3), function(x) {
      if(all(is.na(x))) return(NA)
      else return(sum(x, na.rm=TRUE))
  })

  kmyt_dist <- array(NA, c(M, T, lk))
  lfac_kmyt_dist <- array(0, c(M, T, lk))
  for(i in 1:M) {
    for(t in 1:T) {
      kmyt_dist[i,t,] <- k - ytDist[i,t]
      lfac_kmyt_dist[i, t, ] <- lgamma(kmyt_dist[i, t, ] + 1)
    }
  }

  # Removal
  ytRem <- array(gd$yRem, c(Rrem/T, T, numSites(data)))
  ytRem <- aperm(ytRem, c(3,1,2))
  ytRem <- apply(ytRem, c(1,3), function(x) {
      if(all(is.na(x))) return(NA)
      else return(sum(x, na.rm=TRUE))
  })

  kmyt_rem <- array(NA, c(M, T, lk))
  lfac_kmyt_rem <- array(0, c(M, T, lk))
  for(i in 1:M) {
    for(t in 1:T) {
      kmyt_rem[i,t,] <- k - ytRem[i,t]
      lfac_kmyt_rem[i, t, ] <- lgamma(kmyt_rem[i, t, ] + 1)
    }
  }

  kmyt_rem <- aperm(kmyt_rem, c(3,2,1))
  Kmin <- apply(cbind(apply(ytDist, 1, max, na.rm=T), apply(ytRem, 1, max, na.rm=T)),1,max,na.rm=T)

  if(missing(starts)){
    starts <- rep(0, nP)
    starts[sum(n_param[1:3])+1] <- log(mean(db))
  } else if(length(starts)!=nP){
    stop(paste0("starts must be length ",sum(n_param)), call.=FALSE)
  }

  nll <- function(param){
    nll_distRemoval(param, n_param, gd$yDist, gd$yRem, mixture_code, keyfun,
                    gd$Xlam, A, gd$Xphi, gd$Xrem, gd$Xdist, db, a, t(u), w,
                    k, lfac.k, lfac_kmyt_dist, kmyt_dist,
                    lfac_kmyt_rem, kmyt_rem, Kmin, threads=2)
  }

  opt <- optim(starts, nll, method=method, hessian=se, ...)

  covMat <- invertHessian(opt, nP, se)
  ests <- opt$par
  fmAIC <- 2 * opt$value + 2 * nP

  names(ests) <- pnames

  lam_ind <- 1:n_param[1]
  lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
    estimates = ests[lam_ind],
    covMat = as.matrix(covMat[lam_ind, lam_ind]), invlink = "exp",
    invlinkGrad = "exp")

  estimateList <- unmarkedEstimateList(list(lambda=lamEstimates))

  if(mixture!="P"){
    a_ind <- n_param[1]+1
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = ests[a_ind],
        covMat = as.matrix(covMat[a_ind, a_ind]), invlink = "exp",
        invlinkGrad = "exp")
  }

  if(T>1){
    phi_ind <- (sum(n_param[1:2])+1):(sum(n_param[1:3]))
    estimateList@estimates$phi <- unmarkedEstimate(name = "Availability",
        short.name = "phi", estimates = ests[phi_ind],
        covMat = as.matrix(covMat[phi_ind, phi_ind]), invlink = "logistic",
        invlinkGrad = "logistic.grad")
  }

  rem_ind <- (sum(n_param[1:5])+1):(sum(n_param[1:6]))
  estimateList@estimates$rem <- unmarkedEstimate(name = "Removal",
      short.name = "rem", estimates = ests[rem_ind],
      covMat = as.matrix(covMat[rem_ind, rem_ind]), invlink = "logistic",
      invlinkGrad = "logistic.grad")

  dist_ind <- (sum(n_param[1:3])+1):(sum(n_param[1:4]))
  estimateList@estimates$dist <- unmarkedEstimate(name = "Distance",
      short.name = "dist", estimates = ests[dist_ind],
      covMat = as.matrix(covMat[dist_ind, dist_ind]), invlink = "exp",
      invlinkGrad = "exp")

  if(keyfun=="hazard"){
    sc_ind <- (sum(n_param[1:4])+1)
    estimateList@estimates$scale <- unmarkedEstimate(name = "Hazard-rate (scale)",
        short.name = "scale", estimates = ests[sc_ind],
        covMat = as.matrix(covMat[sc_ind, sc_ind]), invlink = "exp",
        invlinkGrad = "exp")
  }

  new("unmarkedFitDistRemoval", fitType = "distRemoval",
    call = match.call(), formula = as.formula(paste(formlist, collapse="")),
    formlist = formlist, data = data, estimates = estimateList, sitesRemoved = numeric(0),
    AIC = fmAIC, opt = opt, negLogLike = opt$value, nllFun = nll,
    mixture=mixture, K=K, keyfun=keyfun, unitsOut=unitsOut, output=output)

}

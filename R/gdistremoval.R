setClass("unmarkedFrameGDR",
  representation(
    yDistance = "matrix",
    yRemoval = "matrix",
    survey = "character",
    dist.breaks = "numeric",
    unitsIn = "character",
    period.lengths = "numeric"
  ),
  contains="unmarkedMultFrame"
)

unmarkedFrameGDR <- function(yDistance, yRemoval, numPrimary=1,
                                     siteCovs=NULL, obsCovs=NULL,
                                     yearlySiteCovs=NULL, dist.breaks,
                                     unitsIn, period.lengths=NULL){

  if(is.null(period.lengths)){
    period.lengths <- rep(1, ncol(yRemoval)/numPrimary)
  }

  # input checking here eventually
  umf <- new("unmarkedFrameGDR", y=yRemoval, yDistance=yDistance,
             yRemoval=yRemoval, numPrimary=numPrimary, siteCovs=siteCovs,
             obsCovs=obsCovs, yearlySiteCovs=yearlySiteCovs, survey="point",
             dist.breaks=dist.breaks, unitsIn=unitsIn, period.lengths=period.lengths,
             obsToY=diag(ncol(yRemoval)))
  umf <- umf_to_factor(umf)
  umf
}

setAs("unmarkedFrameGDR", "data.frame", function(from){


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


setMethod("getDesign", "unmarkedFrameGDR",
  function(umf, formula, na.rm=TRUE, return.frames=FALSE){

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

  if(return.frames) return(list(sc=sc, ysc=ysc, oc=oc))

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

setClass("unmarkedFitGDR", contains = "unmarkedFitGDS")

gdistremoval <- function(lambdaformula=~1, phiformula=~1, removalformula=~1,
  distanceformula=~1, data, keyfun=c("halfnorm", "exp", "hazard", "uniform"),
  output=c("abund", "density"), unitsOut=c("ha", "kmsq"), mixture=c('P', 'NB', 'ZIP'),
  K, starts, method = "BFGS", se = TRUE, engine=c("C","R"), threads=1, ...){

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

  Jdist <- Rdist / T
  ysum <- array(t(gd$yDist), c(Jdist, T, M))
  ysum <- t(apply(ysum, c(2,3), sum, na.rm=T))

  Kmin = apply(ysum, 1, max, na.rm=T)

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
  umf_new <- data
  umf_new@y <- umf_new@yDistance
  ua <- getUA(umf_new) #in utils.R
  u <- ua$u; a <- ua$a
  A <- rowSums(a)
  switch(data@unitsIn, m = A <- A / 1e6, km = A <- A)
  switch(unitsOut,ha = A <- A * 100, kmsq = A <- A)
  if(output=='abund') A <- rep(1, numSites(data))

  # Removal info---------------------------------------------------------------
  pl <- data@period.lengths

  # Get K----------------------------------------------------------------------
  if(missing(K) || is.null(K)) K <- max(Kmin, na.rm=TRUE) + 40

  if(missing(starts)){
    starts <- rep(0, nP)
    starts[sum(n_param[1:3])+1] <- log(median(db))
  } else if(length(starts)!=nP){
    stop(paste0("starts must be length ",sum(n_param)), call.=FALSE)
  }

  nll <- function(param){
    nll_gdistremoval(param, n_param, gd$yDist, gd$yRem, ysum, mixture_code, keyfun,
                     gd$Xlam, A, gd$Xphi, gd$Xrem, gd$Xdist, db, a, t(u), w, pl,
                     K, Kmin, threads=threads)
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

  new("unmarkedFitGDR", fitType = "gdistremoval",
    call = match.call(), formula = as.formula(paste(formlist, collapse="")),
    formlist = formlist, data = data, estimates = estimateList, sitesRemoved = numeric(0),
    AIC = fmAIC, opt = opt, negLogLike = opt$value, nllFun = nll,
    mixture=mixture, K=K, keyfun=keyfun, unitsOut=unitsOut, output=output)

}

# Methods

setMethod("predict", "unmarkedFitGDR", function(object, type, newdata,
                                                level=0.95, ...){

  type <- match.arg(type, c("lambda", "phi", "rem", "dist"))
  nm <- switch(type, lambda="lam", phi="phi", rem="rem", dist="dist")
  est <- object[ifelse(nm=="lam","lambda",nm)]

  if(missing(newdata)){
    gd <- getDesign(object@data, object@formlist)
    X <- gd[[paste0("X",nm)]]
  } else{
    if(!inherits(newdata, "data.frame")){
      stop("newdata must be a data frame")
    }
    gd <- getDesign(object@data, object@formlist, return.frames=TRUE)
    fname <- switch(type, lambda="lambda", phi="phi", rem="removal", dist="distance")
    covs <- switch(type, lambda="sc", phi="ysc", rem="oc", dist="ysc")
    X <- make_mod_matrix(object@formlist[[paste0(fname,"formula")]],
                         gd[[covs]], newdata=newdata, re.form=NULL)$X
  }

  stats <- t(sapply(1:nrow(X), function(i){
              bt <- backTransform(linearComb(est, X[i,]))
              ci <- confint(bt, level=level)
              c(Predicted=coef(bt), SE=SE(bt), lower=ci[1], upper=ci[2])
            }))
  as.data.frame(stats)
})

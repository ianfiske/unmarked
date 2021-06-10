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
  mixture_code <- switch(mixture, P={1}, NB={2}, ZIP={3})

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

  rem_ind <- (sum(n_param[1:5])+1):(sum(n_param[1:6]))
  estimateList@estimates$rem <- unmarkedEstimate(name = "Removal",
      short.name = "rem", estimates = ests[rem_ind],
      covMat = as.matrix(covMat[rem_ind, rem_ind]), invlink = "logistic",
      invlinkGrad = "logistic.grad")

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

  if(is.null(level)){
    pred <- do.call(est@invlink, list(drop(X %*% est@estimates)))
    names(pred) <- NULL
    return(data.frame(Predicted=pred, SE=NA, lower=NA, upper=NA))
  }

  stats <- t(sapply(1:nrow(X), function(i){
              bt <- backTransform(linearComb(est, X[i,]))
              ci <- confint(bt, level=level)
              c(Predicted=coef(bt), SE=SE(bt), lower=ci[1], upper=ci[2])
            }))
  as.data.frame(stats)
})

setMethod("getP", "unmarkedFitGDR", function(object){

  M <- numSites(object@data)
  T <- object@data@numPrimary
  Jrem <- ncol(object@data@yRemoval)/T
  Jdist <- ncol(object@data@yDistance)/T

  rem <- predict(object, "rem", level=NULL)$Predicted
  rem <- array(rem, c(Jrem, T, M))
  rem <- aperm(rem, c(3,1,2))

  pif <- array(NA, dim(rem))
  int_times <- object@data@period.lengths
  removalPiFun2 <- makeRemPiFun(int_times)
  for (t in 1:T){
    pif[,,t] <- removalPiFun2(rem[,,t])
  }

  phi <- rep(1, M*T)
  if(T>1) phi <- predict(object, "phi", level=NULL)$Predicted
  phi <- matrix(phi, M, T, byrow=TRUE)

  keyfun <- object@keyfun
  sig <- predict(object, "dist", level=NULL)$Predicted
  sig <- matrix(sig, M, T, byrow=TRUE)
  if(keyfun=="hazard") scale <- exp(coef(object, type="scale"))

  db <- object@data@dist.breaks
  a <- u <- rep(NA, Jdist)
  a[1] <- pi*db[2]^2
  for (j in 2:Jdist){
    a[j] <- pi*db[j+1]^2 - sum(a[1:(j-1)])
  }
  u <- a/sum(a)

  cp <- array(NA, c(M, Jdist, T))
  kf <- switch(keyfun, halfnorm=grhn, exp=grexp, hazard=grhaz,
               uniform=NULL)

  for (m in 1:M){
    for (t in 1:T){
      if(object@keyfun == "uniform"){
        cp[m,,t] <- u
      } else {
        for (j in 1:Jdist){
          cl <- call("integrate", f=kf, lower=db[j], upper=db[j+1], sigma=sig[m])
          names(cl)[5] <- switch(keyfun, halfnorm="sigma", exp="rate",
                                 hazard="shape")
          if(keyfun=="hazard") cl$scale=scale
          cp[m,j,t] <- eval(cl)$value * 2*pi / a[j] * u[j]
        }
      }
    }
  }

  #p_rem <- apply(pif, c(1,3), sum)
  #p_dist <- apply(cp, c(1,3), sum)

  out <- list(dist=cp, rem=pif)
  if(T > 1) out$phi <- phi
  out
})

setMethod("fitted", "unmarkedFitGDR", function(object){

  T <- object@data@numPrimary

  lam <- predict(object, "lambda", level=NULL)$Predicted
  gp <- getP(object)
  rem <- gp$rem
  dist <- gp$dist
  if(T > 1) phi <- gp$phi
  p_rem <- apply(rem, c(1,3), sum)
  p_dist <- apply(dist, c(1,3), sum)

  for (t in 1:T){
    rem[,,t] <- rem[,,t] * p_dist[,rep(t, ncol(rem[,,t]))]
    dist[,,t] <- dist[,,t] * p_rem[,rep(t,ncol(dist[,,t]))]
    if(T > 1){
      rem[,,t] <- rem[,,t] * phi[,rep(t, ncol(rem[,,t]))]
      det[,,t] <- dist[,,t] * phi[,rep(t, ncol(dist[,,t]))]
    }
  }

  if(T > 1){
    rem_final <- rem[,,t]
    dist_final <- det[,,t]
    for (t in 1:T){
      rem_final <- cbind(rem_final, rem[,,t])
      dist_final <- cbind(dist_final, dist[,,t])
    }
  } else {
    rem_final <- drop(rem)
    dist_final <- drop(dist)
  }

  ft_rem <- lam * rem_final
  ft_dist <- lam * dist_final
  list(dist=ft_dist, rem=ft_rem)
})

setMethod("residuals", "unmarkedFitGDR", function(object){
  ft <- fitted(object)
  list(dist=object@data@yDistance - ft$dist, rem=object@data@yRemoval-ft$rem)
})

# ranef

setMethod("ranef", "unmarkedFitGDR", function(object){

  M <- numSites(object@data)
  T <- object@data@numPrimary
  K <- object@K
  mixture <- object@mixture

  Rdist <- ncol(object@data@yDistance)
  Jdist <- Rdist / T
  ysum <- array(t(object@data@yDistance), c(Jdist, T, M))
  ysum <- t(apply(ysum, c(2,3), sum, na.rm=T))
  Kmin = apply(ysum, 1, max, na.rm=T)

  lam <- predict(object, "lambda", level=NULL)$Predicted
  if(object@mixture != "P"){
    alpha <- backTransform(object, "alpha")@estimate
  }

  dets <- getP(object)
  phi <- matrix(1, M, T)
  if(T > 1){
    phi <- dets$phi
  }
  cp <- dets$dist
  pif <- dets$rem

  pr <- apply(cp, c(1,3), sum)
  prRem <- apply(pif, c(1,3), sum)

  post <- array(0, c(M, K+1, T))
  colnames(post) <- 0:K
  for (i in 1:M){
    if(mixture=="P"){
      f <- dpois(0:K, lam[i])
    } else if(mixture=="NB"){
      f <- dnbinom(0:K, mu=lam[i], size=alpha)
    } else if(mixture=="ZIP"){
      f <- dzip(0:K, lam[i], alpha)
    }
    g <- rep(1, K+1)
    for (t in 1:T){
      for (k in 1:(K+1)){
        g[k] <- g[k] * dbinom(ysum[i,t], k-1, prob=pr[i,t]*prRem[i,t]*phi[i,t],
                              log=FALSE)
      }
    }
    fg <- f*g
    post[i,,t] <- fg/sum(fg)
  }

  new("unmarkedRanef", post=post)
})


setMethod("simulate", "unmarkedFitGDR", function(object, nsim, seed=NULL, na.rm=FALSE){

  lam <- predict(object, "lambda", level=NULL)$Predicted
  dets <- getP(object)

  if(object@mixture != "P"){
    alpha <- backTransform(object, "alpha")@estimate
  }

  M <- length(lam)
  T <- object@data@numPrimary

  if(T > 1){
    phi <- dets$phi
  } else {
    phi <- matrix(1, M, T)
  }

  Jrem <- dim(dets$rem)[2]
  Jdist <- dim(dets$dist)[2]

  p_dist <- apply(dets$dist, c(1,3), sum)
  p_rem <- apply(dets$rem, c(1,3), sum)

  dist_scaled <- array(NA, dim(dets$dist))
  rem_scaled <- array(NA, dim(dets$rem))
  for (t in 1:T){
    dist_scaled[,,t] <- dets$dist[,,t] / p_dist[,t]
    rem_scaled[,,t] <- dets$rem[,,t] / p_rem[,t]
  }

  p_total <- p_dist * p_rem * phi
  stopifnot(dim(p_total) == c(M, T))

  out <- vector("list", nsim)

  for (i in 1:nsim){

    switch(object@mixture,
      P = N <- rpois(M, lam),
      NB = N <- rnbinom(M, size=alpha, mu=lam),
      ZIP = N <- rzip(M, lam, alpha)
    )

    ydist <- matrix(NA, M, T*Jdist)
    yrem <- matrix(NA, M, T*Jrem)

    for (m in 1:M){
      ysum <- rbinom(T, N[m], p_total[m,])

      ydist_m <- yrem_m <- c()

      for (t in 1:T){
        rem_class <- sample(1:Jrem, ysum[t], replace=TRUE, prob=rem_scaled[m,,t])
        rem_class <- factor(rem_class, levels=1:Jrem)
        yrem_m <- c(yrem_m, as.numeric(table(rem_class)))
        dist_class <- sample(1:Jdist, ysum[t], replace=TRUE, prob=dist_scaled[m,,t])
        dist_class <- factor(dist_class, levels=1:Jdist)
        ydist_m <- c(ydist_m, as.numeric(table(dist_class)))
      }
      stopifnot(length(ydist_m)==ncol(ydist))
      stopifnot(length(yrem_m)==ncol(yrem))

      ydist[m,] <- ydist_m
      yrem[m,] <- yrem_m
    }
    out[[i]] <- list(yRemoval=yrem, yDistance=ydist)
  }
  out
})


setMethod("update", "unmarkedFitGDR",
    function(object, lambdaformula, phiformula, removalformula, distanceformula,
             ..., evaluate = TRUE)
{

    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")

    if(!missing(lambdaformula)){
      call$lambdaformula <- lambdaformula
    }
    if(!missing(phiformula)){
      if(!is.null(call$phiformula)){
        call$phiformula <- phiformula
      }
    }
    if(!missing(removalformula)){
      call$removalformula <- removalformula
    }
    if(!missing(distanceformula)){
      call$distanceformula <- distanceformula
    }

    extras <- match.call(call=sys.call(-1),
                         expand.dots = FALSE)$...
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing])
            call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
            }
        }
    if (evaluate)
        eval(call, parent.frame(2))
    else call
})


setMethod("replaceY", "unmarkedFrameGDR",
          function(object, newY, replNA=TRUE, ...){

      ydist <- newY$yDistance
      stopifnot(dim(ydist)==dim(object@yDistance))
      yrem <- newY$yRemoval
      stopifnot(dim(yrem)==dim(object@yRemoval))

      if(replNA){
        ydist[is.na(object@yDistance)] <- NA
        yrem[is.na(object@yRemoval)] <- NA
      }

      object@yDistance <- ydist
      object@yRemoval <- yrem
      object
})


setMethod("SSE", "unmarkedFitGDR", function(fit, ...){
    r <- sapply(residuals(fit), function(x) sum(x^2, na.rm=T))
    return(c(SSE = sum(r)))
})

setMethod("nonparboot", "unmarkedFitGDR",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
   stop("Not currently supported for unmarkedFitGDR", call.=FALSE)
})


setMethod("plot", c(x = "unmarkedFitGDR", y = "missing"), function(x, y, ...)
{
    r <- residuals(x)
    e <- fitted(x)

    old_mfrow <- graphics::par("mfrow")
    on.exit(graphics::par(mfrow=old_mfrow))
    graphics::par(mfrow=c(2,1))

    plot(e[[1]], r[[1]], ylab="Residuals", xlab="Predicted values",
         main="Distance")
    abline(h = 0, lty = 3, col = "gray")

    plot(e[[2]], r[[2]], ylab="Residuals", xlab="Predicted values",
         main="Removal")
    abline(h = 0, lty = 3, col = "gray")
})

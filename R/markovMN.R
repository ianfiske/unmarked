#' @include utils.R
roxygen()

# TODO:  improve documentation!!
# TODO:  improve mechanism for choosing the detection matrix!

#' Estimate parameters of the general multiseason multistate occpancy model.
#'
#' Site level covariates are currently not implemented for markovMN and so stateformula is ignored.
#' See \link{unmarked} for a discussion of detformula.  See \link{unMarkedFrame} for a description of how to create
#' an unMarkedFrame for supplying data to the argument \code{umf}.
#'
#' \code{K} and \code{phiMatrix} together determine the model form.  Multiple phi matrices
#' may be possible for each K.  Options for phi are:
#'
#'
#' TODO:  describe phi matrices in detail.
#'
#' Currently, selection of the appropriate detection matrix for a given \code{K}
#' is given by \code{arDet}.  If \code{arDet} is TRUE, then the autoregressive flavored
#' detection matrix (reduced form) is chose.
#'
#' TODO:  describe detection matrices in detail.
#'
#'
#' Each freely varying detection parameter can be modeled as a linear function of observation level covariates.
#' \code{detconstraint} is a matrix with number of rows equal to the detection matrix parameters and number of
#' columns equal to the number of covariates (plus 1 for the intercept).  Each column is a constraint vector beginning with 1 that indicates
#' which covariates are restricted to have the same effect on a given matrix parameter.
#'
#' @title Fit the general multistate multiseason occupancy model.
#' @param stateformula right-hand side formula describing covariates of occurence.
#' @param detformula right-hand side formula describing covariates of detection.
#' @param umf unMarkedFrame object that supplies the data (see \link{unMarkedFrame})..
#' @param detconstraint matrix to describe constraints on detection parameters
#' @param phiconstraint vector to describe phi constraints
#' @param psiconstraint vector to describe psi constraints
#' @param K integer maximum level of y
#' @param phiMatrix character describing which phi to use
#' @param EM logical specifying to use EM (TRUE) or BFGS (FALSE)
#' @param psi.init initial values for psi
#' @param phiParms.init initial values for phi parameters
#' @param detParms.init initial values for detection parameters
#' @param get.inits logical, specifying whether or not to search for better initial values
#' @param booststrap.se logical. TRUE calculates standard errors using bootstrap
#' @param B number of bootstrap interations
#' @param trace logical, TRUE show BFGS steps
#' @param arDet TRUE chooses the autoregressive type detection structure
#' @useDynLib unmarked
#' @export
markovMN <-
  function(stateformula = ~ 1, detformula = ~ 1, umf,
           detconstraint = NULL,
           phiconstraint = NULL, psiconstraint = NULL, K,
           phiMatrix = "cumlogit", EM = FALSE, psi.init = NULL,
           phiParms.init = NULL, detParms.init = NULL,
           get.inits = TRUE, bootstrap.se = FALSE, B = 50, trace = FALSE,
           phiMatFun = 0, detVecFun, n.starts = 1, homo.det = FALSE,
           max.bad = 3, min.good = 1, fitStats = TRUE)
{
  ## truncate at K
  umf@y[umf@y > K] <- K

  ## coerce vector detcon for K = 1.
  if(identical(K,1)) {
    if(class(detconstraint) %in% c("numeric","integer")) {
      detconstraint <-  t(as.matrix(detconstraint))
    }
  }

#  if(K == 3) {
#    if(!(phiMatrix %in%
#         c("4state", "cumlogit", "4state4", "4stateAR","logit4","logit4ar")))
#      stop(paste("Inappropriate phiMatrix specified for K =",K))
#  }
#  if(K == 2) {
#    if(!(phiMatrix %in% c("3state", "cumlogit", "logit3")))
#      stop(paste("Inappropriate phiMatrix specified for K =",K))
#  }

  ##################################
  ## section determines appropriate default detconstraint...
  ## needs more investigation.
  nDCP <- length(attr(terms(detformula), "term.labels")) + 1
#  if(arDet)
#    nDMP <- K
#  else
    nDMP <-  K*(K+1)/2

  if(is.null(detconstraint))
    detconstraint <- matrix(rep(1:nDMP,nDCP),nDMP, nDCP)
    #detconstraint <- matrix(1:(nDMP*nDCP), nDMP, nDCP)
#################################


  ## only do NA checking for variables being used!
  ## create reduced detection formula from detconstraint
  to.rm <- which(apply(detconstraint == 0, 2, all))
  if(length(to.rm) != 0 & length(to.rm) != (ncol(detconstraint) - 1)) {
    detformula.red <- eval(parse(text=paste("~ ",
                               paste(attr(terms(detformula), "term.labels")[-(to.rm-1)],
                                     collapse="+"))))
  } else if(length(to.rm) == (ncol(detconstraint) - 1)) {
    detformula.red <- ~1
  } else {
    detformula.red <- detformula
  }

  umf <- handleNA(stateformula, detformula.red, umf)
  y <- umf@y
  J <- umf@obsNum / umf@primaryNum

  M <- nrow(y)
  nY <- ncol(y)/J
  n.det <- sum(apply(y > 0, 1, any, na.rm = TRUE))

  fc <- match.call()
  fc[[1]] <- as.name("markovMN.fit")
  fc$bootstrap.se <- fc$covdata.site <- fc$covdata.obs <- fc$data <-
    fc$B <- NULL
  fc$umf <- as.name("umf")
  fc$J <- as.name("J")
  fc$detconstraint <- as.name("detconstraint")
#  fc$max.bad <- as.name("max.bad")
#  fc$min.good <- as.name("min.good")
  fc$homo.det <- as.name("homo.det")
  fm <- eval(fc)

#  markovMN.fit.args <- names(formals(markovMN.fit))
#  for(a in markovMN.fit.args) {
#    if(!missing(a)) {
#      fc[[a]] <- as.name(a)
#    }
#  }

  if(bootstrap.se) {

    smooth.b <- array(NA, c(K + 1, nY, B))
    psi.b <- matrix(NA, B, K + 1)
    ss.b <- matrix(NA, B, K + 1)

    obsCovs.siteInd <- matrix(1:(M*J*nY),M,J*nY,byrow=T)
    for(b in 1:B) {

      samp <- sample(1:M, M, replace = TRUE)
#      y.b <- y[samp,]
#      obsdata.b <- lapply(obsdata, function(x) x[samp,])
#
#      fc$obsdata <- quote(obsdata.b)
#      fc$y <- quote(y.b)

      umf.b <- umf
      umf.b@y <- umf@y[samp,]
      obsCovsInd.b <- as.vector(t(obsCovs.siteInd[samp,]))
      umf.b@obsCovs <- umf.b@obsCovs[obsCovsInd.b,]
      fc$umf <- quote(umf.b)

      fm.b <- eval(fc)

      if(!(TRUE %in% is.nan(smooth.b[,,b]))) {
        smooth.b[,,b] <- fm.b$smooth
      } else {
        cat("smooth.b contained NaN.\n")
      }
      psi.b[b,] <- fm.b$psi
      ss.b[b,] <- fm.b$ss

      cat(paste("Bootstrap iteration",b,"completed.\n"))
    }
    ss.b.cov <- cov(ss.b, use = "complete.obs")
    psi.b.cov <- cov(psi.b, use = "complete.obs")

    smooth.b.cov <- apply(smooth.b, c(2),
                          function(x) cov(t(x), use = "complete.obs"))
    smooth.b.cov <- array(smooth.b.cov, c(K + 1, K + 1, nY))

    ## also get the time-series style covariances for K=1 only.
    if(identical(K,1)) {
      smooth.mat <- t(smooth.b[2,,])
      fm$smooth.covmat <- cov(smooth.mat, use="complete.obs")
      fm$smooth.cormat <- cor(smooth.mat, use="complete.obs")
    }

    fm$psi.cov <- psi.b.cov
    fm$ss.cov <- ss.b.cov
    fm$smooth.cov <- smooth.b.cov

  }

  fm$n.det <- n.det

  ## compute the projected trajectory
  fm$projected <- matrix(NA, K + 1, nY)
  fm$projected[,1] <- fm$psi
  for(year in 2:nY) {
    fm$projected[,year] <- t(fm$phi) %*%
        fm$projected[,year-1]
  }

  ## check hessian for NaN's, Inf's, 0's, etc.
  ## Their presence should trigger "convergence <- 1"
  if(any(is.na(fm$hessian)) || any(is.infinite(fm$hessian)) ||
      identical(sum(abs(fm$hessian)), 0)) {
    fm$convergence <- 1
  }

  if(fitStats == TRUE) {
    fitStats <- computeFitStats(detMats = fm$detMats, smooth = fm$smooth.sites, y = y, 
			J.it = fm$J.it)
    fm$fitStats <- fitStats
  }

  return(fm)
}


markovMN.fit <- function(stateformula = ~ 1, detformula = ~ 1, umf,
                         detconstraint = NULL,
                         phiconstraint = NULL, psiconstraint = NULL, J, K,
                         phiMatrix = "cumlogit", EM = FALSE, psi.init = NULL,
                         phiParms.init = NULL, detParms.init = NULL,
                         get.inits = TRUE, trace = FALSE,
                         phiMatFun = 0, detVecFun = NULL,
                         homo.det = FALSE, max.bad = 3, min.good = 3)
{
  y <- umf@y

  M <- nrow(y)
  nY <- ncol(y)/J

  designMats <- getDesign(stateformula = stateformula, detformula = detformula, umf)

  V.itj <- designMats$V
  nDCP <- ncol(V.itj)
  detParms <- colnames(V.itj)

#  ## number of free terms in detection matrix.
#  ## these will be modeled by XDet's
#  if(arDet)
#    nDMP <- K
#  else
    nDMP <-  K*(K+1)/2

  if(nrow(detconstraint) != nDMP)
    stop(paste("detconstraint has wrong number of rows.\n
It should have",nDMP))

  ## for performance, remove variables eliminated by detconstraint
  to.rm <- which(apply(detconstraint == 0, 2, all))
  if(length(to.rm) != 0) {
    detconstraint <- detconstraint[,-to.rm]
    V.itj <- V.itj[,-to.rm]
    detParms <- detParms[-to.rm]
    nDCP <- nDCP - length(to.rm)
  }
  if(is.null(dim(detconstraint)))
    detconstraint <- matrix(detconstraint,nDMP,nDCP)

  ## create design matrix for each matrix parameter
  # TODO: make this do the right thing if there are constraints.
  V.itjk <- V.itj %x%  diag(nDMP)
                                        #get a better line here    if(is.null(detconstraint)) detconstraint <- 1:nDMP.un

  fphiMatrix <- paste("f",phiMatrix,sep="")
  ## add "f" to fool switch because it doesn't like characters that start
  ## with numbers
#  nPhiP.un <- switch(fphiMatrix,
#                     fcumlogit = K + 1,
#                     f4state = 12,
#                     f4state4 = 4,
#                     f4stateAR = 10,
#                     f3state = 6,
#                     f2state = 2,
#                     flogit4 = 12,
#                     flogit3 = 6,
#                     flogit2 = 2,
#                     flogit4ar = 6)
	nPhiP.un <- switch(K,  # TODO: this is a placeholder.... not very flexible.
			2, 6, 12)

#  nPhiP.un.table <- matrix(c(
#      1, 0, 2,
#      2, 0, 6,	## K=3, 0: logit3
#      3, 0, 12,  ## 0 = full mat
#      3, 1, 6,   ## 1 = small ord
#      3, 2, 10), ## 2 = ord with interecept
#  5, 3, byrow = TRUE)
#  colnames(nPhiP.un.table) <- c("K", "phiNum", "nPhiP.un")
#
#  nPhiP.un <- nPhiP.un.table[nPhiP.un.table[,2] == phiMatFun &
#          nPhiP.un.table[,1] == K,3]

  nSP.un <- K                # number of parameters for psi vector of initial

  if(is.null(phiconstraint)) phiconstraint <- 1:nPhiP.un
  if(is.null(psiconstraint)) psiconstraint <- 1:nSP.un

#  nPhiP <- switch(fphiMatrix,
#                  fcumlogit = K + 1,
#                  f4state = max(phiconstraint),
#                  f4stateAR = max(phiconstraint),
#                  f4state4 = max(phiconstraint),
#                  f3state = max(phiconstraint),
#                  f2state = max(phiconstraint),
#                  flogit4 = max(phiconstraint),
#                  flogit3 = max(phiconstraint),
#                  flogit2 = max(phiconstraint),
#                  flogit4ar = max(phiconstraint))
  nPhiP <- max(phiconstraint)
  nSP <- max(psiconstraint)

  ## create linked list of parameters
  theta.df <- data.frame(parameter = character(), start = numeric(),
                         end = numeric(), stringsAsFactors = FALSE)

  theta.df <- addParm(theta.df, "phiParms", nPhiP)

  ## convert from new easy-input style detcon to computational format
  prev.col.max <- 0
  detcon.corr <- matrix(0, nDMP, nDCP)
  for(j in 1:nDCP) {
    detcon.corr[,j] <- ifelse(detconstraint[,j] != 0,
                              detconstraint[,j] + prev.col.max, 0)
    prev.col.max <- max(c(detcon.corr[,j],prev.col.max))
  }
  detcon.vec <- as.vector(detcon.corr)
  nDP <- max(detcon.vec)
  nDP.un <- length(detcon.vec)


  H.det <- matrix(0, nDP.un, nDP)
  for(i in 1:length(detcon.vec)){
    H.det[i,detcon.vec[i]] <- 1
  }

  theta.df <- addParm(theta.df, "detParms", nDP)
  nP <- nDP + nSP + nPhiP  # total number of parameters

  ## construct constrain matrix for phi paramters
  H.phi <- matrix(0, nPhiP.un, nPhiP)
  for(i in 1:nPhiP.un){
    H.phi[i, phiconstraint[i]] <- 1
  }

  H.psi <- matrix(0, nSP.un, nSP)
  for(i in 1:nSP.un){
    H.psi[i, psiconstraint[i]] <- 1
  }

  y.itj <- as.numeric(t(y))

  ## reorder X.tjik to be X.itjk
###   t.tjik <- rep(1:nY, each = nDMP *M*J)
###   i.tjik <- rep(rep(1:M, each = nDMP), nY*J)
###   j.tjik <- rep(rep(1:J, each = nDMP*M), nY)
###   k.tjik <- rep(1:nDMP, M * J * nY)

###   XDet.itjk <- XDet.tjik[order(i.tjik, t.tjik, j.tjik, k.tjik),]

  ## replace NA's with 99 before passing to C++
  ## TODO: need better missing data passing mechanism (maybe NaN of Inf?)
  y.itj[is.na(y.itj)] <- 99
  V.itjk[is.na(V.itjk)] <- 9999

  if(is.null(psi.init)) psi.init <- rnorm(nSP)
  if(is.null(phiParms.init)) phiParms.init <- rnorm(nPhiP)
  if(is.null(detParms.init)) detParms.init <- rnorm(nDP)

  # get ragged array indices
  y.it <- matrix(t(y), nY*M, J, byrow = TRUE)
  J.it <- rowSums(!is.na(y.it))

  fm <- findMLE2(y.itj=y.itj, V.itjk = V.itjk, J.it = J.it, nDMP=nDMP, nDP=nDP, nSP=nSP, nSP.un=nSP.un,
                   nPhiP=nPhiP, nP=nP, nDP.un=nDP.un,
                   nPhiP.un=nPhiP.un, H.det=H.det, H.phi=H.phi, H.psi=H.psi, K=K, M=M, J=J, nY=nY,
                   phiMatrix=phiMatrix, EM= EM, psi.init = psi.init,
                   phiParms.init = phiParms.init,
                   detParms.init = detParms.init, get.inits = get.inits,
                   trace = trace,# arDet = arDet,
                   phiMatFun = phiMatFun, detVecFun = detVecFun,
                   homo.det = homo.det, max.bad = max.bad, min.good = min.good)

  phiEsts <- H.phi %*% fm$phiParms
  detEsts <- H.det %*% fm$detParms
  detEsts <- matrix(detEsts, nDMP, nDCP)
  colnames(detEsts) <- detParms
  colnames(detconstraint) <- detParms
  phi <- fm$phi
  psiEsts <- H.psi %*% fm$psiParms
  psi <- exp(c(0,psiEsts))/sum(exp(c(0,psiEsts)))
  ss <- getSS(phi)
  hessian <- fm$hessian

  mle <- c(fm$psiParms, fm$phiParms, fm$detParms)
  parm.names <- c(rep("psi", nSP), rep("phi", nPhiP), rep("det", nDP))
  parm.names <- paste(parm.names, c(1:nSP, 1:nPhiP, 1:nDP), sep="")
  mle.df <- data.frame(names = parm.names, value = mle)
  rownames(mle.df) <- 1:nP

  ## get epsi SE
###     psi.ind <- grep("psi", rawfit.names)
###     psi.theta <- ests[psi.ind]
###     psi.hessian <- fm$hessian[psi.ind,psi.ind]
###     epsi.se <- sqrt(deltaVar(psi.theta, psi.hessian, meanPsi))

###     ## get ess SE
###     phi.ind <- grep("phi", rawfit.names)
###     phi.theta.con <- ests[phi.ind]
###     phi.hessian <- fm$hessian[phi.ind, phi.ind]
###     ess.se <- sqrt(deltaVar(phi.theta.con, phi.hessian, meanSS,
###                             phiconstraint = phiconstraint,
###                             phifun = phiMatrix))

###     ## smoothing SE

#  smooth.mean <- apply(fm$smooth, c(2,3), meanstate)
#  smooth.overall.mean <- rowMeans(smooth.mean)

  smooth <- apply(fm$smooth, c(1,2), mean)
  smooth.sites <- fm$smooth

#  q <- 1
#  w <- 1
#  detMats <- list()
#  for(i in 1:M) {
#    detMats[[i]] <- list()
#    for(t in 1:nY) {
#      detMats[[i]][[t]] <- list()
#      for(j in seq(length=J.it[q])) {
#        detMats[[i]][[t]][[j]] <- matrix(fm$detMats[w:(w+(K+1)^2-1)], K+1,K+1)
#        w = w + (K+1)^2
#      }
#      q = q + 1
#    }
#  }

#  q <- 1
#  w <- 1
#  detMats <- array(NA, c(K+1, K+1, J,nY,M))
#  for(i in 1:M) {
##    detMats[[i]] <- list()
#    for(t in 1:nY) {
##      detMats[[i]][[t]] <- list()
#      for(j in seq(length=J.it[q])) {
#        detMats[,,j,t,i] <- matrix(fm$detMats[w:(w+(K+1)^2-1)], K+1,K+1)
#        w = w + (K+1)^2
#      }
#      q = q + 1
#    }
#  }



###     smooth.cov <- deltaVar(ests, fm$hessian, meanSmooth, y.itj=y.itj,
###                           XDet.itjk=XDet.itjk, nDMP=nDMP, nDCP=nDCP,
###                           nDP=nDP, nDYP=nDP, nSP=nSP, nPhiP=nPhiP, nP=nP,
###                           nDMP.un=nDMP.un, nPhiP.un = nPhiP.un,
###                           H.det = H.det, H.phi = H.phi, K = K,
###                           yearly.det = yearly.det, M=M, J=J, nY=nY,
###                           phiMatrix=phiMatrix)
###     smooth.se <- sqrt(diag(smooth.cov))

  list(mle = mle.df, psiEsts = psiEsts, phiEsts = phiEsts,
       detEsts = detEsts, phi = phi,
       AIC = 2*fm$nll + 2*nP,
       NegLogLike = fm$nll,
       psi = psi, epsi = meanstate(psi), #epsi.se = epsi.se,
       ss = ss, ess = meanstate(ss), #ess.se = ess.se,
       smooth = smooth, #smooth.se = smooth.se,
       #arDet = arDet,
       detform = detformula,
       psiconstraint = psiconstraint,
       phiconstraint = phiconstraint,
       detconstraint = detconstraint,
       phiMatrix = phiMatrix, K = K, hessian = hessian,
       n = M, convergence = fm$convergence, detMats = fm$detMats,
       smooth.sites = smooth.sites, J.it = J.it)
}


findMLE <-
  function(y.itj, V.itjk, J.it, nDMP, nDP, nSP, nSP.un, nPhiP, nP, nDP.un,
           nPhiP.un, H.det, H.phi, H.psi, K, M, J, nY,
           phiMatrix, smooth.only = FALSE, EM,
           psi.init, phiParms.init, detParms.init, get.inits, trace,
           arDet=0,phiMatFun = phiMatFun, detVecFun = detVecFun, homo.det,
           min.good, max.bad)
{
  ncX <- ifelse(is.null(ncol(V.itjk)), 1, ncol(V.itjk))
  x <- .C("findMLE",
          as.integer(y.itj),
          as.double(t(V.itjk)),
          as.integer(t(J.it)),
          as.integer(ncX),
          as.integer(nDMP),
          as.integer(nDP),
          as.integer(nSP),
          as.integer(nSP.un),
          as.integer(nPhiP),
          as.integer(nP),
          as.integer(nDP.un),
          as.integer(nPhiP.un),
          as.double(H.det),
          as.double(H.phi),
          as.double(H.psi),
          as.integer(K),
          as.integer(M),
          as.integer(J),
          as.integer(nY),
          nll = double(1),
          phi = double((K + 1)^2),
          psiParms = double(nSP),
          detParms = double(nDP),
          phiParms = double(nPhiP),
          smooth = double(M*nY*(K + 1)),
          smooth.only = as.integer(smooth.only),
          as.character(phiMatrix),
          hessian = double(nP^2),
          EM = as.integer(EM),
          as.double(psi.init),
          as.double(phiParms.init),
          as.double(detParms.init),
          as.integer(get.inits),
          as.integer(trace),
          as.integer(homo.det),
          as.integer(arDet),
          convergence = integer(1),
          as.integer(K), # replace with detfun in future.
          as.integer(phiMatFun),
          as.integer(min.good),
          as.integer(max.bad),
          detMats = double(sum(J.it)*(K+1)^2)) # store det mat for each obs

  return(list(nll = x$nll, phi = matrix(x$phi,K + 1, K + 1),
              psiParms = x$psiParms, detParms = x$detParms,
              phiParms = x$phiParms,
              smooth = array(x$smooth, c(K + 1, nY, M)),
              hessian = matrix(x$hessian,nP,nP),
              convergence = x$convergence, detMats = x$detMats))
}

#profile1D.C <- function(y.itj, V.itjk, nDMP, nDP, nSP, nSP.un, nPhiP, nP, nDP.un,
#    nPhiP.un, H.det, H.phi, H.psi, K, M, J, nY,
#    phiMatrix,
#    arDet,phiMatFun = phiMatFun, detVecFun = detVecFun,homo.det,
#    profSeq, mle, whichVarying)
#{
#
#}



## phi4state <- function(phiParms) {
##   x <- .C("phi4state", as.double(phiParms), phi = numeric(16))
##   t(matrix(x$phi, 4, 4))
## }

## phi3state <- function(phiParms) {
##   x <- .C("phi3state", as.double(phiParms), phi = numeric(9))
##   t(matrix(x$phi, 3, 3))
## }

## phi2state <- function(phiParms) {
##   x <- .C("phi2state", as.double(phiParms), phi = numeric(4))
##   t(matrix(x$phi, 2, 2))
## }

## phi4state4 <- function(phiParms) {
##   x <- .C("phi4state4", as.double(phiParms), phi = numeric(16))
##   t(matrix(x$phi, 4, 4))
## }

## phi4stateAR <- function(phiParms) {
##   x <- .C("phi4stateAR", as.double(phiParms), phi = numeric(16))
##   t(matrix(x$phi, 4, 4))
## }

## logit4 <- function(phiParms) {
##   x <- .C("phi4logit", as.double(phiParms), phi = numeric(16))
##   t(matrix(x$phi, 4, 4))
## }


# The R version
findMLE2 <-
    function(y.itj, V.itjk, J.it, nDMP, nDP, nSP, nSP.un, nPhiP, nP, nDP.un,
        nPhiP.un, H.det, H.phi, H.psi, K, M, J, nY,
        phiMatrix, smooth.only = FALSE, EM,
        psi.init, phiParms.init, detParms.init, get.inits, trace,
        arDet=0, phiMatFun, detVecFun, homo.det = FALSE,
        max.bad,min.good) {

  alpha <- array(NA, c(K + 1, nY, M))
  beta <- array(NA, c(K + 1, nY, M))
  gamma <- array(NA, c(K + 1, nY, M))

  ## V.arr is array of matrices for dims: 3=i, 4=t, 5=j
  #V.arr <- apply(V.itjk,c(1,2),function(x) x)
  V.arr <- array(t(V.itjk), c(nDP, nDMP, J, nY, M))
  V.arr <- aperm(V.arr, c(2,1,5,4,3))

  y.arr <- array(y.itj, c(J, nY, M))
  y.arr <- aperm(y.arr, c(3:1))
#  library(Matrix)
#  D <- Diagonal(M)
  ## modifies alpha as side effect
  ## returns negloglike
  storage.mode(J.it) <- storage.mode(y.arr) <- storage.mode(K) <- "integer"
  #dyn.load("../../../unmarked_testing/inline_testing2.so")
  
  
  forward <- function(detParms, phi, psi, storeAlpha = FALSE) {
	  
	  negloglike <- 0
	  ##################    version to get all detvec's at once
	  psiSite <- matrix(psi, K + 1, M)
	  #detVec.arr <- array(1, c(K + 1, nY, M))
#    if(homo.det) {
#      detVecMat <- matrix(NA,K+1,K+1)
#      for(i in 1:(K+1)) {
#        detVecMat[,i] <- do.call(detVecFun, list(p = detParms, y = i-1))
#      }
#    } else {
	  mp <- array(V.itjk %*% detParms, c(nDMP, J, nY, M))
#    }
	  
	  for(t in 1:nY) {
		  storage.mode(t) <- "integer"
#browser()
#		.Call("setDetVecs", y.arr, detVec.arr, mp, 
#				J.it[seq(from = t,to = length(J.it)-nY+t, by=nY)], t)
		  detVecs <- .Call("getDetVecs", y.arr, mp, J.it[seq(from = t,to = length(J.it)-nY+t, by=nY)], t, K)
#      for(i in 1:M) {
#        for(j in 1:J) {
#          if(y.arr[i,t,j] != 99) {
#              detVec.arr[,t,i] <- detVec.arr[,t,i] *
#                  do.call(detVecFun, list(p=mp[,j,t,i], y=y.arr[i,t,j]))
#          }
#        }
#      }
		  #psiSite <- psiSite * detVec.arr[,t,]
		  psiSite <- psiSite * detVecs
		  if(storeAlpha) alpha[,t,] <<- psiSite[,]
		  if(t < nY) {
			  psiSite <- phi %*% psiSite
		  } else {
			  negloglike <- negloglike - sum(log(colSums(psiSite)))
		  }
	  }
	  negloglike
  


####################
#
#    mp <- array(V.itjk %*% detParms, c(nDMP, J, nY, M))
#
#    for(i in 1:M) {
#      psiSite <- psi
#      for(t in 1:nY) {
#
#        detVec <- rep(1, K + 1)
#        for(j in 1:J) {
#          if(y.arr[i,t,j] != 99) {
#            detVecObs <- do.call(detVecFun, list(p=mp[,j,t,i], y=y.arr[i,t,j]))
#            detVec <- detVec * detVecObs
#          }
#        }
#        psiSite <- psiSite * detVec
#        if(storeAlpha) alpha[,t,i] <<- psiSite
#        if(t < nY) {
#          psiSite <- phi %*% psiSite
#        } else {
#          negloglike <- negloglike - log(sum(psiSite))
#        }
#
#      }
#    }
#    negloglike
  }

  backward <- function(detParams, phi, psi) {
    for(i in 1:M) {
      backP <- rep(1, K + 1)
      for(t in nY:1) {

        beta[, t, i] <<- backP

        detVec <- rep(1, K + 1)
        for(j in 1:J) {
          if(y.arr[i,t,j] != 99) {
            mp <- V.arr[,,i,t,j] %*% detParams
            #detVecObs <- do.call(detVecFun, list(p=mp, y=y.arr[i,t,j]))
			detVecObs <- .Call("getSingleDetVec", y.arr[i,t,j], mp, K)
			detVec <- detVec * detVecObs
          }
        }

        backP <- t(phi) %*% (detVec * backP)

      }
    }
  }

  getDetMats <- function(detParams, phi, psi) {
	  detMats <- array(NA, c(K+1,K+1,J,nY,M))
	  for(i in 1:M) {
		  for(t in nY:1) {
			  for(j in 1:J) {
				  if(y.arr[i,t,j] != 99) {
					  for(k in 0:K) {
						  mp <- V.arr[,,i,t,j] %*% detParams
						  #detVecObs <- do.call(detVecFun, list(p=mp, y=y.arr[i,t,j]))
						  detMats[,k+1,j,t,i] <- .Call("getSingleDetVec", k, mp, K)
					  }
				  }
			  }
		  }
	  }
	  return(detMats)
  }
  
  
  nll <- function(params) {

    psiParams <- c(0, H.psi %*% params[1:nSP])
    psi <- exp(psiParams)
    psi <- psi/sum(psi)

    phi = do.call(phiMatFun, list(p=H.phi %*% params[(nSP + 1):(nSP + nPhiP)]))
    detParams = H.det %*% params[(nSP + nPhiP+1):(nSP + nDP + nPhiP)]
    forward(detParams, phi, psi)
	
  }

  fmList <- list()
  GF <- 0
  BF <- 0
  run <- 1
  
  while(GF < min.good && BF < max.bad) {
    starts <- matrix(rnorm(nP*20), 20, nP)  # sample space initially for better start
    starts.nll <- apply(starts, 1, nll)
    start <- starts[which.min(starts.nll),]

    fmList[[run]] <- optim(start, nll, method="BFGS",hessian = TRUE,
        control=list(trace=1, maxit=400))

	GF <- GF + 1
	run <- run + 1  ## TODO: track bad fits.
  }
  
  nlls <- sapply(fmList, function(x) x$value)
  fm <- fmList[[which.min(nlls)]]

  mle <- fm$par
  psiParams <- c(0, H.psi %*% mle[1:nSP])
  psi <- exp(psiParams)
  psi <- psi/sum(psi)

  phiParams <- H.phi %*% mle[(nSP+1):(nSP + nPhiP)]
  phi = t(do.call(phiMatFun, list(p=phiParams)))
  detParams = as.vector(H.det %*% mle[(nSP + nPhiP+1):(nSP + nDP + nPhiP)])
 # fm <- genoud(nll, nP)

  ## smoothing
  forward(detParams, phi, psi, storeAlpha = TRUE)
  backward(detParams, phi, psi)
  beta[,,1]
  for(i in 1:M) {
    for(t in 1:nY) {
      gamma[,t,i] <- alpha[,t,i] * beta[,t,i] / sum(alpha[,t,i] * beta[,t,i])
    }
  }

  ## get expected detection matrices
#  mp <- array(V.itjk %*% detParams, c(nDMP, J, nY, M))
#  expectedDetMats <- .Call("getDetMats", y.arr, mp, K)
#  expectedDetMats <- array(expectedDetMats, c(K+1,K+1,J,nY,M))
  detMats <- getDetMats(detParams, phi, psi)
  
  list(nll = fm$value, phi = phi,
      psiParms = mle[1:nSP], detParms = mle[(nSP + nPhiP+1):(nSP + nDP + nPhiP)],
      phiParms = mle[(nSP+1):(nSP + nPhiP)],
      smooth = gamma,
      hessian = fm$hessian,
      convergence = fm$convergence, detMats = detMats)
}

# return (y+1)^th column of detMat
detVecLogit4 <- function(p, y) {
#  p1 = detPars[1]
#  p2 = detPars[2]
#  p3 = detPars[3]
#  p4 = detPars[4]
#  p5 = detPars[5]
#  p6 = detPars[6]

  switch(y+1,
      c(1, exp(p[1])/(1 + exp(p[1])),
          exp(p[2])/(1 + exp(p[2]) + exp(p[3])),
          exp(p[4])/(1 + exp(p[6]) + exp(p[4]) + exp(p[5]))),
      c(0, 1/(1 + exp(p[1])),
          exp(p[3])/(1 + exp(p[2]) + exp(p[3])),
          exp(p[5])/(1 + exp(p[6]) + exp(p[4]) + exp(p[5]))),
      c(0, 0, 1/(1 + exp(p[2]) + exp(p[3])),
          exp(p[6])/(1 + exp(p[6]) + exp(p[4]) + exp(p[5]))),
      c(0,0,0,1/(1 + exp(p[6]) + exp(p[4]) + exp(p[5]))))

#  matrix(c(1, exp(p[1])/(1 + exp(p[1])),
#      exp(p[2])/(1 + exp(p[2]) + exp(p[3])),
#      exp(p[4])/(1 + exp(p[6]) + exp(p[4]) + exp(p[5])),
#  0, 1/(1 + exp(p[1])),
#      exp(p[3])/(1 + exp(p[2]) + exp(p[3])),
#      exp(p[5])/(1 + exp(p[6]) + exp(p[4]) + exp(p[5])),
#  0, 0, 1/(1 + exp(p[2]) + exp(p[3])),
#      exp(p[6])/(1 + exp(p[6]) + exp(p[4]) + exp(p[5])),
#  0,0,0,1/(1 + exp(p[6]) + exp(p[4]) + exp(p[5]))),4,4)[,y+1]

}

phiLogit2 <- function(p) {
  p0 <- p[1]
  p1 <- p[2]
  matrix(plogis(c(-p0,p0,-p1,p1)),2,2)
}

detVecLogit2 <- function(p, y) {

	switch(y+1,
			c(1, exp(p[1])/(1 + exp(p[1]))),
			c(0, 1/(1 + exp(p[1]))))
}

setClass("hmmSpec",
		representation(phiMatFun = "function",
				detVecFun = "function",
				nDMP.un = "numeric",
				nPhiP.un = "numeric"
))

phiLogit4 <- function(p) {
  p0 <- p[1]
  p1 <- p[2]
  p2 <- p[3]
  p3 <- p[4]
  p4 <- p[5]
  p5 <- p[6]
  p6 <- p[7]
  p7 <- p[8]
  p8 <- p[9]
  p9 <- p[10]
  p10 <- p[11]
  p11 <- p[12]

  matrix(c(1/(1 + exp(p0) + exp(p1) + exp(p2)),
          exp(p0)/(1 + exp(p0) + exp(p1) + exp(p2)),
          exp(p1)/(1 + exp(p0) + exp(p1) + exp(p2)),
          exp(p2)/(1 + exp(p0) + exp(p1) + exp(p2)),
          1/(1 + exp(p3) + exp(p4) + exp(p5)),
          exp(p3)/(1 + exp(p3) + exp(p4) + exp(p5)),
          exp(p4)/(1 + exp(p3) + exp(p4) + exp(p5)),
          exp(p5)/(1 + exp(p3) + exp(p4) + exp(p5)),
          1/(1 + exp(p6) + exp(p7) + exp(p8)),
          exp(p6)/(1 + exp(p6) + exp(p7) + exp(p8)),
          exp(p7)/(1 + exp(p6) + exp(p7) + exp(p8)),
          exp(p8)/(1 + exp(p6) + exp(p7) + exp(p8)),
          1/(1 + exp(p9) + exp(p10) + exp(p11)),
          exp(p9)/(1 + exp(p9) + exp(p10) + exp(p11)),
          exp(p10)/(1 + exp(p9) + exp(p10) + exp(p11)),
          exp(p11)/(1 + exp(p9) + exp(p10) + exp(p11))),4,4) # TODO:  make sure this is transposed correctly
}

#' @export
profile1D <-
    function(stateformula = ~ 1, detformula = ~ 1, umf,
        detconstraint = NULL,
        phiconstraint = NULL, psiconstraint = NULL, K,
        phiMatrix = "cumlogit",
        phiMatFun, detVecFun, homo.det = FALSE,
        whichVarying, parmsIn, profSeq)
{
  nDCP <- length(attr(terms(detformula), "term.labels")) + 1
#  if(arDet)
#    nDMP <- K
#  else
    nDMP <-  K*(K+1)/2

  if(is.null(detconstraint))
    detconstraint <- matrix(rep(1:nDMP,nDCP),nDMP, nDCP)
  #detconstraint <- matrix(1:(nDMP*nDCP), nDMP, nDCP)
  #################################


  ## only do NA checking for variables being used!
  ## create reduced detection formula from detconstraint
  to.rm <- which(apply(detconstraint == 0, 2, all))
  if(length(to.rm) != 0 & length(to.rm) != (ncol(detconstraint) - 1)) {
    detformula.red <- eval(parse(text=paste("~ ",
                paste(attr(terms(detformula), "term.labels")[-(to.rm-1)],
                    collapse="+"))))
  } else if(length(to.rm) == (ncol(detconstraint) - 1)) {
    detformula.red <- ~1
  } else {
    detformula.red <- detformula
  }

  y <- umf@y
  M <- nrow(y)
  J <- umf@obsNum / umf@primaryNum
  nY <- ncol(y)/J


  designMats <- getDesign(stateformula = stateformula, detformula = detformula, umf)

  V.itj <- designMats$V
  nDCP <- ncol(V.itj)
  detParms <- colnames(V.itj)

  ## number of free terms in detection matrix.
  ## these will be modeled by XDet's
#  if(arDet)
#    nDMP <- K
#  else
#    nDMP <-  K*(K+1)/2

  if(nrow(detconstraint) != nDMP)
    stop(paste("detconstraint has wrong number of rows.\n
                It should have",nDMP))

  ## for performance, remove variables eliminated by detconstraint
  to.rm <- which(apply(detconstraint == 0, 2, all))
  if(length(to.rm) != 0) {
    detconstraint <- detconstraint[,-to.rm]
    V.itj <- V.itj[,-to.rm]
    detParms <- detParms[-to.rm]
    nDCP <- nDCP - length(to.rm)
  }
  if(is.null(dim(detconstraint)))
    detconstraint <- matrix(detconstraint,nDMP,nDCP)

  ## create design matrix for each matrix parameter
  # TODO: make this do the right thing if there are constraints.
  V.itjk <- V.itj %x%  diag(nDMP)
  #get a better line here    if(is.null(detconstraint)) detconstraint <- 1:nDMP.un

  fphiMatrix <- paste("f",phiMatrix,sep="")
  ## add "f" to fool switch because it doesn't like characters that start
  ## with numbers
  nPhiP.un.table <- matrix(c(
          1, 0, 2,
          2, 0, 6,  ## K=3, 0: logit3
          3, 0, 12,  ## 0 = full mat
          3, 1, 6,   ## 1 = small ord
          3, 2, 10), ## 2 = ord with interecept
      5, 3, byrow = TRUE)
  colnames(nPhiP.un.table) <- c("K", "phiNum", "nPhiP.un")
  nPhiP.un <- nPhiP.un.table[nPhiP.un.table[,2] == phiMatFun &
          nPhiP.un.table[,1] == K,3]
  nSP.un <- K                # number of parameters for psi vector of initial

  if(is.null(phiconstraint)) phiconstraint <- 1:nPhiP.un
  if(is.null(psiconstraint)) psiconstraint <- 1:nSP.un

  nPhiP <- switch(fphiMatrix,
      fcumlogit = K + 1,
      f4state = max(phiconstraint),
      f4stateAR = max(phiconstraint),
      f4state4 = max(phiconstraint),
      f3state = max(phiconstraint),
      f2state = max(phiconstraint),
      flogit4 = max(phiconstraint),
      flogit3 = max(phiconstraint),
      flogit2 = max(phiconstraint),
      flogit4ar = max(phiconstraint))

  nSP <- max(psiconstraint)

  ## create linked list of parameters
  theta.df <- data.frame(parameter = character(), start = numeric(),
      end = numeric(), stringsAsFactors = FALSE)

  theta.df <- addParm(theta.df, "phiParms", nPhiP)

  ## convert from new easy-input style detcon to computational format
  prev.col.max <- 0
  detcon.corr <- matrix(0, nDMP, nDCP)
  for(j in 1:nDCP) {
    detcon.corr[,j] <- ifelse(detconstraint[,j] != 0,
        detconstraint[,j] + prev.col.max, 0)
    prev.col.max <- max(c(detcon.corr[,j],prev.col.max))
  }
  detcon.vec <- as.vector(detcon.corr)
  nDP <- max(detcon.vec)
  nDP.un <- length(detcon.vec)


  H.det <- matrix(0, nDP.un, nDP)
  for(i in 1:length(detcon.vec)){
    H.det[i,detcon.vec[i]] <- 1
  }

  theta.df <- addParm(theta.df, "detParms", nDP)
  nP <- nDP + nSP + nPhiP  # total number of parameters

  ## construct constrain matrix for phi paramters
  H.phi <- matrix(0, nPhiP.un, nPhiP)
  for(i in 1:nPhiP.un){
    H.phi[i, phiconstraint[i]] <- 1
  }

  H.psi <- matrix(0, nSP.un, nSP)
  for(i in 1:nSP.un){
    H.psi[i, psiconstraint[i]] <- 1
  }

  y.itj <- as.numeric(t(y))

  ## reorder X.tjik to be X.itjk
  ###   t.tjik <- rep(1:nY, each = nDMP *M*J)
  ###   i.tjik <- rep(rep(1:M, each = nDMP), nY*J)
  ###   j.tjik <- rep(rep(1:J, each = nDMP*M), nY)
  ###   k.tjik <- rep(1:nDMP, M * J * nY)

  ###   XDet.itjk <- XDet.tjik[order(i.tjik, t.tjik, j.tjik, k.tjik),]

  ## replace NA's with 99 before passing to C++
  ## TODO: need better missing data passing mechanism (maybe NaN of Inf?)
  y.itj[is.na(y.itj)] <- 99
  V.itjk[is.na(V.itjk)] <- 9999

  # get ragged array indices
  y.it <- matrix(t(y), nY*M, J, byrow = TRUE)
  J.it <- rowSums(!is.na(y.it))


  ncX <- ifelse(is.null(ncol(V.itjk)), 1, ncol(V.itjk))
  lengthProfSeq <- length(profSeq)
  x <- .C("profile1D",
      as.integer(y.itj),
      as.double(t(V.itjk)),
      as.integer(t(J.it)),
      as.integer(ncX),
      as.integer(nDMP),
      as.integer(nDP),
      as.integer(nSP),
      as.integer(nSP.un),
      as.integer(nPhiP),
      as.integer(nP),
      as.integer(nDP.un),
      as.integer(nPhiP.un),
      as.double(H.det),
      as.double(H.phi),
      as.double(H.psi),
      as.integer(K),
      as.integer(M),
      as.integer(J),
      as.integer(nY),
      as.double(parmsIn),
      nlls = double(lengthProfSeq),
      lengthProfSeq = as.integer(lengthProfSeq),
      profSeq = as.double(profSeq),
      as.integer(whichVarying),
      as.character(phiMatrix),
      as.integer(homo.det),
      as.integer(0), # ardet
      as.integer(K)) # replace with detFUn

  return(cbind(vals = profSeq, nll = x$nlls))

}


computeFitStats.max <- function(detMats, smooth, y, J.it) {

  dims <- dim(detMats)
  K <- dims[1] - 1
  J <- dims[3]
  nY <- dims[4]
  M <- dims[5]

  ## computes expected cell probabilities for each i,t,j
  computeExpectedObs <- function(detMats, smooth) {
    expProbs <- array(NA,c(K+1,J,nY,M))
    q <- 1
    for(i in 1:M) {
      for(t in 1:nY) {
        psi.t <- smooth[,t,i]
        for(j in seq(length=J.it[q])) {
          expProbs[,j,t,i] <- as.vector(t(psi.t) %*% detMats[,,j,t,i])
        }
        q <- q + 1
      }
    }
    return(expProbs)
  }
  eo <- computeExpectedObs(detMats, smooth)

# given matrix of expected probilities for site i year t to compute the CDF of max
  maxCDF.siteyear <- function(mat) {
    mat[is.na(mat)] <- 0
    CDF.itj <- apply(mat, 2, cumsum)
    CDF.itj[CDF.itj == 0] <- 1
    CDF <- rowProds(CDF.itj)
    if(CDF[1] == 1) CDF <- rep(NA,length(CDF))  ## replace columns of all 1's with NA
    CDF
  }
  maxCDF <- apply(eo, 3:4, maxCDF.siteyear)

  ## convert a discrete CDF to a PDF
  CDFtoPDF <- function(CDF) {
    CDF2 <- c(0,CDF[-length(CDF)])
    CDF - CDF2
  }
  maxPDF <- apply(maxCDF, 2:3, CDFtoPDF)

  MaxDist <- apply(maxPDF, 1:2, sum, na.rm=TRUE)

  ## create y array (multinomial cell style)
  y.arr <- array(t(y), c(J, nY, M))#<- array(NA, c(K+1, J, nY, M))
  y.arr <- apply(y.arr, 1:3, function(x) {
        if(!is.na(x)) {
          v <- numeric(K+1)
          v[x+1] <- 1
          v
        } else {
          rep(NA, K+1)
        }
      })

  X2 <- sum((y.arr - eo)^2 / eo, na.rm = TRUE)


  return(list(MaxDist = MaxDist, eo = eo, X2 = X2))
}

computeFitStats <- function(detMats, smooth, y, J.it) {

	dims <- dim(detMats)
	K <- dims[1] - 1
	J <- dims[3]
	nY <- dims[4]
	M <- dims[5]
		
	## computes expected cell probabilities for each i,t,j
	computeExpectedObs <- function(detMats, smooth) {
		expProbs <- array(NA,c(K+1,J,nY,M))
		q <- 1
		for(i in 1:M) {
			for(t in 1:nY) {
				psi.t <- smooth[,t,i]
				for(j in seq(length=J.it[q])) {
					expProbs[,j,t,i] <- as.vector(t(psi.t) %*% detMats[,,j,t,i])
				}
				q <- q + 1
			}
		}
		return(expProbs)
	}
	eo <- computeExpectedObs(detMats, smooth)
	
	e.counts <- apply(eo, c(1,3), sum, na.rm = TRUE)
	
	y.arr <- array(t(y), c(J, nY, M))
	o.counts <- apply(y.arr, 2, table, useNA = "no")
	
	X2 <- sum((e.counts - o.counts)^2 / e.counts)
	
	list(e.counts = e.counts, o.counts = o.counts, chisq = X2)
}
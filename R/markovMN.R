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
           arDet = FALSE, phiMatFun, detVecFun, n.starts = 1, homo.det = FALSE,
           max.bad = 3, min.good = 1)
{
  ## truncate at K
  umf@y[umf@y > K] <- K

  ## coerce vector detcon for K = 1.
  if(identical(K,1)) {
    if(class(detconstraint) %in% c("numeric","integer")) {
      detconstraint <-  t(as.matrix(detconstraint))
    }
  }

  if(K == 3) {
    if(!(phiMatrix %in%
         c("4state", "cumlogit", "4state4", "4stateAR","logit4","logit4ar")))
      stop(paste("Inappropriate phiMatrix specified for K =",K))
  }
  if(K == 2) {
    if(!(phiMatrix %in% c("3state", "cumlogit", "logit3")))
      stop(paste("Inappropriate phiMatrix specified for K =",K))
  }

  ##################################
  ## section determines appropriate default detconstraint...
  ## needs more investigation.
  nDCP <- length(attr(terms(detformula), "term.labels")) + 1
  if(arDet)
    nDMP <- K
  else
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
  fc$max.bad <- as.name("max.bad")
  fc$min.good <- as.name("max.bad")
  fc$homo.det <- as.name("homo.det")
  fm <- eval(fc)

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

  ## TODO: add checking of hessian for NaN's, Inf's, 0's, etc.
  ## Their presence should trigger "convergence <- 1" catching.


  return(fm)
}


markovMN.fit <- function(stateformula = ~ 1, detformula = ~ 1, umf,
                         detconstraint = NULL,
                         phiconstraint = NULL, psiconstraint = NULL, J, K,
                         phiMatrix = "cumlogit", EM = FALSE, psi.init = NULL,
                         phiParms.init = NULL, detParms.init = NULL,
                         get.inits = TRUE, trace = FALSE,
                         arDet = FALSE, phiMatFun = NULL, detVecFun = NULL,
                         n.starts = 3, homo.det = FALSE, max.bad, min.good)
{
  y <- umf@y

  M <- nrow(y)
  nY <- ncol(y)/J

  designMats <- getDesign(stateformula = stateformula, detformula = detformula, umf)

  V.itj <- designMats$V
  nDCP <- ncol(V.itj)
  detParms <- colnames(V.itj)

  ## number of free terms in detection matrix.
  ## these will be modeled by XDet's
  if(arDet)
    nDMP <- K
  else
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
  nPhiP.un <- switch(fphiMatrix,
                     fcumlogit = K + 1,
                     f4state = 12,
                     f4state4 = 4,
                     f4stateAR = 10,
                     f3state = 6,
                     f2state = 2,
                     flogit4 = 12,
                     flogit3 = 6,
                     flogit2 = 2,
                     flogit4ar = 6)
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

  if(is.null(psi.init)) psi.init <- rnorm(nSP)
  if(is.null(phiParms.init)) phiParms.init <- rnorm(nPhiP)
  if(is.null(detParms.init)) detParms.init <- rnorm(nDP)

  fm <- findMLE(y.itj, V.itjk, nDMP, nDP, nSP, nSP.un,
                   nPhiP, nP, nDP.un,
                   nPhiP.un, H.det, H.phi, H.psi, K, M, J, nY,
                   phiMatrix, EM= EM, psi.init = psi.init,
                   phiParms.init = phiParms.init,
                   detParms.init = detParms.init, get.inits = get.inits,
                   trace = trace, arDet = arDet,
                   phiMatFun = phiMatFun, detVecFun = detVecFun, n.starts = n.starts,
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
       arDet = arDet,
       detform = detformula,
       psiconstraint = psiconstraint,
       phiconstraint = phiconstraint,
       detconstraint = detconstraint,
       phiMatrix = phiMatrix, K = K, hessian = hessian,
       n = M, convergence = fm$convergence)
}


findMLE <-
  function(y.itj, V.itjk, nDMP, nDP, nSP, nSP.un, nPhiP, nP, nDP.un,
           nPhiP.un, H.det, H.phi, H.psi, K, M, J, nY,
           phiMatrix, smooth.only = FALSE, EM,
           psi.init, phiParms.init, detParms.init, get.inits, trace,
           arDet,phiMatFun = phiMatFun, detVecFun = detVecFun,n.starts,homo.det,
           max.bad,min.good)
{
  ncX <- ifelse(is.null(ncol(V.itjk)), 1, ncol(V.itjk))
  x <- .C("findMLE",
          as.integer(y.itj),
          as.double(t(V.itjk)),
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
          nll = numeric(1),
          phi = numeric((K + 1)^2),
          psiParms = double(nSP),
          detParms = double(nDP),
          phiParms = double(nPhiP),
          smooth = numeric(M*nY*(K + 1)),
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
          convergence = integer(1))
#          as.integer(max.bad),
#          as.integer(min.good))

  return(list(nll = x$nll, phi = matrix(x$phi,K + 1, K + 1),
              psiParms = x$psiParms, detParms = x$detParms,
              phiParms = x$phiParms,
              smooth = array(x$smooth, c(K + 1, nY, M)),
              hessian = matrix(x$hessian,nP,nP),
              convergence = x$convergence))
}

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
    function(y.itj, V.itjk, nDMP, nDP, nSP, nSP.un, nPhiP, nP, nDP.un,
        nPhiP.un, H.det, H.phi, H.psi, K, M, J, nY,
        phiMatrix, smooth.only = FALSE, EM,
        psi.init, phiParms.init, detParms.init, get.inits, trace,
        arDet, phiMatFun, detVecFun, n.starts = 3, homo.det = FALSE,
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
  forward <- function(detParms, phi, psi, storeAlpha = FALSE) {

    negloglike <- 0
##################    version to get all detvec's at once
    psiSite <- matrix(psi, K + 1, M)
    detVec.arr <- array(1, c(K + 1, nY, M))
    if(homo.det) {
      detVecMat <- matrix(NA,K+1,K+1)
      for(i in 1:(K+1)) {
        detVecMat[,i] <- do.call(detVecFun, list(p = detParms, y = i-1))
      }
    } else {
      mp <- array(V.itjk %*% detParms, c(nDMP, J, nY, M))
    }
    for(t in 1:nY) {
      for(i in 1:M) {
        for(j in 1:J) {
          if(y.arr[i,t,j] != 99) {
            if(!homo.det) {
              detVec.arr[,t,i] <- detVec.arr[,t,i] *
                  do.call(detVecFun, list(p=mp[,j,t,i], y=y.arr[i,t,j]))
            } else {
              detVec.arr[,t,i] <- detVec.arr[,t,i] * detVecMat[,y.arr[i,t,j]+1]
            }
          }
        }
      }
      psiSite <- psiSite * detVec.arr[,t,]
      if(storeAlpha) alpha[,t,] <<- psiSite[,]
      if(t < nY) {
        psiSite <- phi %*% psiSite
      } else {
        negloglike <- negloglike - sum(log(colSums(psiSite)))
      }
    }




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
    negloglike
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
            detVecObs <- do.call(detVecFun, list(p=mp, y=y.arr[i,t,j]))
            detVec <- detVec * detVecObs
          }
        }

        backP <- t(phi) %*% (detVec * backP)

      }
    }
  }

  nll <- function(params) {

    psiParams <- c(0, H.psi %*% params[1:nSP])
    psi <- exp(psiParams)
    psi <- psi/sum(psi)

    phi = do.call(phiMatFun, list(p=H.phi %*% params[(nSP+1):(nSP + nPhiP)]))
    detParams = H.det %*% params[(nSP + nPhiP+1):(nSP + nDP + nPhiP)]
    forward(detParams, phi, psi)
  }

  fmList <- list()
  for(run in seq(length=n.starts)) {
    starts <- matrix(rnorm(nP*20), 20, nP)  # sample space initially for better start
    starts.nll <- apply(starts, 1, nll)
    start <- starts[which.min(starts.nll),]

    fmList[[run]] <- optim(starts, nll, method="BFGS",hessian = TRUE,
        control=list(trace=6, maxit=300, REPORT=1))
  }
  nlls <- sapply(fmList, function(x) x$nll)
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
  for(i in 1:M) {
    for(t in 1:nY) {
      gamma[,t,i] <- alpha[,t,i] * beta[,t,i] / sum(alpha[,t,i] * beta[,t,i])
    }
  }

  list(nll = fm$value, phi = phi,
      psiParms = mle[1:nSP], detParms = mle[(nSP + nPhiP+1):(nSP + nDP + nPhiP)],
      phiParms = mle[(nSP+1):(nSP + nPhiP)],
      smooth = gamma,
      hessian = fm$hessian,
      convergence = fm$convergence)
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
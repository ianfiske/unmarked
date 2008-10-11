#' @include classes.R
roxygen()

#' Fit the Occupancy model of Royle and Nichols
#'
#' This function fits the standard occupancy model of 
#'
#' @param stateformula Right-hand side formula describing covariates of abundance
#' @param detformula Right-hand side formula describing covariates of detection
#' @param umf unMarkedFrame supplying data to the model.
#' @return unMarkedFit object describing the model fit.
#' @references Royle and Nichols (2003)
#' @examples
#' data(birds)
#' woodthrushUMF <- unMarkedFrame(woodthrush.bin)
#' fm.wood.rn <- occuRN(~1, ~1, y = woodthrush.bin)
#' data(frogs)
#' pcruUM <- unMarkedFrame(pcru.bin)
#' fm <- occu(~1, ~1, pcruUM)
#' @export
occuRN <- 
function(stateformula, detformula, umf)
{

  umf <- handleNA(stateformula, detformula, umf)
  designMats <- getDesign(stateformula, detformula, umf)
  X <- designMats$X; V <- designMats$V
  y <- umf@y
  J <- ncol(y)
  M <- nrow(y)
  K <- 20
  
  occParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nOP <- ncol(X)

  nP <- nDP + nOP
  y.ji <- as.numeric(y)
  y.jik <- rep(y.ji, each = K + 1)
  navec <- is.na(y.jik)
  nd <- ifelse(rowSums(y, na.rm=TRUE) == 0, 1, 0) # no det site indicator
  k <- 0:K
  k.i <- rep(k, M)
  k.ji <- rep(k, M * J)
  
  nll <- function(parms, f = "Poisson")
  {    
    lam.i <- exp(X %*% parms[1 : nOP])
    lam.ik <- rep(lam.i, each = K + 1)
    r.ji <- plogis(V %*% parms[(nOP + 1) : nP])

    r.jik <- rep(r.ji, each = K + 1)
    p.sup <- 1 - (1 - r.jik)^(k.ji)
    cp <- p.sup^y.jik * (1 - p.sup)^(1 - y.jik)
    cp[navec] <- 1
    cp.mat <- matrix(cp, M * (K + 1), J)
    p.ik <- rowProds(cp.mat)

    f.ik <- dpois(k.i, lam.ik)
    dens.mat <- matrix(p.ik * f.ik, M, K + 1, byrow = TRUE)
    dens.integ <- rowSums(dens.mat)
  
    -sum(log(dens.integ))      # multiply likelihood over all sites
  }

  fm <- optim(rep(0, nP), nll, method = "BFGS", hessian = TRUE)
  ests.se <- diag(solve(fm$hessian))
  ests <- fm$par
  fm <- optim(rep(0,nP), nll, method = "BFGS")
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP + 2 * nP * (nP + 1) / (M - nP - 1)
  names(ests) <- c(occParms, detParms)
  names(ests.se) <- c(occParms, detParms)
  umfit <- unMarkedFit(fitType = "occuRN",
                       stateformula = stateformula, detformula = detformula,
                       data = umf, stateMLE = ests[1:nOP],
                       stateSE = ests.se[1:nOP], 
                       detMLE = ests[(nOP + 1) : nP],
                       detSE = ests.se[(nOP + 1): nP], AIC = fmAIC)
  return(umfit)
}

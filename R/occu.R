#' @include classes.R
roxygen()

#' Fit the MacKenzie Occupancy Model
#'
#' This function fits the standard occupancy model of 
#'
#' @param stateformula Right-hand side formula describing covariates of abundance
#' @param detformula Right-hand side formula describing covariates of detection
#' @return unMarkedFit object describing the model fit.
#' @references MacKenzie 2002, Occupancy.
#' @examples
#' data(frogs)
#' pcruUM <- unMarkedFrame(pcru.bin)
#' fm <- occu(~1, ~1, pcruUM)
#' @export
occu <-
function(stateformula, detformula, umf)
{
  umf <- handleNA(stateformula, detformula, umf)
  designMats <- getDesign(stateformula, detformula, umf)
  X <- designMats$X; V <- designMats$V
  y <- umf@y
  J <- ncol(y)
  M <- nrow(y)

  occParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nOP <- ncol(X)

  nP <- nDP + nOP
  yvec <- as.numeric(y)
  navec <- is.na(yvec)
  nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i indicator
  
  nll <- function(parms) {
    psi <- plogis(X %*% parms[1 : nOP])
    pvec <- plogis(V %*% parms[(nOP + 1) : nP])
    cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
    cp[navec] <- 1  # so that NA's don't modify likelihood        
    cpmat <- matrix(cp, M, J) # put back into matrix to multiply appropriately
    loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi)) 
    -sum(loglik)
  }
  
  fm <- optim(rep(0, nP), nll, method = "BFGS", hessian = TRUE)
  ests.se <- diag(solve(fm$hessian))
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP
  names(ests) <- c(occParms, detParms)
  names(ests.se) <- c(occParms, detParms)
  umfit <- unMarkedFit(fitType = "occu",
                       stateformula = stateformula, detformula = detformula,
                       data = umf, stateMLE = ests[1:nOP],
                       stateSE = ests.se[1:nOP], 
                       detMLE = ests[(nOP + 1) : nP],
                       detSE = ests.se[(nOP + 1): nP], AIC = fmAIC)
  return(umfit)
}



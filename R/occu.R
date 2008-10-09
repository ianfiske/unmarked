#' @include classes.R
roxygen()

#' Fit the MacKenzie Occupancy Model
#'
#' This function fits the standard occupancy model of 
#'
#' @param stateformula Right-hand side formula describing covariates of abundance
#' @param detformula Right-hand side formula describing covariates of detection
#' @return unMarkedFit object describing the model fit.
#' @export
occu <-
function(stateformula, detformula, umf)
{

  umf <- handleNA(stateformula, detformula, umf)
  y <- umf@y
  ## Compute detection design matrix
  ## add site Covariates at observation-level
  V.mf <- model.frame(detformula, umf@obsCovs)
  V <- model.matrix(detformula, V.mf) 
  ## Compute state design matrix
  X.mf <- model.frame(stateformula, umf@siteCovs)
  X <- model.matrix(stateformula, X.mf)

  occParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nOP <- ncol(X)
  J <- ncol(y)
  M <- nrow(y)

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
  
  fm <- optim(rep(0, nP), nll, method = "BFGS")
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP
  names(ests) <- c(occParms, detParms)
  list(estimates = ests, AIC = fmAIC)
}



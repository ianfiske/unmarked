#' @include classes.R
roxygen()

## Take an M x J matrix of detection probabilities and return a matrix
## of M x K observation probs
removalPiFun <- function(p){
  M <- nrow(p)
  J <- ncol(p)
  pi <- matrix(NA, M, J)
  pi[,1] <- p[,1]
  for(i in seq(from = 2, length = J - 1)) {
    pi[, i] <- pi[, i-1]*(1 - p[,i-1])
  }
  return(pi)
}

# p is an M x 2 matrix of detection probabilities (site x observer).
# return
doublePiFun <- function(p){
  M <- nrow(p)
  pi <- matrix(NA, M, 3)
  pi[,1] <- p[,1] * (1 - p[,2])
  pi[,2] <- p[,2] * (1 - p[,2])
  pi[,3] <- p[,1] * p[,2]
  return(pi)
}

#' Fit the multinomial-Poisson abundance mixture model
#'
#' This function takes advantage of the closed form of the integrated
#' likelihood when a latent Poisson distribution is assumed for abundance
#' at each site and a multinomial distribution 
#'
#' @param stateformula Right-hand side formula describing covariates of abundance
#' @param detformula Right-hand side formula describing covariates of detection
#' @param piFun Function to define multinomial cell probabilities.
#' @return unMarkedFit object describing the model fit.
#' @export
multinomPois <- 
function(stateformula, detformula, piFun, umf)
{

  umf <- handleNA(stateformula, detformula, umf)
  y <- umf@y
  ## Compute detection design matrix
  V.mf <- model.frame(detformula, umf@obsCovs)
  V <- model.matrix(detformula, V.mf)
 
  ## Compute state design matrix
  X.mf <- model.frame(stateformula, umf@siteCovs)
  X <- model.matrix(stateformula, X.mf)

###   cleaned <- handleNA(arranged, stateformula, detformula)

  lamParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nAP <- ncol(X)

  J <- ncol(y)
  M <- nrow(y)

  nP <- nDP + nAP
  yvec <- as.numeric(y)
  navec <- is.na(yvec)

  nll <- function(parms) {
    lambda <- exp(X %*% parms[1 : nAP])
    p <- plogis(V %*% parms[(nAP + 1) : nP])
    p.matrix <- matrix(p, M, umf@obsNum, byrow = TRUE)
    pi <- do.call(piFun, list(p = p.matrix))
    logLikeSite <- dpois(y, matrix(lambda, M, J) * pi, log = TRUE)
    logLikeSite[navec] <- 0
    -sum(logLikeSite)
  }
  
  fm <- optim(rep(0, nP), nll, method = "BFGS", hessian = TRUE)
  ests.se <- diag(solve(fm$hessian))
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP
  names(ests) <- c(lamParms, detParms)
  names(ests.se) <- c(lamParms, detParms)
  umfit <- unMarkedFit(fitType = "multinomPois",
                       stateformula = stateformula, detformula = detformula,
                       data = umf, stateMLE = ests[1:nAP],
                       stateSE = ests.se[1:nAP], 
                       detMLE = ests[(nAP + 1) : nP],
                       detSE = ests.se[(nAP + 1): nP], AIC = fmAIC)
}

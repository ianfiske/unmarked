#' @include classes.R
#' @include utils.R
roxygen()

## Take an M x J matrix of detection probabilities and return a matrix
## of M x J observation probs
#' @export
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
# returns an M x 3 matrix of row=(1 not 2, 2 not 1, 1 and 2).
#' @export
doublePiFun <- function(p){
  M <- nrow(p)
  pi <- matrix(NA, M, 3)
  pi[,1] <- p[,1] * (1 - p[,2])
  pi[,2] <- p[,2] * (1 - p[,2])
  pi[,3] <- p[,1] * p[,2]
  return(pi)
}


#' Fit the multinomial-Poisson abundance mixture model.
#'
#' This function takes advantage of the closed form of the integrated
#' likelihood when a latent Poisson distribution is assumed for abundance
#' at each site and a multinomial distribution is taken for the observation
#' state. Many common sampling methods can be framed in this context.  For
#' example, double-observer point counts, removal sampling, and distance
#' sampling can all be analyzed with this function by specifying the proper
#' multinomial cell probablilities.  This is done with by supplying the
#' appropriate function (piFun) argument.  removalPiFun() and doublePiFun()
#' are supplied as example cell probability functions.
#'
#' @title Multinomial-Poisson Mixtures
#' @param stateformula Right-hand side formula describing covariates of abundance
#' @param detformula Right-hand side formula describing covariates of detection
#' @param piFun Function to define multinomial cell probabilities.
#' @param umf unMarkedFrame supplying data.
#' @return unMarkedFit object describing the model fit.
#' @author Ian Fiske
#' @keywords models
#' @references
#' Royle, J., Dawson, D., & Bates, S. (2004). Modeling abundance effects in distance sampling. Ecology, 85(6), 1591-1597.
#'
#' Royle, J. A. (2004). Generalized estimators of avian abundance from count survey data. Animal Biodiversity and Conservation, 27(1), 375-386.
#'
#' Royle, J. A., & Dorazio, R. M. (2006). Hierarchical Models of Animal Abundance and Occurrence. Journal Of Agricultural Biological And Environmental Statistics, 11(3), 249.
#' @examples
#' data(ovendata)
#' ovenFrame <- unMarkedFrame(ovendata.list$data,
#'                            siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])))
#' fm1 <- multinomPois(~ ufp + trba, ~ 1, removalPiFun, ovenFrame)
#' fm2 <- multinomPois(~ ufp, ~1, removalPiFun, ovenFrame)
#' fm3 <- multinomPois(~ trba, ~1, removalPiFun, ovenFrame)
#' fm4 <- multinomPois(~ 1, ~1, removalPiFun, ovenFrame)
#' fm4
#' fm1
#' @export
multinomPois <-
function(stateformula, detformula, piFun, umf)
{

  umf <- handleNA(stateformula, detformula, umf)
  designMats <- getDesign(stateformula, detformula, umf)
  X <- designMats$X; V <- designMats$V
  y <- umf@y
  J <- ncol(y)
  M <- nrow(y)

  lamParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nAP <- ncol(X)
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
  tryCatch(covMat <- solve(fm$hessian),
      error=simpleError("Hessian is not invertible.  Try using fewer covariates."))
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP
  names(ests) <- c(lamParms, detParms)

  stateEstimates <- unMarkedEstimate(name = "Abundance",
      estimates = ests[1:nAP],
      covMat = as.matrix(covMat[1:nAP,1:nAP]), invlink = "exp",
      invlinkGrad = "exp")

  detEstimates <- unMarkedEstimate(name = "Detection",
      estimates = ests[(nAP + 1) : nP],
      covMat = as.matrix(covMat[(nAP + 1) : nP, (nAP + 1) : nP]), invlink = "logistic",
      invlinkGrad = "logistic.grad")

  umfit <- unMarkedFit(fitType = "multinomPois",
                       stateformula = stateformula, detformula = detformula,
                       data = umf, stateEstimates = stateEstimates,
                       detEstimates = detEstimates, AIC = fmAIC,
                       hessian = fm$hessian)
  return(umfit)
}

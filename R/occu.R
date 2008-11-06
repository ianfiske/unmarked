#' @include classes.R
#' @include utils.R
roxygen()

#' This function fits the standard occupancy model of MacKenzie et al
#'
#'  See \link{unmarked-package} for detailed descriptions of passing data \code{y},
#'  \code{covdata.site}, and \code{covdata.obs}, and specifying covariates
#'  with \code{stateformula} and \code{detformula}.
#'
#'  \command{occu} fits the traditional occupancy model based on the
#'  binomial mixture models (MacKenzie et al. 2006, Royle and Dorazio
#'  2008).  The occupancy state process (\eqn{z_i}) of site \eqn{i} is
#'  modeled as
#'
#'  \deqn{z_i \sim Bernoulli(\psi_i)}{z_i ~ Bernoulli(psi_i)}
#'
#'  The observation process is modeled as
#'
#'  \deqn{y_{ij} \sim Bernoulli(p_{ij})}{y_ij ~ Bernoulli(p_ij)}
#'
#'  Covariates of \eqn{\psi_i}{psi_i} and \eqn{p_{ij}}{p_ij} are modelled
#'  using the logit link.
#' @title Fit the MacKenzie Occupancy Model
#' @param stateformula right-hand side formula describing covariates of occurence.
#' @param detformula right-hand side formula describing covariates of detection.
#' @param umf unMarkedFrame object that supplies the data (see \link{unMarkedFrame})..
#' @param knownOcc vector of sites that are known to be occupied.
#' @return unMarkedFit object describing the model fit.
#' @references
#' MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege, J. Andrew Royle, and C. A. Langtimm. “Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One.” Ecology 83, no. 8 (2002): 2248-2255.
#' MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.  Amsterdam: Academic Press.  Royle, J. A. and R. Dorazio. (2008) \emph{Book Name}.
#' @author Ian Fiske
#' @examples
#' data(frogs)
#' pcruUM <- unMarkedFrame(pcru.bin)
#' fm <- occu(~1, ~1, pcruUM)
#' fm
#' @keywords models
#' @export
occu <-
function(stateformula, detformula, umf, knownOcc = numeric(0))
{
  umf <- handleNA(stateformula, detformula, umf)
  designMats <- getDesign(stateformula, detformula, umf)
  X <- designMats$X; V <- designMats$V
  y <- truncateToBinary(umf@y)
  J <- ncol(y)
  M <- nrow(y)

  if(nrow(umf@y) != M & length(knownOcc) > 0)
    stop("sites dropped, but knownOcc was specified.")
  
  occParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nOP <- ncol(X)

  nP <- nDP + nOP
  yvec <- as.numeric(t(y))
  navec <- is.na(yvec)
  nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i indicator

  nll <- function(parms) {
    psi <- plogis(X %*% parms[1 : nOP])
    psi[knownOcc] <- 1
    pvec <- plogis(V %*% parms[(nOP + 1) : nP])
    cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
    cp[navec] <- 1  # so that NA's don't modify likelihood        
    cpmat <- matrix(cp, M, J) # put back into matrix to multiply appropriately
    loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi)) 
    -sum(loglik)
  }
  
  fm <- optim(rnorm(nP), nll, method = "BFGS", hessian = TRUE)
  ests.se <- diag(solve(fm$hessian))
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP + 2*nP*(nP + 1)/(M - nP - 1)
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

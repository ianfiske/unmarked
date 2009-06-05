#' @include unmarkedFit.R
#' @include unmarkedEstimate.R
#' @include utils.R
roxygen()

#'  This function estimates the standard occupancy model.
#'
#'  See \link{unmarkedFrame} for a description of how to supply data to the \command{umf}
#'  argument.
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
#' @param umf unmarkedFrame object that supplies the data (see \link{unmarkedFrame})..
#' @param knownOcc vector of sites that are known to be occupied.
#' @return unmarkedFit object describing the model fit.
#' @references
#' MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege, J. Andrew Royle, and C. A. Langtimm. Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. Ecology 83, no. 8 (2002): 2248-2255.
#' MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.  Amsterdam: Academic Press.  Royle, J. A. and R. Dorazio. (2008) \emph{Book Name}.
#' @author Ian Fiske
#' @examples
#' data(frogs)
#' pferUMF <- unmarkedFrame(pfer.bin)
#' fm <- occu(~1, ~1, pferUMF)
#' fm
#' @keywords models
#' @export
occu <-
function(stateformula, detformula, umf, knownOcc = numeric(0), profile = FALSE)
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

  nll <- function(params) {
    psi <- plogis(X %*% params[1 : nOP])
    psi[knownOcc] <- 1
    pvec <- plogis(V %*% params[(nOP + 1) : nP])
    cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
    cp[navec] <- 1  # so that NA's don't modify likelihood
    cpmat <- matrix(cp, M, J, byrow = TRUE) # put back into matrix to multiply appropriately
    loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
    -sum(loglik)
  }

  fm <- optim(rnorm(nP), nll, method = "BFGS", hessian = TRUE)
  tryCatch(covMat <- solve(fm$hessian),
      error=simpleError("Hessian is not invertible.  Try using fewer covariates."))
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP + 2*nP*(nP + 1)/(M - nP - 1)
  names(ests) <- c(occParms, detParms)

  if(profile) {
    profile.matrix <- matrix(NA, nP, 2)
    for(i in seq(length=nP)) {
      profile.matrix[i,] <- profileCI(nll, i, ests, c(-10,10))
          #c(ests[i] - 20*ests.se[i], ests[i] + 20*ests.se[i]))
    }
  }

  state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
      estimates = ests[1:nOP],
      covMat = as.matrix(covMat[1:nOP,1:nOP]), invlink = "logistic",
      invlinkGrad = "logistic.grad")

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
      estimates = ests[(nOP + 1) : nP],
      covMat = as.matrix(covMat[(nOP + 1) : nP, (nOP + 1) : nP]), invlink = "logistic",
      invlinkGrad = "logistic.grad")

  estimateList <- unmarkedEstimateList(list(state=state, det=det))

  umfit <- unmarkedFit(fitType = "occu",
      call = match.call(), data = umf, estimates = estimateList,
      AIC = fmAIC, hessian = fm$hessian, negLogLike = fm$value)

  return(umfit)
}

#' @include classes.R
#' @include utils.R
roxygen()

#' Fit the N-mixture point count model
#'
#' This function fits binomial-Poisson mixture model for spatially replicated point count data.
#'
#'  See \code{\link{unMarkedFrame}} for a description of how to supply by creating
#'  and unMarkedFrame.
#'
#'  This function fits the latent N-mixture model for point count data
#'  (Royle 2004, Kery and Royle 2005).
#'
#'  The latent abundance distribution, \eqn{f(N | \mathbf{\theta})}{f(N |
#'  theta)} can be set as either a Poisson or a negative binomial random
#'  variable, depending on the setting of the \code{mixture} argument.
#'  \code{mixture = "P"} or \code{mixture = "NB"} select the Poisson or
#'  negative binomial distribution respectively.  The mean of \eqn{N_i} is
#'  \eqn{\lambda_i}{lambda_i}.  If \eqn{N_i \sim NB}{N_i ~ NB}, then an
#'  additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
#'  \eqn{\alpha}{alpha} implies higher variance).
#'
#'  The detection process is modeled as binomial: \eqn{y_{ij} \sim
#'  Binomial(N_i, p_{ij})}{y_ij ~ Binomial(N_i, p_ij)}.
#'
#'  Covariates of \eqn{\lambda_i}{lamdba_i} use the log link and
#'  covariates of \eqn{p_{ij}}{p_ij} use the logit link.
#'
#' @param stateformula Right-hand side formula describing covariates of abundance
#' @param detformula Right-hand side formula describing covariates of detection
#' @param umf an unMarkedFrame supplying data to the model.
#' @param K Integer upper index of integration for N-mixture.
#' @param mixture character specifying mixture: either "P" or "NB".
#' @return unmarkedFit object describing the model fit.
#' @author Ian Fiske \email{ianfiske@@gmail.com}
#' @references Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.
#'
#' Kery, M. and Royle, J. A. (2005) Modeling Avaian Abundance from Replicated Counts Using Binomial Mixture Models. \emph{Ecological Applications} 15(4), pp. 1450--1461.
#' @examples
#' data(mallard)
#' mallardUMF <- unMarkedFrame(mallard.y, siteCovs = mallard.site,
#'                            obsCovs = mallard.obs)
#' fm.mallard <- pcount(~ length + elev + forest, ~ ivel+ date + I(date^2),
#'                   mallardUMF)
#' fm.mallard
#' @export
#' @keywords models
pcount <-
function(stateformula, detformula, umf, K = NULL, mixture = "P")
{
  if ((mixture %in% c("P","NB")) == FALSE)
    stop("Mixture familiy not recognized. Please choose 'P' or 'NB'.")

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

  if(is.null(K)) K <- max(y, na.rm = TRUE) + 20
  if(K <= max(y, na.rm = TRUE))
    stop("specified K is too small. Try a value larger than any observation")
  k <- 0:K
  M <- nrow(y)
  J <- ncol(y)
  k.ik <- rep(k, M)
  k.ijk <- rep(k, M*J)

  nP <- nAP + nDP + ifelse(identical(mixture,"NB"),1,0)
  y.ij <- as.numeric(t(y))
  y.ijk <- rep(y.ij, each = K + 1)
  navec <- is.na(y.ijk)
  nd <- ifelse(rowSums(y, na.rm=TRUE) == 0, 1, 0) # I(no detection at site i)

  nll <- function(parms){

    theta.i <- exp(X %*% parms[1 : nAP])
    p.ij <- plogis(V %*% parms[(nAP + 1) : (nAP + nDP)])
    theta.ik <- rep(theta.i, each = K + 1)
    p.ijk <- rep(p.ij, each = K + 1)

    bin.ijk <- dbinom(y.ijk,k.ijk,p.ijk)
    bin.ijk[which(is.na(bin.ijk))] <- 1
    bin.ik.mat <- matrix(bin.ijk, M * (K + 1), J)
    g.ik <- rowProds(bin.ik.mat)

    if(identical(mixture,"P")) {
      f.ik <- dpois(k.ik,theta.ik)
    }
    else if (identical(mixture,"NB")){
      f.ik <- dnbinom(k.ik, mu = theta.ik, size = exp(parms[nP]))
    }
    dens.i.mat <- matrix(f.ik * g.ik, M, K + 1, byrow = TRUE)
    dens.i <- rowSums(dens.i.mat)  # sum over the K

    -sum(log(dens.i))
  }

#  if(identical(mixture,"P")) {
    fm <- optim(rep(0,nP), nll, method="BFGS", hessian = TRUE)
#  }
#  else if (identical(mixture,"NB")){
#    fm <- optim(c(rep(0, nP),1), nll, method="BFGS")
#  }

  ests <- fm$par
  if(identical(mixture,"NB"))
     {nbParm <- "alpha"}
  else
     {nbParm <- character(0)}
  names(ests) <- c(lamParms, detParms, nbParm)
  tryCatch(covMat <- solve(fm$hessian),
      error=simpleError("Hessian is not invertible.  Try using fewer covariates."))
  fmAIC <- 2 * fm$value + 2 * nP

  stateEstimates <- unmarkedEstimate(name = "Abundance",
      estimates = ests[1:nAP],
      covMat = as.matrix(covMat[1:nAP,1:nAP]), invlink = "exp",
      invlinkGrad = "exp")

  detEstimates <- unmarkedEstimate(name = "Detection",
      estimates = ests[(nAP + 1) : nP],
      covMat = as.matrix(covMat[(nAP + 1) : nP, (nAP + 1) : nP]), invlink = "logistic",
      invlinkGrad = "logistic.grad")

  estimateList <- unmarkedEstimateList(list(state=stateEstimates, det=detEstimates))

#  umfit <- unmarkedFit(fitType = "pcount",
#      stateformula = stateformula, detformula = detformula,
#      data = umf, stateEstimates = stateEstimates,
#      detEstimates = detEstimates, AIC = fmAIC, hessian = fm$hessian)

  umfit <- unmarkedFit(fitType = "pcount",
      call = match.call(), data = umf, estimates = estimateList,
      AIC = fmAIC, hessian = fm$hessian)

  return(umfit)
}

#' @include classes.R
#' @include utils.R
roxygen()

#' Fit the Occupancy model of Royle and Nichols
#'
#'  See \link{unmarked} for detailed descriptions of passing data \code{y},
#'  \code{covdata.site}, and \code{covdata.obs}, and specifying covariates
#'  with \code{stateformula} and \code{detformula}.
#'
#'  This function fits the latent abundance mixture model described in
#'  Royle and Nichols (2003).
#'
#'  The latent abundance of site \eqn{i} is modelled as Poisson:
#'
#'  \deqn{N_i \sim Poisson(\lambda_i)}{N_i ~ Poisson(lambda_i)}
#'
#'  The detection of a single individual in site \eqn{i} during sample
#'  \eqn{j} is modelled as Bernoulli:
#'
#'  \deqn{w_{ij} \sim Bernoulli(r_{ij})}{w_ij ~ Bernoulli(r_ij)}.
#'
#'  Thus, the detection probability for a single site is linked to the
#'  detection probability for an individual by
#'
#'  \deqn{p_{ij} = 1 - (1 - r_{ij}) ^ {N_i}}{p_ij = 1 - (1 - r_ij) ^ N_i}
#'
#'  Covariates of \eqn{\lambda_i}{lambda_i} are modelled with the log link
#'  and covariates of \eqn{r_{ij}}{r_ij} are modelled with the logit link.
#'
#' @param stateformula Right-hand side formula describing covariates of abundance
#' @param detformula Right-hand side formula describing covariates of detection
#' @param umf unMarkedFrame supplying data to the model.
#' @return unMarkedFit object describing the model fit.
#' @author Ian Fiske
#' @references
#' Royle, J. A. and Nichols, J. D. (2003) Estimating Abundance from
#' Repeated Presence-Absence Data or Point Counts. \emph{Ecology}, 84(3)
#' pp. 777--790.
#' @examples
#' data(birds)
#' woodthrushUMF <- unMarkedFrame(woodthrush.bin)
#' fm.wood.rn <- occuRN(~1, ~obs, woodthrushUMF)
#' fm.wood.rn
#' @keywords models
#' @export
occuRN <-
function(stateformula, detformula, umf)
{
  umf <- handleNA(stateformula, detformula, umf)
  designMats <- getDesign(stateformula, detformula, umf)
  X <- designMats$X; V <- designMats$V
  y <- truncateToBinary(umf@y)

  J <- ncol(y)
  M <- nrow(y)
  K <- 25

  occParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nOP <- ncol(X)

  nP <- nDP + nOP
  y.ji <- as.vector(y)
  navec <- is.na(y.ji)
  n <- 0:K

  nll <- function(parms, f = "Poisson")
  {

    ## compute individual level detection probabilities
    r.ij <- matrix(plogis(V %*% parms[(nOP + 1) : nP]), M, J, byrow = TRUE)

    ## compute list of detection probabilities along N
    p.ij.list <- lapply(n, function(k) 1 - (1 - r.ij)^k)

    ## compute P(y_{ij} | N) (cell probabilities) along N
    cp.ij.list <- lapply(p.ij.list, function(pmat) pmat^y * (1-pmat)^(1-y))

    ## replace NA cell probabilities with 1.
    cp.ij.list <- lapply(cp.ij.list, function(cpmat) {
      cpmat[navec] <- 1
      cpmat
    })

    ## multiply across J to get P(y_i | N) along N
    cp.in <- sapply(cp.ij.list, rowProds)

    ## compute P(N = n | lambda_i) along i
    lambda.i <- exp(X %*% parms[1 : nOP])
    lambda.in <- sapply(n, function(x) dpois(x, lambda.i))

    ## integrate over P(y_i | N = n) * P(N = n | lambda_i) wrt n
    like.i <- rowSums(cp.in * lambda.in)

    -sum(log(like.i))
  }

  fm <- optim(rep(0, nP), nll, method = "BFGS", hessian = TRUE)
  tryCatch(covMat <- solve(fm$hessian),
      error=simpleError("Hessian is not invertible.  Try using fewer covariates."))
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP + 2 * nP * (nP + 1) / (M - nP - 1)
  names(ests) <- c(occParms, detParms)

  stateEstimates <- unMarkedEstimate(estimates = ests[1:nOP],
      covMat = as.matrix(covMat[1:nOP,1:nOP]), invlink = exp,
      invlinkGrad = exp)

  detEstimates <- unMarkedEstimate(estimates = ests[(nOP + 1) : nP],
      covMat = as.matrix(covMat[(nOP + 1) : nP, (nOP + 1) : nP]), invlink = logistic,
      invlinkGrad = logistic.grad)

  umfit <- unMarkedFit(fitType = "occuRN",
      stateformula = stateformula, detformula = detformula,
      data = umf, stateEstimates = stateEstimates,
      detEstimates = detEstimates, AIC = fmAIC, hessian = fm$hessian)

#  umfit <- unMarkedFit(fitType = "occuRN",
#                       stateformula = stateformula, detformula = detformula,
#                       data = umf, stateMLE = ests[1:nOP],
#                       stateSE = ests.se[1:nOP],
#                       detMLE = ests[(nOP + 1) : nP],
#                       detSE = ests.se[(nOP + 1): nP], AIC = fmAIC)
  return(umfit)
}

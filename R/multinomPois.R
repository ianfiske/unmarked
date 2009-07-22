#' @include classes.R
#' @include utils.R
{}

## Take an M x J matrix of detection probabilities and return a matrix
## of M x J observation probs
#' Compute the cell probabilities for the observation classes
#' in removal sampling.
#'
#' Both p and the returned matrix are M x J for M sites and J sampling occasions.
#'
#' This function and \link{doublePiFun} provide example functions for computing multinomial
#' cell probabilities for \link{multinomPois}.
#'
#' @param p matrix of detection probabilities at each occasion for each site
#' @return matrix of cell probabilties for multinomial cells in removal sampling observation classes.
#' @export
removalPiFun <- function(p){
  M <- nrow(p)
  J <- ncol(p)
  pi <- matrix(NA, M, J)
  pi[,1] <- p[,1]
  for(i in seq(from = 2, length = J - 1)) {
	  pi[, i] <- pi[,i-1] / p[,i-1] * (1-p[,i-1]) * p[,i]
  }
  return(pi)
}

# p is an M x 2 matrix of detection probabilities (site x observer).
# returns an M x 3 matrix of row=(1 not 2, 2 not 1, 1 and 2).
#' Compute the cell probabilities for the observation classes
#' in double observer sampling.
#'
#' This function and \link{removalPiFun} provide example functions for computing multinomial
#' cell probabilities for \link{multinomPois}.
#'
#' @param p M x 2 matrix of detection probabilities at each occasion for each site
#' @return M x 3 matrix of cell probabilties for double observer observation classes.
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
#' appropriate function (piFun) argument.  \link{removalPiFun} and \link{doublePiFun}
#' are supplied as example cell probability functions.
#'
#' @title Multinomial-Poisson Mixtures
#' @param stateformula Right-hand side formula describing covariates of abundance
#' @param detformula Right-hand side formula describing covariates of detection
#' @param piFun Function to define multinomial cell probabilities.
#' @param umf unmarkedFrame supplying data.
#' @return unmarkedFit object describing the model fit.
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
#' ovenFrame <- unmarkedFrame(ovendata.list$data,
#'                            siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])))
#' fm1 <- multinomPois(~ 1 ~ ufp + trba, removalPiFun, ovenFrame)
#' fm2 <- multinomPois(~ 1 ~ ufp, removalPiFun, ovenFrame)
#' fm3 <- multinomPois(~ 1 ~ trba, removalPiFun, ovenFrame)
#' fm4 <- multinomPois(~ 1 ~ 1, removalPiFun, ovenFrame)
#' fm4
#' fm1
#' @export
multinomPois <-
function(formula, piFun, data) # TODO: remove piFun argument here and in examples.
{

	umf <- switch(class(data),
			data.frame = as(data, "unmarkedFrame"),
			unmarkedFrame = data,
			stop("Data is not a data frame or unmarkedFrame."))
	
	obsToY(umf) <- switch(as.character(substitute(piFun)),
			doublePiFun = matrix(c(1, 0, 0, 1, 1, 1), 2, 3),
			removalPiFun = {
				n.samples <- numY(umf)
				mat <- matrix(1, n.samples, n.samples)
				mat[col(mat) < row(mat)] <- 0
				mat
			},
			{
				if(is.null(obsToY(umf)) && missing(obsToY)) {
					stop("obsToY must be supplied in data or as an argument to multinomPois() if using a non-supplied piFun.")
				}
			})
	
	designMats <- getDesign2(formula, umf)
	X <- designMats$X; V <- designMats$V; y <- designMats$y; plotArea <- designMats$plotArea
  
  J <- ncol(y)
	R <- obsNum(umf)
  M <- nrow(y)

  lamParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nAP <- ncol(X)
  nP <- nDP + nAP
  yvec <- as.numeric(y)
  navec <- is.na(yvec)

  nll <- function(parms) {
    lambda <- exp(X %*% parms[1 : nAP]) * plotArea
    p <- plogis(V %*% parms[(nAP + 1) : nP])
    p.matrix <- matrix(p, M, R, byrow = TRUE)
    pi <- do.call(piFun, list(p = p.matrix))
    logLikeSite <- dpois(y, matrix(lambda, M, J) * pi, log = TRUE)
    logLikeSite[navec] <- 0
    -sum(logLikeSite)
  }

  fm <- optim(rep(0, nP), nll, method = "BFGS", hessian = TRUE)
	opt <- fm
  tryCatch(covMat <- solve(fm$hessian),
      error=simpleError("Hessian is not invertible.  Try using fewer covariates."))
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP
  names(ests) <- c(lamParms, detParms)

  stateEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
      estimates = ests[1:nAP],
      covMat = as.matrix(covMat[1:nAP,1:nAP]), invlink = "exp",
      invlinkGrad = "exp")

  detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
      estimates = ests[(nAP + 1) : nP],
      covMat = as.matrix(covMat[(nAP + 1) : nP, (nAP + 1) : nP]), invlink = "logistic",
      invlinkGrad = "logistic.grad")

  estimateList <- unmarkedEstimateList(list(state=stateEstimates,
          det=detEstimates))

  umfit <- unmarkedFit(fitType = "multinomPois",
      call = match.call(), formula = formula, data = umf, estimates = estimateList,
      AIC = fmAIC, opt = opt, negLogLike = fm$value, nllFun = nll)

  return(umfit)
}

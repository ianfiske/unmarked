#' @include unmarkedFit.R
#' @include unmarkedEstimate.R
#' @include utils.R
{}

#'  This function estimates the standard occupancy model of MacKenzie et al (2002).
#'
#'  See {\link{unmarkedFrame} for a description of how to supply data to the \code{umf}
#'  argument.
#'
#'  \code{occu} fits the standard occupancy model based on zero-inflated
#'  binomial models (MacKenzie et al. 2006, Royle and Dorazio
#'  2008).  The occupancy state process (\eqn{z_i}) of site \eqn{i} is
#'  modeled as
#'
#'  \deqn{z_i \sim Bernoulli(\psi_i)}{z_i ~ Bernoulli(psi_i)}
#'
#'  The observation process is modeled as
#'
#'  \deqn{y_{ij} \sim Bernoulli(p_{ij})}{y_ij ~ Bernoulli(p_ij)}
#'
#'  Covariates of \eqn{\psi_i}{psi_i} and \eqn{p_{ij}}{p_ij} are modeled
#'  using the logit link according to the \code{formula} argument.  The formula is a double right-hand sided formula
#' like \code{~ detform ~ occform} where \code{detform} is a formula for the detection process and \code{occform} is a 
#' formula for the partially observed occupancy state.  See \link{formula} for details on constructing model formulae
#'  in \R. 
#' @title Fit the MacKenzie Occupancy Model
#' @param formula double right-hand side formula describing covariates of detection and occupancy in that order.
#' @param data an unmarkedFrameOccu object (see \link{unmarkedFrame})..
#' @param knownOcc vector of sites that are known to be occupied.
#' @return unmarkedFitOccu object describing the model fit.
#' @references
#' MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege, J. Andrew Royle, and C. A. Langtimm. Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. Ecology 83, no. 8 (2002): 2248-2255. \cr
#' MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.  Amsterdam: Academic Press.  Royle, J. A. and R. Dorazio. (2008) \emph{Book Name}.
#' @author Ian Fiske
#' @examples
#' data(frogs)
#' pferUMF <- unmarkedFrameOccu(pfer.bin)
#' 
#' # add some fake covariates for illustration
#' siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))
#' 
#' # observation covariates are in site-major, observation-minor order
#' obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))
#' 
#' fm <- occu(~ obsvar1 ~ 1, pferUMF)
#' fm
#' 
#' confint(fm, type='det', method = 'normal')
#' confint(fm, type='det', method = 'profile')
#' 
#' # estimate detection effect at obsvars=0.5
#' lc <- linearComb(fm['det'],c(1,0.5))
#' lc
#' 
#' # transform this to probability (0 to 1) scale and get confidence limits
#' btlc <- backTransform(lc)
#' btlc
#' confint(btlc)
#' @keywords models
#' @export
occu <-
function(formula, data, knownOcc = numeric(0), starts)
{
	if(!is(data, "unmarkedFrameOccu")) stop("Data is not an unmarkedFrameOccu object.")
		
  designMats <- getDesign2(formula, data)
	X <- designMats$X; V <- designMats$V; y <- designMats$y
	
  y <- truncateToBinary(y)
  J <- ncol(y)
  M <- nrow(y)

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

	if(missing(starts)) starts <- rnorm(nP)
  fm <- optim(starts, nll, method = "BFGS", hessian = TRUE)
	opt <- fm
  tryCatch(covMat <- solve(fm$hessian),
      error=function(x) simpleError("Hessian is not invertible.  Try using fewer covariates."))
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP + 2*nP*(nP + 1)/(M - nP - 1)
  names(ests) <- c(occParms, detParms)

  state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
      estimates = ests[1:nOP],
      covMat = as.matrix(covMat[1:nOP,1:nOP]), invlink = "logistic",
      invlinkGrad = "logistic.grad")

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
      estimates = ests[(nOP + 1) : nP],
      covMat = as.matrix(covMat[(nOP + 1) : nP, (nOP + 1) : nP]), invlink = "logistic",
      invlinkGrad = "logistic.grad")

  estimateList <- unmarkedEstimateList(list(state=state, det=det))

  umfit <- new("unmarkedFitOccu", fitType = "occu",
      call = match.call(), formula = formula, data = data, sitesRemoved = designMats$removed.sites, 
			estimates = estimateList,
      AIC = fmAIC, opt = opt, negLogLike = fm$value, nllFun = nll)

  return(umfit)
}

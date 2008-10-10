#' @include classes.R
roxygen()

#' Fit the N-mixture point count model 
#'
#' This function fits the standard occupancy model of 
#'
#' @param stateformula Right-hand side formula describing covariates of abundance
#' @param detformula Right-hand side formula describing covariates of detection
#' @param umf an unMarkedFrame supplying data to the model.
#' @return unMarkedFit object describing the model fit.
#' @references (Royle 2004)
#' @examples
#' data(mallard)
#' mallardUMF <- unMarkedFrame(mallard.y, siteCovs = mallard.site,
#'                            obsCovs = mallard.obs)
#' fm.nmx1 <- pcount(~ length + elev + forest, ~ ivel+ date + I(date^2),
#'                   mallardUMF)
#' @export
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
  k.jik <- rep(k, M*J)

  nP <- nAP + nDP + ifelse(identical(mixture,"NB"),1,0)
  y.ji <- as.numeric(y)
  y.jik <- rep(y.ji, each = K + 1)
  navec <- is.na(y.jik)
  nd <- ifelse(rowSums(y, na.rm=TRUE) == 0, 1, 0) # I(no detection at site i)

  nll <- function(parms){

    theta.i <- exp(X %*% parms[1 : nAP])
    p.ji <- plogis(V %*% parms[(nAP + 1) : (nAP + nDP)])
    theta.ik <- rep(theta.i, each = K + 1)
    p.jik <- rep(p.ji, each = K + 1)

    bin.jik <- dbinom(y.jik,k.jik,p.jik)
    bin.jik[which(is.na(bin.jik))] <- 1
    bin.ik.mat <- matrix(bin.jik, M * (K + 1), J)
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
  ests.se <- diag(solve(fm$hessian))
  fmAIC <- 2 * fm$value + 2 * nP
  names(ests.se) <- c(lamParms, detParms)
  umfit <- unMarkedFit(fitType = "pcount",
                       stateformula = stateformula, detformula = detformula,
                       data = umf, stateMLE = ests[1:nAP],
                       stateSE = ests.se[1:nAP], 
                       detMLE = ests[(nAP + 1) : nP],
                       detSE = ests.se[(nAP + 1): nP], AIC = fmAIC)
  return(umfit)
}

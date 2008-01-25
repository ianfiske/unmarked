# model allows covariates at the observation level stateFormula and
# detFormula are as above obsdata can either be an array of M x J
# matrices or a list of M x J matrices sitedata is optional, as all
# site vectors can be in the obsdata list

occu <-
function(stateformula, detformula,
         data = list(y = y, covdata.site = covdata.site, covdata.obs = covdata.obs),
         y, covdata.site = NULL, covdata.obs = NULL)
{

  arranged <- arrangeData(data)
  cleaned <- handleNA(arranged, stateformula, detformula)
  y <- cleaned$y
  sitedata <- cleaned$covdata.site
  obsdata <- cleaned$covdata.obs

  design <- getDesign(stateformula = stateformula, detformula = detformula,
    y = y, sitedata = sitedata, obsdata = obsdata)  
  nOP <- design$nOP
  nDP <- design$nDP
  XDet <-design$XDet
  XOcc <- design$XOcc
  occParms <- design$occParms
  detParms <- design$detParms

  J <- ncol(y)
  M <- nrow(y)

  nP <- nDP + nOP
  yvec <- as.numeric(y)
  navec <- is.na(yvec)
  nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i indicator
  
  nll <- function(parms) {

    psi <- plogis(XOcc %*% parms[1:nOP])
      pvec <- plogis(XDet %*% parms[(nOP+1):nP])
      cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
      cp[navec] <- 1  # so that NA's don't modify likelihood        
      cpmat <- matrix(cp,M,J) # put back into matrix to multiply appropriately
      loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi)) 
      sum(-1 * loglik)
  }
  
  fm <- optim(rep(0, nP), nll, method = "BFGS")
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP
  names(ests) <- c(occParms, detParms)
  list(estimates = ests, AIC = fmAIC)
}



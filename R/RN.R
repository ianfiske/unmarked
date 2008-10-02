# fit the RN model of Royle and Nichols (2003)

occuRN <- 
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
  occParms[1] <- "lamconst"
  detParms[1] <- "rconst"
  
  K <- 20
  M <- nrow(y)
  J <- ncol(y)

  nP <- nDP + nOP
  y.ji <- as.numeric(y)
  y.jik <- rep(y.ji, each = K + 1)
  navec <- is.na(y.jik)
  nd <- ifelse(rowSums(y, na.rm=TRUE) == 0, 1, 0) # no det site indicator
  k <- 0:K
  k.i <- rep(k, M)
  k.ji <- rep(k, M * J)
  
  nll <- function(parms, f = "Poisson")
  {    
    lam.i <- exp(XOcc %*% parms[1 : nOP])
    lam.ik <- rep(lam.i, each = K + 1)
    r.ji <- plogis(XDet %*% parms[(nOP + 1) : nP])

    r.jik <- rep(r.ji, each = K + 1)
    p.sup <- 1 - (1 - r.jik)^(k.ji)
    cp <- p.sup^y.jik * (1 - p.sup)^(1 - y.jik)
    cp[navec] <- 1
    cp.mat <- matrix(cp, M * (K + 1), J)
    p.ik <- rowProds(cp.mat)

    f.ik <- dpois(k.i, lam.ik)
    dens.mat <- matrix(p.ik * f.ik, M, K + 1, byrow = TRUE)
    dens.integ <- rowSums(dens.mat)
  
    nll <- - sum(log(dens.integ))      # multiply likelihood over all sites
    nll
  }

  fm <- optim(rep(0,nP), nll, method = "BFGS")
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP + 2 * nP * (nP + 1) / (M - nP - 1)
  names(ests) <- c(occParms,detParms)
  list(estimates = ests, AIC_c = fmAIC)
}

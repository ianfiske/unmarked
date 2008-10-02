# fit the N-mixture point count model (Royle 2004)
pcount <-
function(stateformula, detformula,
         data = list(y = y, covdata.site = covdata.site, covdata.obs = covdata.obs),
         y, covdata.site = NULL, covdata.obs = NULL, K = NULL, mixture = "P")
{ 
  if ((mixture %in% c("P","NB")) == FALSE) stop("Mixture familiy not recognized. Please choose \"P\" or \"NB\".")
  
  arranged <- arrangeData(data)
  cleaned <- handleNA(arranged, stateformula, detformula)
  y <- cleaned$y
  sitedata <- cleaned$covdata.site
  obsdata <- cleaned$covdata.obs

  design <- getDesign(stateformula = stateformula, detformula = detformula,
    y = y, sitedata = sitedata, obsdata = obsdata)
  nNP <- design$nOP
  nDP <- design$nDP
  XDet <-design$XDet
  XN <- design$XOcc
  NParms <- design$occParms
  detParms <- design$detParms
  NParms[1] <- "lamconst"
  detParms[1] <- "pconst"

  if(is.null(K)) K <- max(y, na.rm = TRUE) + 20
  if(K <= max(y, na.rm = TRUE))
    stop("specified K is too small. Try a value larger than any observation") 
  k <- 0:K
  M <- nrow(y)
  J <- ncol(y)
  k.ik <- rep(k, M)
  k.jik <- rep(k, M*J)

  nP <- nNP + nDP + ifelse(identical(mixture,"NB"),1,0)
  y.ji <- as.numeric(y)
  y.jik <- rep(y.ji, each = K + 1)
  navec <- is.na(y.jik)
  nd <- ifelse(rowSums(y, na.rm=TRUE) == 0, 1, 0) # I(no detection at site i)

  nll <- function(parms){

    theta.i <- exp(XN %*% parms[1 : nNP])
    p.ji <- plogis(XDet %*% parms[(nNP + 1) : (nNP + nDP)])
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
  
  if(identical(mixture,"P")) {
    fm <- optim(rep(0,nP), nll, method="BFGS")
  }
  else if (identical(mixture,"NB")){
    fm <- optim(c(rep(0, nNP + nDP),1), nll, method="BFGS")
  }
  ests <- fm$par
  fmAIC <- 2*fm$value + 2*nP

  if(identical(mixture,"NB"))
     {nbParm <- "alpha"}
  else
     {nbParm <- character(0)}
  names(ests) <- c(NParms, detParms, nbParm)
  list(estimates = ests, AIC = fmAIC)
}

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

mnMix <-
function(stateformula = ~ 1, detformula = ~ 1,
         data = list(y = y, covdata.site = covdata.site, covdata.obs = covdata.obs),
         y, covdata.site = NULL, covdata.obs = NULL, constraint = NULL)
{

  arranged <- arrangeData(data)
  cleaned <- handleNA(arranged, stateformula, detformula)
  y <- cleaned$y
  sitedata <- cleaned$covdata.site
  obsdata <- cleaned$covdata.obs
  
  design <- getDesign(stateformula = stateformula, detformula = detformula,
    y = y, sitedata = sitedata, obsdata = obsdata)
  nSP <- design$nOP
  nDCP <- design$nDP
  XDet <-design$XDet
  XN <- design$XOcc
  NParms <- design$occParms
  detParms <- design$detParms
  NParms[1] <- "lamconst"
  detParms[1] <- "pconst"

  M <- nrow(y)
  J <- ncol(y)

  K <- max(y, na.rm = TRUE)
  con <- constraint
  nDMP.un <- K*(K+1)/2
  if (is.null(con)) con = c(1:K, rep(K + 1, nDMP.un - K))

  # create design matrix with alpha_k intercept for each k=1,2,...,K
  XDet <- XDet[, -1] # remove intercept
  nDCP <- nDCP - 1
  XDet.jik <- XDet %x% rep(1, K)  # repeat rows of X, each = K
  k.diag <- rep(1, M * J) %x% diag(K) # add intercepts the alpha_k intercepts
  XDet.jik <- cbind(k.diag, XDet.jik)

  nDMP <- max(con)    
  nDP <- nDCP + nDMP
  nSP <- K
  nP <- nDP + nSP

  # construct constraint equation
  H <- matrix(0,nDMP.un,nDMP)
  for(i in 1:nDMP.un){
    H[i,con[i]] <- 1
  }
 
  # compute indices
  det.row <- c(rep(NA,K), rep(3:(K+1), times = 1 : (K - 1))) 
  det.col <- c(rep(NA,K), sequence(1 : (K - 1)) + 1)
  arr.offset <- 0:(M * J - 1) * (K + 1)^2
  diag.els <- K*(1:K) + (2:(K+1))*2 - 1
  lower.els <- which((lower.tri(matrix(1,K+1,K+1))),arr.ind=T)
  lower.els <- lower.els[order(lower.els[,1]),] # reorder indices by row
  lower.els <- (lower.els[,2]-1)*(K+1) + lower.els[,1]  # compute vector indices
  diag.els.arr <- rep(arr.offset, each = K) + diag.els  # get offsetted indices
  lower.els.arr <- rep(arr.offset, each = sum(1:K)) + lower.els
  
  # vectorized version of detMatrix
  detMatrix <- function(dPars) {
    nmats <- nrow(dPars)#nrow(dPars)
    detMat <- lower.tri(matrix(1,(K +1),(K+1)),diag = TRUE) %x%
      array(1,c(1,1,nmats))
    # put the p's in the detMats
    detMat[diag.els.arr] <- t(dPars[, 1:K])
    detMat[lower.els.arr] <- 1 - t(dPars[,rep(1:K,times=1:K)])
    # put beta's in the mats
    for(i in (K+1):nDMP.un){
      detMat[det.row[i], det.col[i],] <- dPars[,i] *
        detMat[det.row[i], det.col[i],]
      detMat[det.row[i], 1:(det.col[i] - 1),] <- (1 - dPars[,i]) * 
        detMat[det.row[i], 1:(det.col[i] - 1),]
    }
    detMat
  }

  y.ji <- as.numeric(y)  
  y.jik <- rep(y.ji, each = K + 1)
  K.jik <- rep(0:K, M*J)
  
  nll <- function(parms) {
    # recover full parameters
    dPars <- H %*% parms[1:nDMP]

    # recover parameters
    alpha <- dPars[1 : K]
    beta <- plogis(dPars[(K + 1) : nDMP.un])
    b <- if(nDCP > 0) {parms[(nDMP + 1) : nDP]}
         else {NULL}
    psi <- parms[(nDP + 1) : nP]
    psi <- exp(c(0,psi))/sum(exp(c(0,psi)))
    
    beta.ji.mat <- matrix(beta, M*J, length(beta), byrow=TRUE)
 
    # model parms 
    p.jik <- plogis(XDet.jik %*% c(alpha, b))

    # get detMat paramters (rows are alphas then betas)
    p.ji.k.mat <- matrix(p.jik, M*J, K, byrow=TRUE)
    detMat.pars <- cbind(p.ji.k.mat, beta.ji.mat)

    detMats.ji <- detMatrix(detMat.pars)
    detMats.ji <- matrix(detMats.ji, M * J, (K + 1)^2, byrow = T)
    detMats.jik <- detMats.ji %x% rep(1, K + 1)
    
    rcumsums <- (K+1)^2*(0:(nrow(detMats.jik)-1))
    detMats.jik <- as.numeric(t(detMats.jik))
    f.jik <- detMats.jik[(y.jik * (K + 1) + K.jik + 1) + rcumsums]

    f.ik <- matrix(f.jik, M*(K+1), J)
    f.ik[is.na(f.ik)] <- 1
    
    psi.ik <- rep(psi,M)
    g.ik <- rowProds(f.ik)*psi.ik
    g.i <- rowSums(matrix(g.ik, M, K + 1, byrow = TRUE))
    -sum(log(g.i))
  }
  
  fm <- optim(rep(0,nP),nll, method = "BFGS")
  ests <- fm$par

  names(ests) <- c(paste(rep("detmat",nDMP),1:nDMP, sep=""), 
    eval(if(nDCP > 0) paste(rep("b", nDCP), 1:nDCP, sep="") else NULL),
    paste(rep("psi",nSP), 1:nSP, sep=""))
  list(estimates = ests, AIC = 2*fm$value + 2*nP)
}

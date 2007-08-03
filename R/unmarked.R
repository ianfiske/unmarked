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
  cleaned <- handleNA(arranged)
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


# fit the RN model of Royle and Nichols (2003)

occuRN <- 
function(stateformula, detformula,
         data = list(y = y, covdata.site = covdata.site, covdata.obs = covdata.obs),
         y, covdata.site = NULL, covdata.obs = NULL)
{

  arranged <- arrangeData(data)
  cleaned <- handleNA(arranged)
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
  cleaned <- handleNA(arranged)
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
  cleaned <- handleNA(arranged)
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

markovMN <-
function(stateformula = ~ 1, detformula = ~ 1,
         data = list(y = y, covdata.site = covdata.site, covdata.obs = covdata.obs),
         y, covdata.site = NULL, covdata.obs = NULL, detconstraint = NULL,
         phiconstraint = NULL, J)
{
  arranged <- arrangeData(data)
  cleaned <- handleNA(arranged)
  y <- cleaned$y
  sitedata <- cleaned$covdata.site
  obsdata <- cleaned$covdata.obs

  M <- nrow(y)
  nY <- ncol(y)/J
  K <- max(y, na.rm = TRUE)
  
  design <- getDesign(stateformula = stateformula, detformula = detformula,
                      y = y, sitedata = sitedata, obsdata = obsdata)
  nSP <- design$nOP
  nDCP <- design$nDP
  XN <- design$XOcc
  NParms <- design$occParms
  detParms <- design$detParms
  NParms[1] <- "lamconst"
  detParms[1] <- "pconst"
  XDet.tji <- design$XDet
  nDCP <- nDCP - 1
  
  XDet.tji <- XDet.tji[, -1] # remove intercept
  XDet.tjik <- XDet.tji %x% rep(1, K)  # repeat rows of X, each = K
  k.diag <- rep(1, M * J * nY) %x% diag(K) # add intercepts the alpha_k
  yr.int <- diag(nY) %x% rep(1, M * K * J)  # add intercepts for gma_k
  XDet.tjik <- cbind(yr.int, k.diag, XDet.tjik)
  
  con <- detconstraint
  nDMP.un <- K*(K+1)/2
  if (is.null(con)) con <- c(1:K, rep(K + 1, nDMP.un - K))

  phicon <- phiconstraint
  # default phi constraint is different s's, equal growth, and equal reg
  if(is.null(phicon)) phicon <- c(1:(K+1), rep(K + 2, K),
                                  rep(K + 3, K + 2))
  nPhiP.un <- K * (K + 1)
  
  nDMP <- max(con)
  nDP <- nDCP + nDMP + nY
  nSP <- K
  nPhiP <- max(phicon)
  nP <- nDP + nSP + nPhiP

  # construct constraint matrix.  note that here we constrain beta's and
  # alphas, whereas the paper talks about constraining p's and beta's
  H <- matrix(0,nDMP.un,nDMP)
  for(i in 1:nDMP.un){
    H[i,con[i]] <- 1
  }

  # construct constrain matrix for phi paramters
  H.phi <- matrix(0, nPhiP.un, nPhiP)
  for(i in 1:nPhiP.un){
    H.phi[i, phicon[i]] <- 1
  }

  # compute indices for vectorization in detMatrix:
  # det.row and .col give indices for betas in detection matrix
  det.row <- c(rep(NA,K), rep(3:(K+1), times = 1 : (K - 1))) 
  det.col <- c(rep(NA,K), sequence(1 : (K - 1)) + 1)
  # compute offset for each matrix
  arr.offset <- 0:(M * J * nY - 1) * (K + 1)^2
  # compute indices of diagonal elements
  diag.els <- K*(1:K) + (2:(K+1))*2 - 1
  # compute indices of lower triangular indices and
  # sort them by row because that is how the (1-p_i) terms are repeated
  # and translate them into vector indices
  lower.els <- which((lower.tri(matrix(1,K+1,K+1))),arr.ind=T)
  lower.els <- lower.els[order(lower.els[,1]),]
  lower.els <- (lower.els[,2]-1)*(K+1) + lower.els[,1]
  # compute indices for entire set of matrices with offsets
  diag.els.arr <- rep(arr.offset, each = K) + diag.els
  lower.els.arr <- rep(arr.offset, each = sum(1:K)) + lower.els
  
  # creates detection matrix from list of detection parameters
  detMatrix <- function(dPars) {
    nmats <- nrow(dPars)#nrow(dPars)
    detMat <- lower.tri(matrix(1, K + 1, K + 1),diag = TRUE) %x%
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

  # create process matrix from vector of paramters
  # WARNING: this is currently frog (K=3)-specific
  # needs to be modified for general or other use
  phiMatrix <- function(phiPars) {
    phi <- matrix(NA, K + 1, K + 1)
    s <- phiPars[1 : (K + 1)]
    g <- phiPars[(K + 2) : (K + 4)]
    r <- phiPars[(K + 5) : 12]
    diag(phi) <- s
    phi[1,2:4] <- (1-s[1])
    phi[2,c(1,3,4)] <- (1-s[2])
    phi[3,c(1,2,4)] <- (1-s[3])
    phi[4,1:3] <- (1-s[4])

    # r = r_10, r_21, r_20, r_32, r_31
    # g = g_01, g_02, g_12
    phi[1,2:4] <- phi[1,2:4] * c(g[1], 1-g[1], 1-g[1]) *
      c(1, g[2], 1 - g[2])
    phi[2,c(1,3,4)] <- phi[2,c(1,3,4)] * c(1-g[3], g[3], 1-g[3]) *
      c(r[1], 1, 1 - r[1])
    phi[3,c(1,2,4)] <- phi[3,c(1,2,4)] * c(1-r[2], r[2], 1 - r[2]) *
      c(r[3], 1, 1 - r[3])
    phi[4,1:3] <- phi[4,1:3] * c(1-r[4], 1-r[4], r[4]) *
      c(1-r[5], r[5], 1)
    
    phi    
  }
  
  y.tji <- as.numeric(y)  
  y.tjik <- rep(y.tji, each = K + 1)
  K.tjik <- rep(0:K, M * J * nY)
  t.tjik <- rep(1:nY, each = (K+1) *M*J)
  i.tjik <- rep(rep(1:M, each = K + 1), nY*J)
  j.tjik <- rep(rep(1:J, each = (K + 1)*M), nY)
  total <- nY*J*M*(K+1)
  ind.k.y.tjik <- matrix(c(K.tjik + 1, y.tjik + 1, 1:total),total, 3)

  # parms: alpha's, beta's, b's, gamma's, psi's, phi's
  iteration <- 1
  nll <- function(parms) {
    # recover full parameters: alpha's and beta's
    dPars <- H %*% parms[1:nDMP]

    # recover all parameters
    alpha <- dPars[1 : K]
    beta <- plogis(dPars[(K + 1) : nDMP.un])
    b <- if(nDCP > 0) {parms[(nDMP + 1) : (nDMP + nDCP)]}
    gma <- parms[(nDMP + nDCP + 1) : nDP]

    psi <- parms[(nDP + 1) : (nDP + K)]
    psi <- exp(c(0,psi))/sum(exp(c(0,psi)))

    phiPars <- H.phi %*% plogis(parms[(nDP + K  + 1) : nP])
    phi <- phiMatrix(phiPars)

    
    beta.tji <- matrix(beta, nY * M * J, length(beta), byrow=TRUE)

    # model detection parms (gammas, alphas, and bs)
    p.tjik <- plogis(XDet.tjik %*% c(gma, alpha, b))
#    p.tjik[is.na(p.tjik)] <- 1 # CONSIDER THIS STEP FUTHER!!!!!

    # Get detmat Paramters (Rows Are Alphas Then Betas)
    p.tji.k <- matrix(p.tjik, nY * M * J, K, byrow=TRUE)
    detMat.pars <- cbind(p.tji.k, beta.tji)

    detMats.tji <- detMatrix(detMat.pars)
    detMats.tjik <- detMats.tji %x% array(1,c(1,1,K+1))
    fy.tjik <- detMats.tjik[ind.k.y.tjik]
    fy.tjik[is.na(y.tjik)] <- 1  # RECONSIDER THIS STEP!!!!

    fy.ik.j.t <- array(fy.tjik, c(M * (K + 1), J, nY))
    fy.ik.t.j <- aperm(fy.ik.j.t, c(1,3,2))
    fy.tik.j <- matrix(fy.ik.t.j, nY * M * (K + 1))
    fy.tik <- rowProds(fy.tik.j, na.rm = TRUE)

    # compute products D(p(y_it)) * phi_t for t = 1,..., T-1
    fy.k.1.ti <- array(fy.tik[1:(M * (K + 1) * (nY - 1))],
                       c(K + 1, 1, (nY - 1) * M))
    fy.k.k.ti <- fy.k.1.ti %x% t(rep(1, K + 1)) # repeat columns
    phi.ti <- array(1,c(1,1, M * (nY - 1))) %x% phi
    phi.ti.prod <- fy.k.k.ti * phi.ti

    # try brute force:
    psi.t <- array(psi, c(1, K + 1, M))
    for(t in 1 : (nY - 1)){
      for(i in 1 : M) {
        psi.t[,,i] <- psi.t[,,i] %*% phi.ti.prod[,,(t-1)*M + i]
      }
    }

    fy.k.1.Ti <- array(fy.tik[(M * (K + 1) * (nY - 1) + 1) :
                              (M * (K + 1) * nY)],
                       c(K + 1, 1, M))
    l.i <- numeric(M)
    for(i in 1:M){
      l.i[i] <- psi.t[,,i] %*% fy.k.1.Ti[,,i]
    }

    nLL <- -sum(log(l.i))
    print(sprintf("%i: %f",get("iteration",parent.frame(3)), nLL))
    eval.parent(quote(iteration <- iteration + 1), 3)
    nLL
  }

  fm <- optim(rep(0,nP),nll, method = "BFGS", hessian = TRUE,
              control = list(trace = 6))
  ests <- fm$par

  DMP <- H %*% ests[1:nDMP]
  psi <- ests[(nDP + 1) : (nDP + K)]
  psi <- exp(c(0,psi))/sum(exp(c(0,psi)))

  phiPars <- H.phi %*% plogis(ests[(nDP + K  + 1) : nP])
  phi <- phiMatrix(phiPars)

  list(alpha = DMP[1:K], beta = plogis(DMP[(K+1): nDMP.un]),
       b = if(nDCP > 0) {ests[(nDMP + 1) : (nDMP + nDCP)]},
       gamma = ests[(nDMP + nDCP + 1) : nDP], psi = psi,
       phiPars = as.numeric(phiPars), phi = phi, hessian = fm$hessian,
       AIC = 2*fm$value + 2*nP)
}

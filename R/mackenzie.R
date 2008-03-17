mackenzie <-
function(stateformula = ~ 1, detformula = ~ 1,
         data = list(y = y, covdata.site = covdata.site,
         covdata.obs = covdata.obs), y, covdata.site = NULL,
         covdata.obs = NULL, J,  yearly_alpha = FALSE)
{

  arranged <- arrangeData(data)
  cleaned <- handleNA(arranged, stateformula, detformula)
  y <- cleaned$y
  sitedata <- cleaned$covdata.site
  obsdata <- cleaned$covdata.obs

  M <- nrow(y)
  nY <- ncol(y)/J
  K <- 1 # max of states

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
  nDCP <- nDCP - 1 # number of detection covariates

  XDet.tji <- XDet.tji[, -1] # remove intercept
#  XDet.tjik <- XDet.tji %x% rep(1, K)  # repeat rows of X, each = K
  if(yearly_alpha) {
      yr.int <- diag(nY) %x% rep(1, M * J)  # add intercepts for gma_k
  } else {
      yr.int <- rep(1, M*J*nY)
  }
  XDet.tji <- cbind(yr.int, XDet.tji)
  

   nDMP <- 0  # number of independent detection matrix parameters among
##                                         # the p's and beta's
   nDP <- nDCP + nDMP + ifelse(yearly_alpha,nY,1)  # number of detection parameters including
                                        # covariates (b's), alpha's
   nSP <- K                # number of parameters for psi vector of initial
                                         # latent abundances
   nPhiP <- 2    # number of independent transition matrix pars
   nP <- nDP + nSP + nPhiP  # total number of parameters


# note that now there is only one dPar (p)  
p.inds <- matrix(c(rep(2,M*J*nY),rep(2,M*J*nY),1:(M*J*nY)),M*J*nY,3)
np.inds <- matrix(c(rep(2,M*J*nY),rep(1,M*J*nY),1:(M*J*nY)),M*J*nY,3)
                
detMatrix <- function(dPars) {
   nmats <- nrow(dPars)
   detMat <- array(1,c(1,1,nmats)) %x% matrix(c(1,NA,0,NA),2,2)
   detMat[p.inds] <- dPars
   detMat[np.inds] <- 1-dPars
   detMat
}


# phiPars is (epsilon, gamma)
phiMatrix <- function(phiPars) {
    epsilon <- phiPars[1]
    gamma <- phiPars[2]
    phi <- matrix(c(1-gamma,epsilon,gamma,1-epsilon), 2, 2)
    phi
}
  
  y.tji <- as.numeric(y)  
  y.tjik <- rep(y.tji, each = K + 1)
  K.tjik <- rep(0:K, M * J * nY)
  t.tjik <- rep(1:nY, each = (K+1) *M*J)
  i.tjik <- rep(rep(1:M, each = K + 1), nY*J)
  j.tjik <- rep(rep(1:J, each = (K + 1)*M), nY)
  total <- nY*J*M*(K+1)
  # calc indices of columns of detection matrices
  ind.k.y.tjik <- matrix(c(K.tjik + 1, y.tjik + 1, 1:total),total, 3)

  # parms: alpha, b's, gamma's, psi's, phi's
  iteration <- 1
  nll <- function(parms) {

    # recover detection parameters
    # see Royle and Link 2005 for further explanations of alpha's,
    # beta's, and b's.  Here is the general equation for the stasis
    # (diagonal) detection term in detection matrix:
    # p_tjik = logistic(gma_t + alpha_k + sum_l^L(b_l*x_l))
    #alpha <- dPars[1 : K]    # "stasis" detection intercepts for each state
    #beta <- plogis(dPars[(K + 1) : nDMP]) # underdetection parameters
    if(nDCP > 0) { # detection parameters for covariates
        b <- parms[(nDMP + 1) : (nDMP + nDCP)]
    }

    # recover alpha's
    alpha <- parms[(nDMP + nDCP + 1) : nDP]

    # recover the initial latent abundance vector
    psi <- parms[(nDP + 1) : (nDP + K)]
    psi <- exp(c(0,psi))/sum(exp(c(0,psi)))

    # get transition matrix
    phiPars <- plogis(parms[(nDP + K  + 1) : nP])
    phi <- phiMatrix(phiPars)

    # model detection parms (alphas, and bs)
    p.tji <- plogis(XDet.tji %*% c(alpha, b))
    p.tji[is.na(p.tji)] <- 1 # CONSIDER THIS STEP FUTHER!!!!!

##     # Get detmat Paramters (Rows Are p's Then Betas)
##     detMat.pars <- matrix(p.tjik, nY * M * J, K, byrow=TRUE)
    detMat.pars <- p.tji

    detMats.tji <- detMatrix(detMat.pars)
    detMats.tjik <- detMats.tji %x% array(1,c(1,1,K+1))
    fy.tjik <- detMats.tjik[ind.k.y.tjik] # pull columns of detection matrices
    fy.tjik[is.na(y.tjik)] <- 1  # RECONSIDER THIS STEP!!!!

    fy.ik.j.t <- array(fy.tjik, c(M * (K + 1), J, nY))
    fy.ik.t.j <- aperm(fy.ik.j.t, c(1,3,2))
    fy.tik.j <- matrix(fy.ik.t.j, nY * M * (K + 1), J)
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

  fm <- optim(rep(0,nP),nll, method = "BFGS", hessian = TRUE)
  ests <- fm$par

  DMP <- ests[1:nDMP]
  psi <- ests[(nDP + 1) : (nDP + K)]
  psi <- exp(c(0,psi))/sum(exp(c(0,psi)))

  phiPars <- plogis(ests[(nDP + K  + 1) : nP])
  phi <- phiMatrix(phiPars)

  list(b = if(nDCP > 0) {ests[(nDMP + 1) : (nDMP + nDCP)]},
       alpha = ests[(nDMP + nDCP + 1) : nDP], psi = psi,
       phiPars = as.numeric(phiPars), phi = phi, hessian = fm$hessian,
       AIC = 2*fm$value + 2*nP)
}

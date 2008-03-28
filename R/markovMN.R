markovMN <-
function(stateformula = ~ 1, detformula = ~ 1,
         data = list(y = y, covdata.site = covdata.site, covdata.obs = covdata.obs),
         y, covdata.site = NULL, covdata.obs = NULL, detconstraint = NULL,
         phiconstraint = NULL, J)
{

  arranged <- arrangeData(data)
  cleaned <- handleNA(arranged, stateformula, detformula)
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

  XDet.tji <- as.matrix(XDet.tji[, -1]) # remove intercept
###  XDet.tjik <- XDet.tji %x% rep(1, K)  # repeat rows of X, each = K
  XDet.tjik <- matrix(rep(XDet.tji,each=K),M*J*nY*K,ncol(XDet.tji)) # test version
  k.diag <- rep(1, M * J * nY) %x% diag(K) # add intercepts the alpha_k

  ## add intercepts for gma_k
  ## to keep X full column rank, use first year as baseline
  yr.int <- diag(nY - 1)
  yr.int <- rbind(rep(0, nY - 1), yr.int)
  yr.int <- yr.int %x% rep(1, M * K * J)
  
  XDet.tjik <- cbind(yr.int, k.diag, XDet.tjik)

  con <- detconstraint
  nDMP.un <- K*(K+1)/2
  if (is.null(con)) con <- c(1:K, rep(K + 1, nDMP.un - K))

  phicon <- phiconstraint
  # default phi constraint is different s's, equal growth, and equal reg
  if(is.null(phicon)) phicon <- c(1:(K+1), rep(K + 2, K),
                                  rep(K + 3, K + 2))
  nPhiP.un <- K * (K + 1)  # number of transition matrix parameters in
                                        # an unconstrained matrix
  nDMP <- max(con)  # number of independent detection matrix parameters among
                                        # the p's and beta's
  nDP <- nDCP + nDMP + nY - 1  # number of detection parameters including
                                        # covariates (b's), alpha's, beta's,
                                        # year effects (gamma's)
  nSP <- K                # number of parameters for psi vector of initial
                                        # latent abundances
  nPhiP <- max(phicon)    # number of independent transition matrix pars
  nP <- nDP + nSP + nPhiP  # total number of parameters

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
  nmats <- M*J*nY
  detMat <- lower.tri(matrix(1, K + 1, K + 1),diag = TRUE) %x%
      array(1,c(1,1,nmats))

  # creates detection matrix from list of detection parameters
  detMatrix <- function(dPars) {
    # put the p's in the detMats
    detMat[diag.els.arr] <- t(dPars[, 1:K])
    detMat[lower.els.arr] <- 1 - t(dPars[,rep(1:K,times=1:K)])
    # put beta's in the mats
    for(i in (K+1):nDMP.un) {
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
  #ind.k.y.tjik <- matrix(c(K.tjik + 1, y.tjik + 1, 1:total),total, 3)
  ind.k.y.tjik2 <- matrix(c(K.tjik + 1, y.tjik + 1,
                            rep(1:total,each = K + 1)), total, 3)

  # parms: alpha's, beta's, b's, gamma's, psi's, phi's
  iteration <- 1
  nll <- function(parms) {
    # compute detection parameters, alpha's and beta's, from constrained
    # parameters
    dPars <- H %*% parms[1:nDMP]

    # recover detection parameters
    # see Royle and Link 2005 for further explanations of alpha's,
    # beta's, and b's.  Here is the general equation for the stasis
    # (diagonal) detection term in detection matrix:
    # p_tjik = logistic(gma_t + alpha_k + sum_l^L(b_l*x_l))
    alpha <- dPars[1 : K]    # "stasis" detection intercepts for each state
    beta <- plogis(dPars[(K + 1) : nDMP.un]) # underdetection parameters
    if(nDCP > 0) {
        b <- parms[(nDMP + 1) : (nDMP + nDCP)]
    } # "stasis"

    ## detection parameters for covariates
    gma <- parms[(nDMP + nDCP + 1) : nDP] # stasis detection effect for year

    # recover the initial latent abundance vector
    psi <- parms[(nDP + 1) : (nDP + K)]
    psi <- exp(c(0,psi))/sum(exp(c(0,psi)))

    # get transition matrix
    phiPars <- H.phi %*% plogis(parms[(nDP + K  + 1) : nP])
    phi <- phiMatrix(phiPars)

    # create matrix of repeated covariate parameters for vectorization
    beta.tji <- matrix(beta, nY * M * J, length(beta), byrow=TRUE)

    # model detection parms (gammas, alphas, and bs)
    p.tjik <- plogis(XDet.tjik %*% c(gma, alpha, b))
#    p.tjik[is.na(p.tjik)] <- 1 # CONSIDER THIS STEP FUTHER!!!!!

    # Get detmat Paramters (Rows Are Alphas Then Betas)
    p.tji.k <- matrix(p.tjik, nY * M * J, K, byrow=TRUE)
    detMat.pars <- cbind(p.tji.k, beta.tji)

    detMats.tji <- detMatrix(detMat.pars)
    #detMats.tjik <- detMats.tji %x% array(1,c(1,1,K+1))
    #fy.tjik <- detMats.tjik[ind.k.y.tjik]
    fy.tjik <- detMats.tji[ind.k.y.tjik2]
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

    fy.k.1.Ti <- array(fy.tik[(M * (K + 1) * (nY - 1) + 1) :
                              (M * (K + 1) * nY)],
                       c(K + 1, 1, M))

    nLL <- forward(M,nY,psi,phi.ti.prod,fy.k.1.Ti)

##     # try brute force:
##     psi.t <- array(psi, c(1, K + 1, M))
##     for(t in 1 : (nY - 1)){
##       for(i in 1 : M) {
##         psi.t[,,i] <- psi.t[,,i] %*% phi.ti.prod[,,(t-1)*M + i]
##       }
##     }

##     l.i <- numeric(M)
##     for(i in 1:M){
##       l.i[i] <- psi.t[,,i] %*% fy.k.1.Ti[,,i]
##     }

##     nLL <- -sum(log(l.i))

    print(sprintf("%i: %f",get("iteration",parent.frame(3)), nLL))
    eval.parent(quote(iteration <- iteration + 1), 3)
    nLL
  }

  fm <- optim(rep(0,nP),nll, method = "BFGS", hessian = TRUE)
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


forward <-
function(M,nY,psi,phi.ti.prod,fy.k.1.Ti)
{

    .C("forward", as.integer(M), as.integer(nY), as.double(rep(psi,M)),
   as.double(phi.ti.prod), as.double(fy.k.1.Ti),
   nLL = double(1))$nLL

}

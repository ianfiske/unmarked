markovMN <-
function(stateformula = ~ 1, detformula = ~ 1,
         data = list(y = y, covdata.site = covdata.site, covdata.obs = covdata.obs),
         y, covdata.site = NULL, covdata.obs = NULL, detconstraint = NULL,
         phiconstraint = NULL, J, K, yearly.det = TRUE)
{

  arranged <- arrangeData(data)
  cleaned <- handleNA(arranged, stateformula, detformula)
  y <- cleaned$y
  sitedata <- cleaned$covdata.site
  obsdata <- cleaned$covdata.obs

  M <- nrow(y)
  nY <- ncol(y)/J
  
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
  XDet.tjik <- matrix(rep(XDet.tji,each=K),M*J*nY*K,ncol(XDet.tji)) # test version
  k.diag <- rep(1, M * J * nY) %x% diag(K) # add intercepts the alpha_k


  if(yearly.det) {
      ## add intercepts for gma_k
      ## to keep X full column rank, use first year as baseline
      yr.int <- diag(nY - 1)
      yr.int <- rbind(rep(0, nY - 1), yr.int)
      yr.int <- yr.int %x% rep(1, M * K * J)
      
      XDet.tjik <- cbind(yr.int, k.diag, XDet.tjik)
  } else {
      XDet.tjik <- cbind(k.diag, XDet.tjik)
  }
  
  con <- detconstraint
  nDMP.un <- K*(K+1)/2
  if (is.null(con)) con <- c(1:K, rep(K + 1, nDMP.un - K))

  phicon <- phiconstraint
  # default phi constraint is different s's, equal growth, and equal reg
  #if(is.null(phicon)) phicon <- c(1:(K+1), rep(K + 2, K),
  #                                rep(K + 3, K + 2))
  nPhiP.un <- K * (K + 1)  # number of transition matrix parameters in
                                        # an unconstrained matrix
  if(is.null(phicon)) phicon <- 1:nPhiP.un
  nDMP <- max(con)  # number of independent detection matrix parameters among
                                        # the p's and beta's
  nDYP <- ifelse(yearly.det, nY - 1, 0)

  nDP <- nDCP + nDMP + nDYP # number of detection parameters including
                                        # covariates (b's), alpha's, beta's,
                                        # year effects (gamma's)
  nSP <- K                # number of parameters for psi vector of initial
                                        # latent abundances
  nPhiP <- max(phicon)    # number of independent transition matrix pars
  nP <- nDP + nSP + nPhiP  # total number of parameters

  # construct constraint matrix.  note that here we constrain beta's and
  # alphas, whereas the paper talks about constraining p's and beta's
  H.det <- matrix(0,nDMP.un,nDMP)
  for(i in 1:nDMP.un){
    H.det[i,con[i]] <- 1
  }

  # construct constrain matrix for phi paramters
  H.phi <- matrix(0, nPhiP.un, nPhiP)
  for(i in 1:nPhiP.un){
    H.phi[i, phicon[i]] <- 1
  }

  y.itj <- as.numeric(t(y))
  ## replace NA's with 99 before passing to C++
  y.itj[is.na(y.itj)] <- 99
  
  ## reorder X.tjik to be X.itjk
  t.tjik <- rep(1:nY, each = K *M*J)
  i.tjik <- rep(rep(1:M, each = K), nY*J)
  j.tjik <- rep(rep(1:J, each = K*M), nY)
  k.tjik <- rep(1:K, M * J * nY)
  XDet.itjk <- XDet.tjik[order(i.tjik, t.tjik, j.tjik, k.tjik),]

  fm <- findMLE(y.itj, XDet.itjk, nDMP, nDCP, nDP, nDYP, nSP, nPhiP, nP, nDMP.un,
         nPhiP.un, H.det, H.phi, K, yearly.det, M, J, nY)
  ests <- fm$mle
  
  DMP <- H.det %*% ests[1:nDMP]
  psi <- ests[(nDP + 1) : (nDP + K)]
  psi <- exp(c(0,psi))/sum(exp(c(0,psi)))
  phiPars <- H.phi %*% plogis(ests[(nDP + K  + 1) : nP])

  if(yearly.det) {
      gamma <- ests[(nDMP + nDCP + 1) : nDP]
  } else {
      gamma <- NULL
  }
  
  list(alpha = DMP[1:K], beta = plogis(DMP[(K+1): nDMP.un]),
       b = if(nDCP > 0) {ests[(nDMP + 1) : (nDMP + nDCP)]},
       gamma = gamma,
       psi = psi, phiPars = as.numeric(phiPars), phi = fm$phi,
       hessian = fm$hessian, AIC = 2*fm$nll + 2*nP)
}


forward <-
function(M,nY,psi,phi.ti.prod,fy.k.1.Ti)
{

    .C("forward", as.integer(M), as.integer(nY), as.double(rep(psi,M)),
   as.double(phi.ti.prod), as.double(fy.k.1.Ti),
   nLL = double(1))$nLL

}

findMLE <-
function(y.itj, XDet.itjk, nDMP, nDCP, nDP, nDYP, nSP, nPhiP, nP, nDMP.un,
         nPhiP.un, H.det, H.phi, K, yearly.det, M, J, nY)
{
    x <- .C("findMLE",
            as.integer(y.itj),
            as.double(t(XDet.itjk)),
            ncol(XDet.itjk),
            as.integer(nDMP),
            as.integer(nDCP),
            as.integer(nDP),
            as.integer(nDYP),
            as.integer(nSP),
            as.integer(nPhiP),
            as.integer(nP),
            as.integer(nDMP.un),
            as.integer(nPhiP.un),
            as.double(H.det),
            as.double(H.phi),
            as.integer(K),
            as.integer(yearly.det),
            as.integer(M),
            as.integer(J),
            as.integer(nY),
            mle = numeric(nP),
            hessian = numeric(nP^2),
            nll = numeric(1),
            phi = numeric((K + 1)^2))

    return(list(mle = x$mle, hessian = matrix(x$hessian,nP,nP),
                nll = x$nll, phi = matrix(x$phi,K + 1, K + 1)))
}

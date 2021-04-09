
#' Fit the N-mixture point count model

pcount <- function(formula, data, K, mixture = c("P", "NB", "ZIP"), starts,
                   method = "BFGS", se = TRUE,
                   engine = c("C", "R", "TMB"), threads = 1, ...)
{

    # Check argument validity--------------------------------------------------
    mixture <- match.arg(mixture, c("P", "NB", "ZIP"))
    mixture_code <- switch(mixture, P = {1}, NB = {2}, ZIP = {3})

    if(!is(data, "unmarkedFramePCount"))
        stop("Data is not an unmarkedFramePCount object.")

    engine <- match.arg(engine, c("C", "R", "TMB"))
    if(any(sapply(split_formula(formula), has_random))) engine <- "TMB"
    if(identical(mixture, "ZIP") & engine == "R")
        stop("ZIP mixture not available for R engine")

    # Generate design matrices-------------------------------------------------
    designMats <- getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    Z_state <- designMats$Z_state; Z_det <- designMats$Z_det

    # Set up parameter names and indices---------------------------------------
    lamParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nAP <- ncol(X)
    lamIdx <- 1:nAP
    pIdx <- (nAP+1):(nAP+nDP)
    n_param <- c(nAP, nDP, ifelse(mixture != "P", 1, 0))
    nP <- sum(n_param)
    nbParm <- switch(mixture, NB = "alpha", ZIP = "psi", P = character(0))

    # Handle K (number of possible abundance values to marginalize over)-------
    if(missing(K)) {
        K <- max(y, na.rm = TRUE) + 100
        warning("K was not specified and was set to ", K, ".")
    }
    if(K <= max(y, na.rm = TRUE))
        stop("specified K is too small. Try a value larger than any observation")

    #Minimum observed abundance at each site: Used by C++ and TMB
    Kmin <- apply(y, 1, function(x) max(x, na.rm=TRUE))

    # Specify negative log likelihood functions--------------------------------
    if(identical(engine, "R")) {
        k <- 0:K
        M <- nrow(y)
        J <- ncol(y)
        k.ik <- rep(k, M)
        k.ijk <- rep(k, M*J)
        y.ij <- as.numeric(t(y))
        y.ijk <- rep(y.ij, each = K + 1)
        navec <- is.na(y.ijk)
        ijk <- expand.grid(k = 0:K, j = 1:J, i = 1:M)
        ijk.to.ikj <- with(ijk, order(i, k, j))
        nll <- function(parms) {
            theta.i <- exp(X %*% parms[lamIdx] + X.offset)
            p.ij <- plogis(V %*% parms[pIdx] + V.offset)
            theta.ik <- rep(theta.i, each = K + 1)
            p.ijk <- rep(p.ij, each = K + 1)

            bin.ijk <- dbinom(y.ijk,k.ijk,p.ijk)
            bin.ijk[which(is.na(bin.ijk))] <- 1
            bin.ik.mat <- matrix(bin.ijk[ijk.to.ikj], M * (K + 1), J,
                                 byrow = TRUE)
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
    } else if(identical(engine, "C")) {
        nll <- function(parms) {
          nll_pcount(parms, n_param, y, X, V, X.offset, V.offset, K, Kmin,
                     mixture_code, threads)
        }
    }

    # Fit model in C or R------------------------------------------------------
    if(engine %in% c("C","R")){
      if(missing(starts)) starts <- rep(0, nP)
      if(length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

      fm <- optim(starts, nll, method=method, hessian=se, ...)

      ests <- fm$par
      names(ests) <- c(lamParms, detParms, nbParm)
      covMat <- invertHessian(fm, nP, se)
      fmAIC <- 2 * fm$value + 2 * nP
      tmb_mod <- NULL

      # Organize fixed-effect estimates
      state_coef <- list(ests=ests[lamIdx], cov=as.matrix(covMat[lamIdx,lamIdx]))
      det_coef <- list(ests=ests[pIdx], cov=as.matrix(covMat[pIdx, pIdx]))

      if(mixture %in% c("NB", "ZIP")){
        scale_coef <- list(ests=ests[nP], cov=as.matrix(covMat[nP,nP]))
      }

      # No random effects in C or R engines
      state_rand_info <- det_rand_info <- list()

    # Fit model in TMB---------------------------------------------------------
    } else if(engine == "TMB"){

      # Set up TMB input data
      p_form <- as.formula(formula[[2]])
      lam_form <- as.formula(paste("~", formula[3], sep=""))
      ngv_state <- get_group_vars(lam_form)
      nrand_state <- get_nrandom(lam_form, siteCovs(data))
      ngv_det <- get_group_vars(p_form)
      nrand_det <- get_nrandom(p_form, obsCovs(data))

      tmb_dat <- list(y=y, K=K, Kmin=Kmin, mixture=mixture_code,
                   X_state=X, Z_state=designMats$Z_state, offset_state=X.offset,
                   n_group_vars_state=ngv_state, n_grouplevels_state=nrand_state,
                   X_det=V, Z_det=designMats$Z_det, offset_det=V.offset,
                   n_group_vars_det=ngv_det, n_grouplevels_det=nrand_det)

      # Define TMB parameter dimensions
      tmb_param <- list(beta_state=rep(0,ncol(X)), b_state=rep(0,sum(nrand_state)),
                        lsigma_state=rep(0,ngv_state),
                        beta_det=rep(0,ncol(V)), b_det=rep(0,sum(nrand_det)),
                        lsigma_det=rep(0,ngv_det))
      if(mixture_code > 1){
        tmb_param <- c(tmb_param, list(beta_scale=0))
      }

      # Specify which parameters are random effects
      rand_ef <- NULL
      if(has_random(lam_form)) rand_ef <- c(rand_ef, "b_state")
      if(has_random(p_form)) rand_ef <- c(rand_ef, "b_det")

      # Fit model in TMB
      if(missing(starts)) starts <- NULL
      tmb_out <- fit_TMB("tmb_pcount", tmb_dat, tmb_param, rand_ef,
                         starts=starts, method, ...)
      tmb_mod <- tmb_out$TMB
      fm <- tmb_out$opt
      fmAIC <- tmb_out$AIC
      nll <- tmb_mod$fn

      # Organize fixed-effect estimate from TMB output
      state_coef <- get_coef_info(tmb_out$sdr, "state", lamParms, lamIdx)
      det_coef <- get_coef_info(tmb_out$sdr, "det", detParms, pIdx)

      if(mixture_code > 1){
        scale_coef <- get_coef_info(tmb_out$sdr, "scale", nbParm, nP)
      }

      # Organize random-effect estimates from TMB output
      state_rand_info <- get_randvar_info(tmb_out$sdr, "state",
                                          lam_form, siteCovs(data))
      det_rand_info <- get_randvar_info(tmb_out$sdr, "det", p_form, obsCovs(data))

    }

    # Create unmarkedEstimates-------------------------------------------------
    stateEstimates <- unmarkedEstimate(
        name="Abundance", short.name="lam",
        estimates = state_coef$ests, covMat = state_coef$cov, fixed=1:nAP,
        invlink = "exp", invlinkGrad = "exp",
        randomVarInfo=state_rand_info)

    detEstimates <- unmarkedEstimate(
        name = "Detection", short.name = "p",
        estimates = det_coef$ests, covMat = det_coef$cov, fixed=1:nDP,
        invlink = "logistic", invlinkGrad = "logistic.grad",
        randomVarInfo=det_rand_info)

    estimateList <- unmarkedEstimateList(list(state=stateEstimates,
                                              det=detEstimates))

    if(identical(mixture,"NB")) {
        estimateList@estimates$alpha <- unmarkedEstimate(
            name="Dispersion", short.name = "alpha",
            estimates = scale_coef$ests, covMat = scale_coef$cov, fixed=1,
            invlink = "exp", invlinkGrad = "exp", randomVarInfo=list())
    }

    if(identical(mixture,"ZIP")) {
        estimateList@estimates$psi <- unmarkedEstimate(
            name="Zero-inflation", short.name = "psi",
            estimates = scale_coef$ests, covMat = scale_coef$cov, fixed=1,
            invlink = "logistic", invlinkGrad = "logistic.grad", randomVarInfo=list())
    }

    # Create unmarkedFit object------------------------------------------------
    umfit <- new("unmarkedFitPCount", fitType="pcount", call=match.call(),
                 formula = formula, data = data,
                 sitesRemoved = designMats$removed.sites,
                 estimates = estimateList, AIC = fmAIC, opt = fm,
                 negLogLike = fm$value,
                 nllFun = nll, K = K, mixture = mixture, TMB=tmb_mod)

    return(umfit)
}

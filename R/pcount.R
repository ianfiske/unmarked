
#' Fit the N-mixture point count model

pcount <- function(formula, data, K, mixture = c("P", "NB", "ZIP"), starts,
                   method = "BFGS", se = TRUE,
                   engine = c("C", "R", "TMB"), threads = 1, ...)
{

    mixture <- match.arg(mixture, c("P", "NB", "ZIP"))
    if(!is(data, "unmarkedFramePCount"))
        stop("Data is not an unmarkedFramePCount object.")
    engine <- match.arg(engine, c("C", "R", "TMB"))
    if(any(sapply(split_formula(formula), has_random))) engine <- "TMB"
    if(identical(mixture, "ZIP") & engine %in% c("R","TMB"))
        stop("ZIP mixture not available for R or TMB engines")

    designMats <- getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if (is.null(X.offset))
        X.offset <- rep(0, nrow(X))
    if (is.null(V.offset))
        V.offset <- rep(0, nrow(V))
    NAmat <- is.na(y)

    J <- ncol(y)
    M <- nrow(y)

    lamParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nAP <- ncol(X)

    if(missing(K)) {
        K <- max(y, na.rm = TRUE) + 100
        warning("K was not specified and was set to ", K, ".")
    }
    if(K <= max(y, na.rm = TRUE))
        stop("specified K is too small. Try a value larger than any observation")
    k <- 0:K
    lk <- K+1
    M <- nrow(y)
    J <- ncol(y)
    k.ik <- rep(k, M)
    k.ijk <- rep(k, M*J)

    n_param <- c(nAP, nDP, ifelse(mixture != "P", 1, 0))
    nP <- sum(n_param)
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))
    nbParm <- switch(mixture,
                       NB = "alpha",
                       ZIP = "psi",
                       P = character(0))

    #Used by C++ and TMB
    Kmin <- apply(y, 1, function(x) max(x, na.rm=TRUE))
    Kmin <- Kmin[!1:length(Kmin) %in% designMats$removed.sites]
    mixture_code <- switch(mixture, P = {1}, NB = {2}, ZIP = {3})

    if(identical(engine, "R")) {
        y.ij <- as.numeric(t(y))
        y.ijk <- rep(y.ij, each = K + 1)
        navec <- is.na(y.ijk)
        ijk <- expand.grid(k = 0:K, j = 1:J, i = 1:M)
        ijk.to.ikj <- with(ijk, order(i, k, j))
        nll <- function(parms) {
            theta.i <- exp(X %*% parms[1 : nAP] + X.offset)
            p.ij <- plogis(V %*% parms[(nAP + 1) : (nAP + nDP)] + V.offset)
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
    } else if(identical(engine, "TMB")){

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

      tmb_param <- list(beta_state=rep(0,ncol(X)), b_state=rep(0,sum(nrand_state)),
                        lsigma_state=rep(0,ngv_state),
                        beta_det=rep(0,ncol(V)), b_det=rep(0,sum(nrand_det)),
                        lsigma_det=rep(0,ngv_det))
      if(mixture_code == 2){
        tmb_param <- c(tmb_param, list(beta_scale=0))
      }

      rand_ef <- NULL
      if(has_random(lam_form)) rand_ef <- c(rand_ef, "b_state")
      if(has_random(p_form)) rand_ef <- c(rand_ef, "b_det")

      #old_threads <- TMB::openmp()
      #on.exit(TMB::openmp(old_threads))
      #TMB::openmp(threads)

      tmb_mod <- TMB::MakeADFun(data = c(model = "tmb_pcount", tmb_dat),
                            parameters = tmb_param,
                            random= rand_ef,
                            silent=TRUE,
                            DLL = "unmarked_TMBExports")

      nfixed <- length(unlist(tmb_param[c("beta_state","beta_det",
                        "lsigma_state","lsigma_det")]))
      if(mixture_code > 1) nfixed <- nfixed + 1

      if(missing(starts)) starts <- rep(0, nfixed)
      if(length(starts) != nfixed){
        stop(paste("The number of starting values should be", nfixed))
      }
      fm <- optim(starts, fn=tmb_mod$fn, gr=tmb_mod$gr, method=method, ...)

      tmb_sum <- TMB::sdreport(tmb_mod)
      par_names <- names(tmb_sum$par.fixed)
      if(is.null(par_names)) par_names <- 1:length(tmb_sum$par.fixed)


      is_fixed <- !grepl("lsigma",par_names)
      ests <- tmb_sum$par.fixed[is_fixed]
      names(ests) <- c(lamParms, detParms, nbParm)
      covMat <- tmb_sum$cov.fixed[is_fixed,is_fixed]

      state_est <- c(ests[1:nAP], get_b_vector(tmb_mod, "state"))
      state_cov <- get_joint_cov(tmb_mod, "state")
      det_est <- c(ests[(nAP+1):(nAP+nDP)], get_b_vector(tmb_mod, "det"))
      det_cov <- get_joint_cov(tmb_mod, "det") #it is inefficient to do this twice

      if(mixture_code > 1){
        scale_est <- ests[nP]
        scale_cov <- get_joint_cov(tmb_mod, "scale")
      }

      nll <- tmb_mod$fn

      fmAIC <- 2 * fm$value + 2 * nfixed #+ 2*nP*(nP + 1)/(M - nP - 1)

      state_rand_info <- det_rand_info <- list()

      if(ngv_state > 0){
        state_sigmas <- grepl("lsigma_state", par_names)
        re_est <- tmb_sum$par.fixed[state_sigmas]
        re_names <- sigma_names(lam_form, siteCovs(data))
        re_covMat = as.matrix(tmb_sum$cov.fixed[state_sigmas,state_sigmas])

        state_rand_info <- get_randvar_info(re_names, re_est, re_covMat,
                                            lam_form, siteCovs(data))
      }
      if(ngv_det > 0){
        det_sigmas <- grepl("lsigma_det", par_names)
        re_est <- tmb_sum$par.fixed[det_sigmas]
        re_names <- sigma_names(p_form, obsCovs(data))
        re_covMat = as.matrix(tmb_sum$cov.fixed[det_sigmas,det_sigmas])

        det_rand_info <- get_randvar_info(re_names, re_est, re_covMat,
                                          p_form, obsCovs(data))
      }

    } else {
        nll <- function(parms) {
          nll_pcount(parms, n_param, y, X, V, X.offset, V.offset, K, Kmin,
                     mixture_code, threads)
        }
    }

    if(engine != "TMB"){
      if(missing(starts)) starts <- rep(0, nP)
      fm <- optim(starts, nll, method=method, hessian=se, ...)

      ests <- fm$par
      names(ests) <- c(lamParms, detParms, nbParm)
      covMat <- invertHessian(fm, nP, se)
      fmAIC <- 2 * fm$value + 2 * nP
      state_est <- ests[1:nAP]
      state_cov <- as.matrix(covMat[1:nAP,1:nAP])
      det_est <- ests[(nAP+1) : (nAP + nDP)]
      det_cov <- as.matrix(covMat[(nAP + 1):(nAP + nDP),
                                  (nAP + 1):(nAP + nDP)])
      if(mixture %in% c("NB", "ZIP")){
        scale_est <- ests[nP]
        scale_cov <- as.matrix(covMat[nP, nP])
      }
      state_rand_info <- det_rand_info <- list()
      tmb_mod <- NULL
    }

    stateName <- "Abundance"

    stateEstimates <- unmarkedEstimate(
        name=stateName, short.name="lam",
        estimates = state_est, covMat = state_cov, fixed=1:nAP,
        invlink = "exp", invlinkGrad = "exp",
        randomVarInfo=state_rand_info)

    detEstimates <- unmarkedEstimate(
        name = "Detection", short.name = "p",
        estimates = det_est, covMat = det_cov, fixed=1:nDP,
        invlink = "logistic", invlinkGrad = "logistic.grad",
        randomVarInfo=det_rand_info)

    estimateList <- unmarkedEstimateList(list(state=stateEstimates,
                                              det=detEstimates))

    if(identical(mixture,"NB")) {
        estimateList@estimates$alpha <- unmarkedEstimate(
            name="Dispersion", short.name = "alpha",
            estimates = scale_est, covMat = scale_cov, fixed=1,
            invlink = "exp", invlinkGrad = "exp", randomVarInfo=list())
    }

    if(identical(mixture,"ZIP")) {
        estimateList@estimates$psi <- unmarkedEstimate(
            name="Zero-inflation", short.name = "psi",
            estimates = scale_est, covMat = scale_cov, fixed=1,
            invlink = "logistic", invlinkGrad = "logistic.grad", randomVarInfo=list())
    }

    umfit <- new("unmarkedFitPCount", fitType="pcount", call=match.call(),
                 formula = formula, data = data,
                 sitesRemoved = designMats$removed.sites,
                 estimates = estimateList, AIC = fmAIC, opt = fm,
                 negLogLike = fm$value,
                 nllFun = nll, K = K, mixture = mixture, TMB=tmb_mod)

    return(umfit)
}

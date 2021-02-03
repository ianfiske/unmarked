
#  Fit the occupancy model of MacKenzie et al (2002).

occu <- function(formula, data, knownOcc = numeric(0),
                 linkPsi = c("logit", "cloglog"), starts, method = "BFGS",
                 se = TRUE, engine = c("C", "R", "TMB"), threads=1, ...) {

    if(!is(data, "unmarkedFrameOccu"))
        stop("Data is not an unmarkedFrameOccu object.")

    engine <- match.arg(engine, c("C", "R", "TMB"))
    if(any(sapply(split_formula(formula), has_random))) engine <- "TMB"
    if(length(knownOcc)>0 & engine == "TMB"){
      stop("TMB engine does not support knownOcc argument", call.=FALSE)
    }

    linkPsi <- match.arg(linkPsi, c("logit","cloglog"))

    designMats <- getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    removed <- designMats$removed.sites
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if(is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
    }
    if(is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
    }

    y <- truncateToBinary(y)
    J <- ncol(y)
    M <- nrow(y)

    ## convert knownOcc to logical so we can correctly to handle NAs.
    knownOccLog <- rep(FALSE, numSites(data))
    knownOccLog[knownOcc] <- TRUE
    if(length(removed)>0)
        knownOccLog <- knownOccLog[-removed]

    occParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nOP <- ncol(X)

    nP <- nDP + nOP

    yvec <- as.numeric(t(y))
    navec <- is.na(yvec)
    nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i

    ## need to add offsets !!!!!!!!!!!!!!
    ## and fix bug causing crash when NAs are in V

    linkFunc <- plogis
    invlink <- "logistic"
    linkGrad <- "logistic.grad"
    if(linkPsi == "cloglog"){
      linkFunc <- cloglog
      invlink <- "cloglog"
      linkGrad <- "cloglog.grad"
    }

    if(identical(engine, "C")) {
        nll <- function(params) {
            beta.psi <- params[1:nOP]
            beta.p <- params[(nOP+1):nP]
            .Call("nll_occu",
                  yvec, X, V, beta.psi, beta.p, nd, knownOccLog, navec,
                  X.offset, V.offset, linkPsi,
                  PACKAGE = "unmarked")
        }
    } else if(identical(engine, "TMB")){

      p_form <- as.formula(formula[[2]])
      psi_form <- as.formula(paste("~", formula[3], sep=""))
      ngv_state <- get_group_vars(psi_form)
      nrand_state <- get_nrandom(psi_form, siteCovs(data))
      ngv_det <- get_group_vars(p_form)
      nrand_det <- get_nrandom(p_form, obsCovs(data))

      tmb_dat <- list(y=y, no_detect=nd, link=ifelse(linkPsi=="cloglog",1,0),
                   X_state=X, Z_state=designMats$Z_state, offset_state=X.offset,
                   n_group_vars_state=ngv_state, n_grouplevels_state=nrand_state,
                   X_det=V, Z_det=designMats$Z_det, offset_det=V.offset,
                   n_group_vars_det=ngv_det, n_grouplevels_det=nrand_det)

      tmb_param <- list(beta_state=rep(0,ncol(X)), b_state=rep(0,sum(nrand_state)),
                  lsigma_state=rep(0,ngv_state),
                  beta_det=rep(0,ncol(V)), b_det=rep(0,sum(nrand_det)),
                  lsigma_det=rep(0,ngv_det))

      rand_ef <- NULL
      if(has_random(psi_form)) rand_ef <- c(rand_ef, "b_state")
      if(has_random(p_form)) rand_ef <- c(rand_ef, "b_det")

      old_threads <- TMB::openmp()
      on.exit(TMB::openmp(old_threads))
      TMB::openmp(threads)

      tmb_mod <- TMB::MakeADFun(data = c(model = "tmb_occu", tmb_dat),
                            parameters = tmb_param,
                            random= rand_ef,
                            silent=TRUE,
                            DLL = "unmarked_TMBExports")

      nfixed <- length(unlist(tmb_param[c("beta_state","beta_det",
                        "lsigma_state","lsigma_det")]))
      if(missing(starts)) starts <- rep(0, nfixed)
      if(length(starts) != nfixed){
        stop(paste("The number of starting values should be", nfixed))
      }
      fm <- optim(starts, fn=tmb_mod$fn, gr=tmb_mod$gr)

      tmb_sum <- TMB::sdreport(tmb_mod)
      par_names <- names(tmb_sum$par.fixed)
      if(is.null(par_names)) par_names <- 1:length(tmb_sum$par.fixed)


      is_fixed <- !grepl("lsigma",par_names)
      ests <- tmb_sum$par.fixed[is_fixed]
      names(ests) <- c(occParms, detParms)
      covMat <- tmb_sum$cov.fixed[is_fixed,is_fixed]

      state_est <- c(ests[1:nOP], get_b_vector(tmb_mod, "state"))
      state_cov <- get_joint_cov(tmb_mod, "state")
      det_est <- c(ests[(nOP+1):nP], get_b_vector(tmb_mod, "det"))
      det_cov <- get_joint_cov(tmb_mod, "det") #it is inefficient to do this twice

      nll <- tmb_mod$fn

      fmAIC <- 2 * fm$value + 2 * nfixed #+ 2*nP*(nP + 1)/(M - nP - 1)

      state_rand_info <- det_rand_info <- list()

      if(ngv_state > 0){
        state_sigmas <- grepl("lsigma_state", par_names)
        re_est <- tmb_sum$par.fixed[state_sigmas]
        re_names <- sigma_names(psi_form, siteCovs(data))
        re_covMat = as.matrix(tmb_sum$cov.fixed[state_sigmas,state_sigmas])

        state_rand_info <- get_randvar_info(re_names, re_est, re_covMat,
                                            psi_form, siteCovs(data))
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
        nll <- function(params) {
            psi <- linkFunc(X %*% params[1 : nOP] + X.offset)
            psi[knownOccLog] <- 1
            pvec <- plogis(V %*% params[(nOP + 1) : nP] + V.offset)
            cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
            cp[navec] <- 1 # so that NA's don't modify likelihood
            cpmat <- matrix(cp, M, J, byrow = TRUE) #
            loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
            -sum(loglik)
        }
    }

    if(engine != "TMB"){
      if(missing(starts)) starts <- rep(0, nP)
      if(length(starts) != nP){
        stop(paste("The number of starting values should be", nP))
      }
      fm <- optim(starts, nll, method = method, hessian = se, ...)
      covMat <- invertHessian(fm, nP, se)
      ests <- fm$par
      names(ests) <- c(occParms, detParms)
      tmb_mod <- NULL
      state_rand_info <- det_rand_info <- list()
      state_est <- ests[1:nOP]
      state_cov <- as.matrix(covMat[1:nOP,1:nOP])
      det_est <- ests[(nOP+1):nP]
      det_cov <- as.matrix(covMat[(nOP+1):nP, (nOP+1):nP])
      fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)
    }


    state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                              estimates = state_est,
                              covMat = state_cov,
                              fixed = 1:nOP,
                              invlink = invlink,
                              invlinkGrad = linkGrad,
                              randomVarInfo=state_rand_info
                              )

    det <- unmarkedEstimate(name = "Detection", short.name = "p",
                            estimates =det_est,
                            covMat = det_cov,
                            fixed = 1:nDP,
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad",
                            randomVarInfo=det_rand_info
                           )

    estimateList <- unmarkedEstimateList(list(state=state, det=det))
    umfit <- new("unmarkedFitOccu", fitType = "occu", call = match.call(),
                 formula = formula, data = data,
                 sitesRemoved = designMats$removed.sites,
                 estimates = estimateList, AIC = fmAIC, opt = fm,
                 negLogLike = fm$value,
                 nllFun = nll, knownOcc = knownOccLog, TMB=tmb_mod)

    return(umfit)
}

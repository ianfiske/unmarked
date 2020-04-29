
nmixTTD <- function(lambdaformula=~1, detformula=~1, data, K=25,
                       ttdDist=c("exp", "weibull"),
                       starts, method = "BFGS",
                       se = TRUE, engine = c("C","R"), ...) {

  #Check arguments-------------------------------------------------------------
  if(!is(data, "unmarkedFrameOccuTTD")){
    stop("Data is not an unmarkedFrameOccuTTD object.")
  }

  engine <- match.arg(engine, c("C", "R"))
  ttdDist <- match.arg(ttdDist, c("exp","weibull"))

  formula <- list(lambdaformula, ~1, ~1, detformula)
  formula <- as.formula(paste(unlist(formula),collapse=" "))

  #Process input data----------------------------------------------------------
  designMats <- getDesign(data, formula)
  #V = detection; W = abundance
  V <- designMats$V; W <- designMats$W
  y <- designMats$y
  removed <- designMats$removed.sites

  N <- nrow(y)
  R <- ncol(y)
  T <- 1 #data@numPrimary
  J <- R / T

  #Reformat data for likelihood
  yvec <- as.numeric(t(y))
  naflag <- as.numeric(is.na(yvec))
  surveyLength <- data@surveyLength
  if(length(removed>0)) surveyLength <- surveyLength[-removed,]
  ymax <- as.numeric(t(surveyLength))
  delta <- as.numeric(yvec<ymax)

  #Organize parameters---------------------------------------------------------
  detParms <- colnames(V); nDP <- ncol(V)
  abunParms <- colnames(W); nAP <- ncol(W)
  abun_inds <- 1:nAP

  det_inds <- (nAP+1):(nAP+nDP)

  parms <- c(abunParms, detParms)
  if(ttdDist == "weibull") parms <- c(parms, "k")
  nP <- length(parms)

  #Likelihood functions--------------------------------------------------------

  nll_R <- function(params){

    #Get abundance and detection parameters
    lamN <- exp(W %*% params[abun_inds])
    lamP <- exp(V %*% params[det_inds])

    pK <- sapply(0:K, function(k) dpois(k, lamN))

    #Simplified version of Garrard et al. 2013 eqn 5
    #Extended to Weibull
    if(ttdDist=='weibull'){
      shape <- exp(params[nP])

      e_lamt <- sapply(0:K, function(k){
        lam <- k*lamP
        ( shape*lam*(lam*yvec)^(shape-1) )^delta * exp(-1*(lam*yvec)^shape)
      })

    } else {
      #Exponential
      e_lamt <- sapply(0:K, function(k) (lamP*k)^delta * exp(-lamP*k*yvec))
    }

    get_Py <- function(e_lamt, delta){
      sum_delt <- as.numeric(sum(delta, na.rm=T)>0)

      out <- rep(NA, length=K+1)
    }

    #Begin likelihood calculation
    lik <- rep(NA,N)
    ystart <- 1

    for (n in 1:N){

      yend <- ystart+J-1
      pT <- rep(NA,length=K+1)
      pT[1] <- 1 - max(delta[ystart:yend], na.rm=T)
      for (k in 1:K){
        elamt_sub <- e_lamt[ystart:yend, k+1]
        pT[k+1] <- prod(elamt_sub[!is.na(elamt_sub)])
      }
      ystart <- ystart + J

      lik[n] <- pK[n,] %*% pT
    }
    -sum(log(lik))
  }

  nll_C <- function(params){
    .Call("nll_nmixTTD",
          params, yvec, delta, W, V,
          range(abun_inds)-1, range(det_inds)-1,
          ttdDist, N, J, K, naflag,
          PACKAGE = "unmarked")
  }

  nll <- nll_C
  if(engine == "R") nll <- nll_R

  #Run optim()-----------------------------------------------------------------
  if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))
  if(missing(starts)) starts <- rep(0, nP)

  fm <- optim(starts, nll, method = method, hessian = se, ...)
  if(se) {
    tryCatch(covMat <- solve(fm$hessian),
             error=function(x) stop(simpleError(paste("Hessian is singular.",
                                                      "Try providing starting values or using fewer covariates."))))
  } else {
    covMat <- matrix(NA, nP, nP)
  }

  #Build output object---------------------------------------------------------
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)
  names(ests) <- parms

  abun <- unmarkedEstimate(name = "Abundance", short.name = "lamN",
                          estimates = ests[abun_inds],
                          covMat = as.matrix(covMat[abun_inds,abun_inds]),
                          invlink = "exp",
                          invlinkGrad = "exp")

  det <- unmarkedEstimate(name = "Detection", short.name = "lamP",
                          estimates = ests[det_inds],
                          covMat = as.matrix(covMat[det_inds,det_inds]),
                          invlink = "exp",
                          invlinkGrad = "exp")


  estimateList <- unmarkedEstimateList(list(abun = abun, det=det))

  #Add Weibull shape parameter if necessary
  if(ttdDist=="weibull"){
    estimateList@estimates$shape <- unmarkedEstimate(name = "Weibull shape",
                                                     short.name = "k", estimates = ests[nP],
                                                     covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
                                                    invlinkGrad = "exp")
  }

  umfit <- new("unmarkedFitNmixTTD", fitType = "nmixTTD",
               call = match.call(),
               formula = formula,
               lambdaformula = lambdaformula,
               detformula = detformula,
               data = data, sitesRemoved = removed,
               estimates = estimateList,
               AIC = fmAIC, opt = fm, negLogLike = fm$value,
               nllFun = nll)

  return(umfit)
}

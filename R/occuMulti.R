occuMulti <- function(detformulas, stateformulas,  data, maxOrder, starts,
                      method='BFGS', se=TRUE, engine=c("C","R"), silent=FALSE, ...){

  #Format input data-----------------------------------------------------------
  #Check data object
  if(!inherits(data, "unmarkedFrameOccuMulti"))
    stop("Data must be created with unmarkedFrameOccuMulti()")

  #Check engine
  engine <- match.arg(engine, c("C", "R"))

  #Formula format
  if(missing(detformulas)){
    if(!silent){
      warning("No detection formulas specified; set to intercept-only")
    }
    detformulas <- rep('~1',length(data@ylist))
  }
  all_forms <- c(detformulas, stateformulas)
  all_forms <- all_forms[!all_forms %in% c("0","~0")]
  check_no_support(lapply(all_forms, as.formula))

  #Get design matrices and indices
  designMats <- getDesign(data, detformulas, stateformulas, maxOrder, warn=!silent)

  #Don't think there is a better way...
  N <- designMats$N; S <- designMats$S; J <- designMats$J; M <- designMats$M
  nF <- designMats$nF; nP <- designMats$nP; nOP <- designMats$nOP
  fStart <- designMats$fStart; fStop <- designMats$fStop
  dStart <- designMats$dStart; dStop <- designMats$dStop
  yStart <- designMats$yStart; yStop <- designMats$yStop
  dmF <- designMats$dmF; dmOcc <- designMats$dmOcc; dmDet <- designMats$dmDet
  y <- designMats$y; z <- designMats$z; Iy0 <- designMats$Iy0
  paramNames <- designMats$paramNames; fixed0 <- designMats$fixed0

  dmF <- Matrix::Matrix(dmF, sparse=TRUE) # convert to sparse matrix
  #Only transpose once
  t_dmF <- Matrix::t(dmF)
  #----------------------------------------------------------------------------

  #Likelihood function in R----------------------------------------------------
  nll_R <- function(params){

    #psi
    f <- matrix(NA,nrow=N,ncol=nF)
    index <- 1
    for (i in 1:nF){
      if(fixed0[i]){
        f[,i] <- 0
      } else {
        f[,i] <- dmOcc[[index]] %*% params[fStart[index]:fStop[index]]
        index <- index + 1
      }
    }
    psi <- exp(f %*% t_dmF)
    psi <- psi/rowSums(psi)

    #p
    p <- matrix(NA,nrow=nrow(y),ncol=S)
    for (i in 1:S){
      p[,i] <- plogis(dmDet[[i]] %*% params[dStart[i]:dStop[i]])
    }

    prdProbY <- matrix(NA,nrow=N,ncol=M)
    for (i in 1:N){
      inds <- yStart[i]:yStop[i]
      #column dot product of y and log(p) at site i
      cdp <- exp(colSums(y[inds,,drop=F] * log(p[inds,,drop=F])) +
           colSums((1-y[inds,,drop=F]) * log(1-p[inds,,drop=F])))
      prbSeq <- z * tcrossprod(rep(1,M),cdp) +
        (1-z) * tcrossprod(rep(1,M),Iy0[i,])
      #prod of all p at site i given occupancy states m
      prdProbY[i,] <- apply(prbSeq,1,prod)
    }

    #neg log likelihood
    -sum(log(rowSums(psi*prdProbY)))
  }
  #----------------------------------------------------------------------------

  #Likelihood function in C----------------------------------------------------
  nll_C <- function(params) {
    .Call("nll_occuMulti",
          fStart-1, fStop-1, t_dmF, dmOcc, params, dmDet, dStart-1, dStop-1,
          y, yStart-1, yStop-1, Iy0, as.matrix(z), fixed0,
          PACKAGE = "unmarked")
  }
  #----------------------------------------------------------------------------

  #Run optim()-----------------------------------------------------------------
  if(engine=="R"){
    nll <- nll_R
  } else {
    nll <- nll_C
  }

  if(missing(starts)) starts <- rep(0, nP)
  fm <- optim(starts, nll, method = method, hessian = se, ...)
  covMat <- invertHessian(fm, nP, se)

  fmAIC <- 2 * fm$value + 2 * nP
  #----------------------------------------------------------------------------

  #Format output---------------------------------------------------------------
  ests <- fm$par
  names(ests) <- paramNames

  state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                            estimates = ests[1:nOP],
                            covMat = as.matrix(covMat[1:nOP,1:nOP]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
                          estimates = ests[(nOP + 1) : nP],
                          covMat = as.matrix(covMat[(nOP + 1) : nP,
                                                      (nOP + 1) : nP]),
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  estimateList <- unmarkedEstimateList(list(state=state, det=det))

  umfit <- new("unmarkedFitOccuMulti", fitType = "occuMulti", call = match.call(),
                detformulas = detformulas, stateformulas = stateformulas,
                formula = ~1, data = data,
                #sitesRemoved = designMats$removed.sites,
                estimates = estimateList, AIC = fmAIC, opt = fm,
                negLogLike = fm$value, nllFun = nll_R)

  umfit
}

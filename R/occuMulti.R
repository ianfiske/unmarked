occuMulti <- function(detformulas, stateformulas,  data, starts,
                      method='BFGS', se=TRUE, engine="R", ...){
  
  #Format input data-----------------------------------------------------------
  #Check data object
  if(!class(data) == "unmarkedFrameOccuMulti")
    stop("Data must be created with unmarkedFrameOccuMulti()")
  
  #Check engine
  if(!engine=="R")
    stop("Only R engine is available for now")

  #Generate some indices
  S <- length(data@ylist) # of species
  z <- expand.grid(rep(list(1:0),S))[,S:1] # z matrix
  M <- nrow(z) # of possible z states
  dmF <- model.matrix(as.formula(paste0("~.^",S,"-1")),z) # f design matrix
  nF <- ncol(dmF) # of f parameters
  J <- ncol(data@ylist[[1]]) # of samples
  N <- nrow(data@ylist[[1]]) # of sites
  
  #Check formulas
  if(!(is.list(stateformulas)&is.list(detformulas)))
    stop("Formulas must be provided as lists")
  if(length(stateformulas) != nF)
    stop(paste(nF,"formulas are required in stateformulas list"))
  if(length(detformulas) != S)
    stop(paste(S,"formulas are required in detformulas list"))
  stateformulas <- lapply(stateformulas,as.formula)
  detformulas <- lapply(detformulas,as.formula)

  #Design matrices + parameter counts
  #For f/occupancy
  fInd <- c()
  dmOcc <- lapply(seq_along(stateformulas),function(i){
                    out <- model.matrix(stateformulas[[i]],data@siteCovs)
                    colnames(out) <- paste('f',i,'_',colnames(out),sep='')
                    fInd <<- c(fInd,rep(i,ncol(out)))
                    out
          })
  fStart <- c(1,1+which(diff(fInd)!=0))
  fStop <- c(fStart[2:length(fStart)]-1,length(fInd)) 
  occParams <- unlist(lapply(dmOcc,colnames))
  nOP <- length(occParams)
  
  #For detection
  dInd <- c()
  dmDet <- lapply(seq_along(detformulas),function(i){
                    out <- model.matrix(detformulas[[i]],data@obsCovs)
                    colnames(out) <- paste('sp',i,'_',colnames(out),sep='')
                    dInd <<- c(dInd,rep(i,ncol(out)))
                    out
          })
  dStart <- c(1,1+which(diff(dInd)!=0)) + nOP
  dStop <- c(dStart[2:length(dStart)]-1,length(dInd)+nOP) 
  detParams <- unlist(lapply(dmDet,colnames))
  nD <- length(detParams)
  
  #Combined
  paramNames <- c(occParams,detParams)
  nP <- length(paramNames) 

  #Re-format ylist
  index <- 1
  ylong <- lapply(data@ylist, function(x) {
                   colnames(x) <- 1:J
                   x <- cbind(x,site=1:N,species=index)
                   index <<- index+1
                   x
          })
  ylong <- as.data.frame(do.call(rbind,ylong))
  ylong <- melt(ylong,id.vars=c("site","species"),variable.name='sample')
  ylong <- dcast(ylong, site + sample ~ species)
  
  #Start-stop indices for sites
  yStart <- c(1,1+which(diff(ylong$site)!=0))
  yStop <- yStart + J - 1

  y <- as.matrix(subset(ylong,select=-c(site,sample)))

  #Indicator matrix for no detections at a site
  Iy0 <- do.call(cbind,
                 lapply(data@ylist, function(x) as.numeric(rowSums(x)==0)))
  #----------------------------------------------------------------------------

  #Likelihood function in R----------------------------------------------------
  nll_R <- function(params){
    
    #psi
    f <- matrix(NA,nrow=N,ncol=nF)
    for (i in 1:nF){
      f[,i] <- dmOcc[[i]] %*% params[fStart[i]:fStop[i]]
    }
    psi <- exp(f %*% t(dmF))
    psi <- psi/rowSums(psi)

    #p
    p <- matrix(NA,nrow=N*J,ncol=S)
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
    
    log_lik <- sum(log(rowSums(psi*prdProbY)))
    -log_lik
  }
  #----------------------------------------------------------------------------

  #Run optim()-----------------------------------------------------------------
  if(missing(starts)) starts <- rep(0, nP)
  fm <- optim(starts, nll_R, method = method, hessian = se, ...)

  if(se) {
    tryCatch(covMat <- solve(fm$hessian),
      error=function(x) stop(simpleError("Hessian is singular.
        Try providing starting values or using fewer covariates.")))
  } else { covMat <- matrix(NA, nP, nP) }

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

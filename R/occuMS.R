occuMS <- function(detformulas, stateformulas, data, starts,
                   method='BFGS', se=TRUE, engine=c("R"), silent=FALSE, ...){

  #Format input data-----------------------------------------------------------
  #Check data object
  if(!class(data) == "unmarkedFrameOccuMS")
    stop("Data must be created with unmarkedFrameOccuMS()")
 
  #Check engine
  engine <- match.arg(engine, c("R"))

  #Get design matrices and other info
  gd <- getDesign(data,stateformulas,detformulas)

  y <- gd$y
  N <- nrow(y)
  J <- ncol(y)
  S <- umf@numStates
  npsi <- S-1 #Number of free psi values 
  np <- S * (S-1) / 2 #Number of free p values
  sind <- gd$state_ind
  dind <- gd$det_ind
  nSP <- gd$nSP
  nP <- length(gd$param_names)

  #Index guide used to organize p values
  guide <- matrix(NA,nrow=S,ncol=S)
  guide <- lower.tri(guide,diag=T)
  guide[,1] <- FALSE
  guide <- which(guide,arr.ind=T) 
  #----------------------------------------------------------------------------

  #Likelihood function in R----------------------------------------------------
  #calc p(y) | p
  get_ph <- function(y, probs){

    out <- rep(1,S)
    for (j in 1:J){
      sdp <- matrix(0,nrow=S,ncol=S)
      sdp[guide] <- probs[j,]
      #sdp[guide] <- probs #constant across time
      sdp[,1] <- 1 - rowSums(sdp)
      for (s in 1:S){
        out[s] <- out[s] * sdp[s,y[j]+1]
      }
    }
    out
  }

  nll_R <- function(params){
    
    #Get psi values
    psi <- matrix(NA,nrow=N,ncol=S)
    for(i in 1:npsi){
      psi[,(i+1)] <- plogis(gd$dm_state[[i]] %*% params[sind[i,1]:sind[i,2]])
    }
    psi[,1] <- 1-apply(psi,1,sum,na.rm=T)

    #Get p values
    p <- matrix(NA,nrow=N*J,ncol=np)
    for(i in 1:np){
      p[,i] <- plogis(gd$dm_det[[i]] %*% params[dind[i,1]:dind[i,2]])
    }
  
    lik <- rep(NA,N)
    pstart <- 1
    for (n in 1:N){
      pend <- pstart+J-1
      ph <- get_ph(y[n,], p[pstart:pend,])
      lik[n] <- psi[n,] %*% ph
      pstart <- pstart + J
    }
    -sum(log(lik))
  }
  #----------------------------------------------------------------------------

  #Run optim()-----------------------------------------------------------------
  
  #Try to start params as close to 0 as possible, but negative enough that
  #we don't get initial probabilities = 0 (resulting in loglik=Inf)
  get_inits <- function(){
    out <- rep(0,nP)
    #which params are intercepts?
    ints <- c(sind[,1],dind[,1])
    
    #Calculate reasonable guess for something that will work
    fp <- S * (S-1) / 2 #Number of free p values
    val <- qlogis(1/(fp+2))

    out[ints] <- val
    out
  }
  if(missing(starts)) starts <- get_inits()
  
  fm <- suppressWarnings(
          optim(starts, nll_R, method=method, hessian = se, ...))
  
  if(se) {
    tryCatch(covMat <- solve(fm$hessian),
      error=function(x) stop(simpleError("Hessian is singular.
        Try providing starting values or using fewer covariates.")))
  } else { covMat <- matrix(NA, nP, nP) }

  fmAIC <- 2 * fm$value + 2 * nP
  #----------------------------------------------------------------------------

  #Format output---------------------------------------------------------------
  ests <- fm$par
  names(ests) <- gd$param_names

  state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                            estimates = ests[1:nSP],
                            covMat = as.matrix(covMat[1:nSP,1:nSP]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
                          estimates = ests[(nSP + 1) : nP],
                          covMat = as.matrix(covMat[(nSP + 1) : nP,
                                                      (nSP + 1) : nP]),
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  estimateList <- unmarkedEstimateList(list(state=state, det=det))
  
  umfit <- new("unmarkedFitOccuMS", fitType = "occuMS", call = match.call(),
                detformulas = detformulas, stateformulas = stateformulas,
                formula = ~1, data = data,
                sitesRemoved = gd$removed.sites,
                estimates = estimateList, AIC = fmAIC, opt = fm,
                negLogLike = fm$value, nllFun = nll_R)

  umfit
}

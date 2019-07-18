occuMS <- function(detformulas, psiformulas, phiformulas=NULL, data, 
                   parameterization=c('multinomial','condbinom'),
                   starts, method='BFGS', se=TRUE, engine=c("C","R"), 
                   silent=FALSE, ...){

  #Format input data-----------------------------------------------------------
  #Check data object
  if(!class(data) == "unmarkedFrameOccuMS")
    stop("Data must be created with unmarkedFrameOccuMS()")
  
  #Check engine
  engine <- match.arg(engine, c("C","R"))
  
  #Check parameterization
  parameterization <- match.arg(parameterization, c('multinomial','condbinom')) 
  if(parameterization=='condbinom'&data@numStates!=3){
    stop("Conditional binomial parameterization requires exactly 3 occupancy states")
  }

  #Get design matrices and other info
  gd <- getDesign(data,psiformulas,phiformulas,detformulas,parameterization)

  y <- gd$y
  N <- nrow(y)
  R <- ncol(y)
  T <- data@numPrimary
  J <- R / T
  S <- data@numStates
  npsi <- S-1 #Number of free psi values
  nphi <- S^2 - S #Number of free phi values 
  np <- S * (S-1) / 2 #Number of free p values
  sind <- gd$state_ind
  dind <- gd$det_ind
  pind <- gd$phi_ind
  nSP <- gd$nSP
  NPP <- gd$nPP
  nP <- length(gd$param_names)

  #Index guide used to organize p values
  guide <- matrix(NA,nrow=S,ncol=S)
  guide <- lower.tri(guide,diag=T)
  guide[,1] <- FALSE
  guide <- which(guide,arr.ind=T) 
  #----------------------------------------------------------------------------

  #Likelihood function in R----------------------------------------------------
  
  #Build parameter from design matrix and beta
  get_param <- function(dm_list, params, ind){
    m <- length(dm_list)
    l <- nrow(dm_list[[1]])
    out <- matrix(NA, nrow=l, ncol=m)
    for (i in 1:m){
      out[,i] <- dm_list[[i]] %*% params[ind[i,1]:ind[i,2]]
    }
    out
  }

  #calc p(y) | p
  get_ph <- function(y, probs){

    out <- rep(1,S)
    for (j in 1:J){
      if(is.na(y[j])) next

      sdp <- get_sdp(probs[j,])
      for (s in 1:S){
        out[s] <- out[s] * sdp[s,y[j]+1]
      }
    }
    out
  }

  get_psi <- function(rp, prm){
    if(prm=='multinomial'){
      # [ 1 - psi_1:psi_m, psi_1:psi_m ]
      #Multinomial logit link
      rp <- cbind(1,exp(rp))
      return(rp/rowSums(rp))

    } else if(prm=='condbinom'){
      # [ 1-psi, psi * (1-R), psi * R ]
      rp <- plogis(rp)
      return( cbind( 1-rp[,1], rp[,1]*(1-rp[,2]), rp[,1]*rp[,2] ) )
    }
  }

  get_phi_mult <- function(rp){
    rp <- matrix(rp, nrow=S)
    rp <- cbind(1,exp(rp))
    rp/rowSums(rp)
  }

  get_phi_condbin <- function(rp){
    rp <- matrix(rp, nrow=S)
    rp <- plogis(rp)
    cbind(1-rp[,1], rp[,1]*(1-rp[,2]), rp[,1]*rp[,2])
  }

  get_sdp_mult <- function(probs){
    sdp <- matrix(0,nrow=S,ncol=S)
    #Multinomial logit
    sdp[guide] <- exp(probs)
    sdp[,1] <- 1
    sdp <- sdp/rowSums(sdp)
    sdp
  }

  get_sdp_condbin <- function(probs){
    #probs order is p_1, p_2, delta
    probs <- plogis(probs)
    sdp <- matrix(0,nrow=S,ncol=S)
    sdp[1,1] <- c(1)
    sdp[2,1:2] <- c( 1-probs[1], probs[1])
    sdp[3,] <- c( 1-probs[2], probs[2]*(1-probs[3]), probs[2]*probs[3])
    sdp
  }
  
  #Set correct function to get sdp (why did I call this sdp?)
  #To save repeated conditional checks for correct function
  get_sdp <- get_sdp_mult
  get_phi <- get_phi_mult
  if(parameterization == 'condbinom'){
    get_sdp <- get_sdp_condbin
    get_phi <- get_phi_condbin
  }

  nll_R <- function(params){
    
    #Get psi values
    raw_psi <- get_param(gd$dm_state, params, sind) 
    psi <- get_psi(raw_psi, parameterization)

    if(T>1){
      raw_phi <- get_param(gd$dm_phi, params, pind)
    }

    #Get p values
    p <- get_param(gd$dm_det, params, dind)

    lik <- rep(NA,N)
    pstart <- 1
    phistart <- 1
    for (n in 1:N){
      
      phi_prod <- diag(S) 
      if(T>1){
        for (t in 1:(T-1)){
          pend <- pstart+J-1
          phiend <- phistart+nphi-1 #correct?
          D_ph <- diag(get_ph(y[n,], p[pstart:pend,]))
          phi_t <- get_phi(raw_phi[phistart:phiend])
          phi_prod <- phi_prod %*% ( D_ph %*% phi_t )
          
          phistart <- phistart + 1
          pstart <- pstart + J
        }
      } 

      pend <- pstart+J-1
      ph_T <- get_ph(y[n,], p[pstart:pend,])
      pstart <- pstart + J
          
      lik[n] <- psi[n,] %*% phi_prod %*% ph_T
    }
    -sum(log(lik))
  }
  #----------------------------------------------------------------------------
  
  #Likelihood function in C++--------------------------------------------------
  naflag <- is.na(y)
  nll_C <- function(params){
    .Call("nll_occuMS",
          params, y, gd$dm_state, gd$dm_det, sind-1, dind-1, parameterization,
          S, J, N, naflag, guide-1,
          PACKAGE = "unmarked")
  }
  #----------------------------------------------------------------------------

  #Run optim()-----------------------------------------------------------------
  
  nll <- nll_C
  #Choose function
  if(engine=="R"){
    nll <- nll_R
  }

  if(missing(starts)) starts <- rep(0, nP)
  fm <- optim(starts, nll, method=method, hessian = se, ...)
  
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
                parameterization = parameterization,
                formula = ~1, data = data,
                sitesRemoved = gd$removed.sites,
                estimates = estimateList, AIC = fmAIC, opt = fm,
                negLogLike = fm$value, nllFun = nll_R)

  umfit
}

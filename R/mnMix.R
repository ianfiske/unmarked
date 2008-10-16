#' @title Fit Multinomial-Multinomial Mixture Models
#' 
#' Fit latent abundance models to categorical data such as frog
#' calling data (Royle and Link 2005).  Currently, covariates can only be
#' supplied to the detection process (as in Royle and Link 2005).
#'
#' See \link{unmarked} for detailed descriptions of passing data \code{y},
#'   \code{covdata.site}, and \code{covdata.obs}, and specifying covariates
#'   with \code{stateformula} and \code{detformula}.
#'   
#'   This function fits the model described in Royle and Link (2005).  This
#' function assumes that a site has a latent categorical state, taking on
#' values from 0 to K.  Thus, there are K + 1 multinomial states.  Here, we
#' estimate the probabilities \eqn{\mathbf{\psi}}{\bold{psi}}, of falling
#' into the K + 1 categories.
#' 
#' Further, the latent states are observed with error, and the detection
#' process also follows a multinomial distribution according to equation
#' (3) on page 2508 of Royle and Link (2005).  The user may specify
#' covariates of the \eqn{p}'s, which are modeled as in equation (4) on
#' p. 2509.
#' 
#' The \code{constraint} argument specifies which constraints are placed on
#' the multinomial detection parameters to reduce model complexity.  The
#' same constraint vector format is used as described on page 2508 of Royle
#' and Link (2005).  The default is to allow the \eqn{p_k}'s (\eqn{k = 1,
#' 2, \ldots, K}{k = 1, 2, ..., K}) separate parameters, but constrain the
#' the \eqn{\beta}{beta}'s to be all equal.
#' 
#' @param stateformula formula for covariates of occurrance.
#' @param detformula formula for covariates of detection. 
#' @param umf unMarkedFrame supplying data.
#' @param constraint vector to describe which detection parameters are
#'     the same (see \emph{Details}}
#' @return still need to convert to UMfit
#' @export
#' @author Ian Fiske \email{ianfiske@@gmail.com}
#' @examples
#' data(gf)
#' gfUMF <- unMarkedFrame(gf.data, obsCovs = gf.obs)
#' fm.mmx1 <- mnMix(~ 1, ~ samp1 + samp3 + temp, con=c(1,2,2,3,3,3), gfUMF)
#' @keywords models
mnMix <-
function(stateformula = ~ 1, detformula = ~ 1, umf, constraint = NULL)
{

  umf <- handleNA(stateformula, detformula, umf)
  designMats <- getDesign(stateformula, detformula, umf)
  X <- designMats$X; V <- designMats$V
  y <- umf@y
  J <- ncol(y)
  M <- nrow(y)
  K <- max(y, na.rm = TRUE)
  
  nSP <- ncol(X)
  nDCP <- ncol(V)
  NParms <- colnames(X)
  detParms <- colnames(V)

  con <- constraint
  nDMP.un <- K*(K+1)/2
  if (is.null(con)) con = c(1:K, rep(K + 1, nDMP.un - K))

  # create design matrix with alpha_k intercept for each k=1,2,...,K
  nDCP <- nDCP - 1
  V <- V[,-1]
  V.ji <- V[order(rep(1:J, M), rep(1:M, each = J)),]
  V.jik <- V.ji %x% rep(1, K)  # repeat rows of X, each = K
  k.diag <- rep(1, M * J) %x% diag(K) # add intercepts the alpha_k intercepts
  V.jik <- cbind(k.diag, V.jik)

  nDMP <- max(con)    
  nDP <- nDCP + nDMP
  nSP <- K
  nP <- nDP + nSP

  # construct constraint equation
  H <- matrix(0,nDMP.un,nDMP)
  for(i in 1:nDMP.un){
    H[i,con[i]] <- 1
  }
 
  # compute indices
  det.row <- c(rep(NA,K), rep(3:(K+1), times = 1 : (K - 1))) 
  det.col <- c(rep(NA,K), sequence(1 : (K - 1)) + 1)
  arr.offset <- 0:(M * J - 1) * (K + 1)^2
  diag.els <- K*(1:K) + (2:(K+1))*2 - 1
  lower.els <- which((lower.tri(matrix(1,K+1,K+1))),arr.ind=T)
  lower.els <- lower.els[order(lower.els[,1]),] # reorder indices by row
  lower.els <- (lower.els[,2]-1)*(K+1) + lower.els[,1]  # compute vector indices
  diag.els.arr <- rep(arr.offset, each = K) + diag.els  # get offsetted indices
  lower.els.arr <- rep(arr.offset, each = sum(1:K)) + lower.els
  
  # vectorized version of detMatrix
  detMatrix <- function(dPars) {
    nmats <- nrow(dPars)#nrow(dPars)
    detMat <- lower.tri(matrix(1,(K +1),(K+1)),diag = TRUE) %x%
      array(1,c(1,1,nmats))
    # put the p's in the detMats
    detMat[diag.els.arr] <- t(dPars[, 1:K])
    detMat[lower.els.arr] <- 1 - t(dPars[,rep(1:K,times=1:K)])
    # put beta's in the mats
    for(i in (K+1):nDMP.un){
      detMat[det.row[i], det.col[i],] <- dPars[,i] *
        detMat[det.row[i], det.col[i],]
      detMat[det.row[i], 1:(det.col[i] - 1),] <- (1 - dPars[,i]) * 
        detMat[det.row[i], 1:(det.col[i] - 1),]
    }
    detMat
  }

  y.ji <- as.numeric(y)  
  y.jik <- rep(y.ji, each = K + 1)
  K.jik <- rep(0:K, M*J)
  
  nll <- function(parms) {
    # recover full parameters
    dPars <- H %*% parms[1:nDMP]

    # recover parameters
    alpha <- dPars[1 : K]
    beta <- plogis(dPars[(K + 1) : nDMP.un])
    b <- if(nDCP > 0) {parms[(nDMP + 1) : nDP]}
         else {NULL}
    psi <- parms[(nDP + 1) : nP]
    psi <- exp(c(0,psi))/sum(exp(c(0,psi)))
    
    beta.ji.mat <- matrix(beta, M*J, length(beta), byrow=TRUE)
 
    # model parms 
    p.jik <- plogis(V.jik %*% c(alpha, b))

    # get detMat paramters (rows are alphas then betas)
    p.ji.k.mat <- matrix(p.jik, M*J, K, byrow=TRUE)
    detMat.pars <- cbind(p.ji.k.mat, beta.ji.mat)

    detMats.ji <- detMatrix(detMat.pars)
    detMats.ji <- matrix(detMats.ji, M * J, (K + 1)^2, byrow = T)
    detMats.jik <- detMats.ji %x% rep(1, K + 1)
    
    rcumsums <- (K+1)^2*(0:(nrow(detMats.jik)-1))
    detMats.jik <- as.numeric(t(detMats.jik))
    f.jik <- detMats.jik[(y.jik * (K + 1) + K.jik + 1) + rcumsums]

    f.ik <- matrix(f.jik, M*(K+1), J)
    f.ik[is.na(f.ik)] <- 1
    
    psi.ik <- rep(psi,M)
    g.ik <- rowProds(f.ik)*psi.ik
    g.i <- rowSums(matrix(g.ik, M, K + 1, byrow = TRUE))
    -sum(log(g.i))
  }
  
  fm <- optim(rep(0,nP),nll, method = "BFGS")
  ests <- fm$par

  names(ests) <- c(paste(rep("detmat",nDMP),1:nDMP, sep=""), 
    eval(if(nDCP > 0) paste(rep("b", nDCP), 1:nDCP, sep="") else NULL),
    paste(rep("psi",nSP), 1:nSP, sep=""))
  list(estimates = ests, AIC = 2*fm$value + 2*nP)
}

sim2state <- function(colParms, extParms, detParms, psiParms, M, J, nY)
{
  nDP <- 1
  
  ## generate covariates
  nCov <- length(detParms) - 1
  if(nCov > 0) {
    XDet <- array(rnorm(nCov*M*J*nY), c(nCov,J,nY,M))
    detMatParms <- apply(XDet, 2:4, function(x) detParms %*% c(1,x))
  } else {
    XDet <- array(1, c(1,J,nY,M))
    detMatParms <- array(detParms, c(nDP, J, nY, M))
  }
  nSP <- length(psiParms) - 1
  W <- matrix(rnorm(nSP * M), M, nSP)
  psis <- plogis(cbind(rep(1, M), W) %*% psiParms)
  
  nPhiP <- length(colParms) - 1
  X <- array(rnorm(nPhiP*M*(nY)), c(nPhiP,nY,M))
  X2 <- X[,1:(nY-1),]
  colP <- plogis(t(apply(X, 2:3, function(x) colParms %*% c(1,x))))
  extP <- plogis(t(apply(X, 2:3, function(x) extParms %*% c(1,x))))

  ## generate first year's data
  x <- matrix(0, M, nY)
  x[,1] <- rbinom(M, 1, as.vector(psis))

  phis <- array(NA,c(2,2,nY-1,M)) #array of phis for each
  for(i in 1:M) {
    for(t in 1:(nY-1)) {
      phis[,,t,i] <- matrix(c(1-colP[i,t], colP[i,t], extP[i,t], 1-extP[i,t])) 
    }
  }
  for(i in 1:M) {
    for(t in 2:nY) {
      x[i,t] <- rbinom(1, 1, phis[2,x[i,t-1]+1,t-1,i])
    }
  }
  detMatParms <- aperm(detMatParms, c(3,1,2))
  detProbs <- plogis(detMatParms)
  
  y <- array(NA, c(M, J, nY))
  for(t in 1:nY) {
    y[,,t] <- rbinom(M*J, 1, x[,t]*detProbs[,,t])
  }
  
  XDet <- aperm(XDet,c(1,4,2,3))
  
  covdata.obs <- list()
  for(i in seq(length = dim(XDet)[1])) {
    expr <- substitute(covdata.obs$a <- matrix(XDet[i,,,], M, J*nY),
                       list(a = as.name(paste("cov",i,sep=""))))
    eval(expr)
  }
  
  y.mat <- y[,,1]
  for(i in 2:dim(y)[3]) {
    y.mat <- cbind(y.mat,y[,,i])
  }
  
  umf <- unmarkedMultFrame(y = y.mat, obsCovs = covdata.obs, numPrimary = nY, 
                           yearlySiteCovs = data.frame(cov1 = as.matrix(X)),
                           siteCovs = data.frame(cov1 = as.matrix(W)))
  
  return(list(x = matrix(x, M, nY), umf = umf))
}

test.colext.basic <- function() {

  M <- 50
  J <- 3
  nY <- 7

  psiParms <- c(0, -.2)
  colParms <- c(1,0.5)
  extParms <- c(1,0.4)
  detParms <- c(0, 0.8)

  set.seed(100)
  umf <- sim2state(colParms, extParms, detParms, psiParms, M, J, nY)$umf

  fm <- colext(~ cov1, ~cov1, ~ cov1, ~cov1, data = umf, starts = rep(0, 8))
  
  checkEqualsNumeric(coef(fm),
                     structure(c(-0.209392987585546, 0.472039779152342, 0.467454437820381, 
                                 0.399848309801933, 1.37253041645620, 0.706253195837798,
                                 -0.0149285886623921, 1.00835934154367),
                               .Names = c("psi(Int)", "psi(cov1)", "col(Int)", "col(cov1)",
                                 "ext(Int)", "ext(cov1)", "p(Int)", "p(cov1)")))

  checkEqualsNumeric(vcov(fm), 
                     structure(c(0.101821158185043, -0.0130262575534283, 0.00427456077735448, 
                                 0.00102197851578201, -0.00169134438851322, -0.000860571080969271, 
                                 -0.00388342597747938, -0.000868167091441282, -0.0130262575534283, 
                                 0.0935157485478675, -0.00166233740601025, -0.00103420458839316, 
                                 -0.00119903150370200, -0.000711971191607812, 0.000617005177732459, 
                                 7.89104003883332e-05, 0.00427456077735448, -0.00166233740601025, 
                                 0.0627554710799782, 0.0187724800045270, 0.00260881746585571, 
                                 0.00195201080484985, -0.0144021609828392, -0.0047842656701683, 
                                 0.00102197851578201, -0.00103420458839316, 0.0187724800045270, 
                                 0.0472735675186338, -0.00101763615348208, 0.00699422569135389, 
                                 -0.0051314496518188, -0.00171127300499731, -0.00169134438851322, 
                                 -0.00119903150370200, 0.00260881746585571, -0.00101763615348208, 
                                 0.071545296610352, 0.0204719778501774, 0.00254187347525365,
                                 0.000872384295660537, 
                                 -0.000860571080969271, -0.000711971191607812, 0.00195201080484985, 
                                 0.00699422569135389, 0.0204719778501774, 0.0724978390404245, 
                                 -0.00144495947502659, -0.000267114047910721, -0.00388342597747938, 
                                 0.000617005177732458, -0.0144021609828392, -0.0051314496518188, 
                                 0.00254187347525365, -0.00144495947502659, 0.0176244301586298, 
                                 0.00316532741727148, -0.000868167091441282, 7.89104003883332e-05, 
                                 -0.0047842656701683, -0.00171127300499731, 0.000872384295660537, 
                                 -0.000267114047910721, 0.00316532741727148, 0.0177794904007023
                                 ), .Dim = c(8L, 8L), .Dimnames = list(c("psi(Int)", "psi(cov1)", 
                                                        "col(Int)", "col(cov1)", "ext(Int)",
                                                        "ext(covp1)", "p(Int)", "p(cov1)"
                                                        ), c("psi(Int)", "psi(cov1)", "col(Int)",
                                                             "col(cov1)", "ext(Int)", 
                                                             "ext(cov1)", "p(Int)", "p(cov1)"))))

  checkEqualsNumeric(smoothed(fm),
                     structure(c(0.542947816209451, 0.457052183790549, 0.498135703399816, 
                                 0.501864296600184, 0.580101796528788, 0.419898203471212,
                                 0.557518381422918, 0.442481618577082, 0.624070008194787,
                                 0.375929991805213, 0.618378335592774, 0.381621664407226,
                                 0.491855090691774, 0.508144909308226),
                               .Dim = c(2L, 7L), .Dimnames = list(c("unoccupied", "occupied"),
                                                   c("1", "2", "3", "4", "5", "6", "7"))))


  checkEquals(projected(fm),
              structure(c(0.542946096641997, 0.457053903358003, 0.562649226061001, 
                          0.437350773938999, 0.558068591095918, 0.441931408904082, 0.561823040039478, 
                          0.438176959960522, 0.558873506544472, 0.441126493455528, 0.558775666826067, 
                          0.441224333173933, 0.561832468044821, 0.438167531955179),
                        .Dim = c(2L, 7L), .Dimnames = list(c("unoccupied", "occupied"),
                                            c("1", "2", "3", "4", "5", "6", "7"))))
  
}

test.colext.getP <- function() {
  checkException(p <- getP(fm))
}


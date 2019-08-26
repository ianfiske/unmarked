test.unmarkedFrameOccuTTD <- function() {

  set.seed(123)
  N <- 100
  psi <- 0.4
  lam <- 7
  Tmax <- 10

  z <- rbinom(N, 1, psi)
  y <- rexp(N, 1/lam)
  y[z==0] <- Tmax
  y[y>Tmax] <- Tmax
  
  sc <- as.data.frame(matrix(rnorm(N*2),ncol=2))
  oc <- as.data.frame(matrix(rnorm(N*2),ncol=2))

  umf <- unmarkedFrameOccuTTD(y=y, surveyLength=Tmax, siteCovs=sc, obsCovs=oc)
  
  checkEqualsNumeric(getY(umf), y) 
  checkEqualsNumeric(dim(getY(umf)), c(100,1))
  checkEqualsNumeric(siteCovs(umf), sc)
  checkEqualsNumeric(obsCovs(umf), oc)

  checkEqualsNumeric(umf@numPrimary, 1)
  checkEqualsNumeric(umf@surveyLength, matrix(Tmax, 100, 1))
  checkEquals(class(umf)[1], "unmarkedFrameOccuTTD")
  
  hd <- head(umf)
  checkEqualsNumeric(as(hd, 'data.frame'), as(umf, 'data.frame')[1:10,])
  
  umf_sub <- umf[c(1,3),]
  checkEqualsNumeric(as(umf_sub, 'data.frame'), as(umf, 'data.frame')[c(1,3),])
  
  checkException(umf[,2])

  sl_bad <- c(10,10)
  checkException(unmarkedFrameOccuTTD(y, sl_bad))

  ## Multiple observers
  y <- cbind(y,y)
  oc <- as.data.frame(matrix(rnorm(N*2*2),ncol=2))
  tm <- cbind(rep(10,N),rep(5,N))
  umf <- unmarkedFrameOccuTTD(y=y, surveyLength=tm, siteCovs=sc, obsCovs=oc)
  
  checkEqualsNumeric(getY(umf), y) 
  checkEqualsNumeric(dim(getY(umf)), c(100,2))
  checkEqualsNumeric(obsCovs(umf), oc)

  checkEqualsNumeric(umf@numPrimary, 1)
  checkEqualsNumeric(umf@surveyLength, tm)
  checkException(umf[,2])

  ## Multiple primary periods
  umf <- unmarkedFrameOccuTTD(y=y, surveyLength=tm, siteCovs=sc, 
                              yearlySiteCovs=oc, numPrimary=2)
  
  checkEqualsNumeric(yearlySiteCovs(umf), oc)
  checkEqualsNumeric(umf@numPrimary, 2)
  umf_sub <- umf[,2]
  checkEqualsNumeric(getY(umf_sub), y[,2,drop=F])
  
  y <- rexp(N, 1/lam)
  y <- cbind(y,y,y)
  checkException(unmarkedFrameOccuTTD(y,Tmax,numPrimary=2))
}

test.occuTTD.singleseason <- function(){

  #One observer------------------------------------
  set.seed(123)
  N <- 500; J <- 1

  #Simulate occupancy
  scovs <- data.frame(elev=c(scale(runif(N, 0,100))),
                          forest=runif(N,0,1),
                          wind=runif(N,0,1))
  beta_N <- c(-0.69, 0.71, -0.5)
  lambda_N <- exp(cbind(1, scovs$elev, scovs$forest) %*% beta_N)
  abun <- rpois(N, lambda_N)
  z <- as.numeric(abun>0)

  #Simulate detection
  Tmax <- 10
  beta_lam <- c(1.6, -0.2, 0.7)
  rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
  ttd <- rexp(N, 1/rate)
  ttd[z==0] <- Tmax
  ttd[ttd>Tmax] <- Tmax

  #Build UMF
  umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)

  fitR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="R")

  fitC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="C")

  checkEqualsNumeric(coef(fitR), coef(fitC))
  checkEqualsNumeric(coef(fitC), c(-0.5619,0.8624,-0.7457,1.4808,
                                   0.0211,0.4841), tol=1e-4)
  checkEqualsNumeric(coef(fitC), c(beta_N, beta_lam), tol=0.3)

  #Check weibull
  fitR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='cloglog', 
                  ttdDist='weibull',engine="R")

  fitC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='cloglog', ttdDist='weibull',
                  engine="C")

  checkEqualsNumeric(coef(fitR), coef(fitC))
  checkEqualsNumeric(coef(fitC), c(-0.6636,0.8476,-0.7021,1.2973,
                                   0.04483,0.5944,0.0979), tol=1e-4)

  #Check missing value handling
  ttd_na <- ttd; ttd_na[1] <- NA
  umf_na <- unmarkedFrameOccuTTD(y=ttd_na, Tmax, siteCovs=scovs)

  fit_naR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf_na, linkPsi='cloglog', 
                  ttdDist='weibull',engine="R")

  fit_naC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf_na, linkPsi='cloglog', ttdDist='weibull',
                  engine="C")

  checkEqualsNumeric(coef(fit_naR), coef(fit_naC))
  checkEqualsNumeric(fit_naC@AIC, 1225.399, tol=1e-4)
  checkEqualsNumeric(fit_naC@sitesRemoved, 1)

  #Two observers-----------------------------------
  set.seed(123)
  ocovs <- data.frame(obs=rep(c('A','B'),N))
  Tmax <- 10
  rateB <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam - 0.5)
  rate2 <- as.numeric(t(cbind(rate, rateB)))
  ttd <- rexp(N*2, 1/rate2)
  ttd[z==0] <- Tmax
  ttd[ttd>Tmax] <- Tmax
  ttd <- matrix(ttd, nrow=N, byrow=T)

  umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, 
                              siteCovs=scovs, obsCovs=ocovs)

  fitR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind+obs,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="R")

  fitC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind+obs,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="C")

  checkEqualsNumeric(coef(fitR), coef(fitC), tol=1e-6)
  
  #Check predict
  checkEqualsNumeric(as.numeric(predict(fitC, 'psi')[4,]), 
                     c(0.9203, 0.06757,0.381416,0.9999), tol=1e-4)
  checkEqualsNumeric(as.numeric(predict(fitC, 'det')[1,]),
                     c(29.8237, 3.1514,24.2447,36.6866), tol=1e-4)
  checkEqualsNumeric(as.numeric(predict(fitC, 'det', censor=T)[1,]),
                     c(10, NA, NA, NA))
  checkException(predict(fitC, 'col'))

  #Check getP
  gp <- getP(fitC)
  checkEqualsNumeric(dim(gp), c(500,2))
  checkEqualsNumeric(gp[1,], c(29.8237, 23.34352),tol=1e-4)

  #Check fitted
  checkException(fitted(fitC))
  #Check residuals
  checkException(residuals(fitC))

  #Check site is retained when only one observation is missing
  ttd_na <- ttd; ttd_na[1,1] <- NA
  umf_na <- unmarkedFrameOccuTTD(y=ttd_na, Tmax, siteCovs=scovs,
                                 obsCovs=ocovs)

  fit_naR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind+obs,
                  data=umf_na, linkPsi='cloglog', 
                  ttdDist='weibull',engine="R")

  fit_naC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind+obs,
                  data=umf_na, linkPsi='cloglog', ttdDist='weibull',
                  engine="C")

  checkEqualsNumeric(coef(fit_naR), coef(fit_naC))
  checkEqualsNumeric(fit_naC@sitesRemoved, numeric(0))
  
  #Check site is removed when both obs are NA
  ttd_na <- ttd; ttd_na[1,] <- NA
  umf_na <- unmarkedFrameOccuTTD(y=ttd_na, Tmax, siteCovs=scovs,
                                 obsCovs=ocovs)

  fit_naC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind+obs,
                  data=umf_na, linkPsi='cloglog', 
                  ttdDist='weibull',engine="C")
  checkEqualsNumeric(fit_naC@sitesRemoved, 1)

  #Check logit link
  set.seed(123)
  #Simulate occupancy
  scovs <- data.frame(elev=c(scale(runif(N, 0,100))),
                          forest=runif(N,0,1),
                          wind=runif(N,0,1))
  beta_N <- c(-0.69, 0.71, -0.5)
  psi <- plogis(cbind(1, scovs$elev, scovs$forest) %*% beta_N)
  z <- rbinom(N, 1, psi)

  #Simulate detection
  Tmax <- 10
  beta_lam <- c(1.6, -0.2, 0.7)
  rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
  ttd <- rexp(N, 1/rate)
  ttd[z==0] <- Tmax
  ttd[ttd>Tmax] <- Tmax

  umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)

  fitR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='logit', ttdDist='exp',engine="R")

  fitC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='logit', ttdDist='exp',engine="C")

  checkEqualsNumeric(coef(fitR), coef(fitC))
  checkEqualsNumeric(coef(fitC), c(-0.4646,0.4718,-0.8557,
                                   2.1669,-0.2539,-0.5139), tol=1e-4)

}

test.occuTTD.dynamic <- function(){

  set.seed(123)
  #Simulate initial occupancy
  N <- 5000; J <- 1; T <- 2
  scovs <- data.frame(elev=c(scale(runif(N, 0,100))),
                          forest=runif(N,0,1),
                          wind=runif(N,0,1))
  beta_psi <- c(-0.69, 0.71, -0.5)
  psi <- plogis(cbind(1, scovs$elev, scovs$forest) %*% beta_psi)
  z <- matrix(NA, N, T)
  z[,1] <- rbinom(N, 1, psi)

  #Col/ext process
  ysc <- data.frame(forest=rep(scovs$forest, each=T), 
                    elev=rep(scovs$elev, each=T))
  c_b0 <- -0.4; c_b1 <- 0.3
  gam <- plogis(c_b0 + c_b1 * scovs$forest)

  e_b0 <- -0.7; e_b1 <- 0.4
  ext <- plogis(e_b0 + e_b1 * scovs$elev)

  for (i in 1:N){
    for (t in 1:(T-1)){
      if(z[i,t]==1){
        #ext
        z[i,t+1] <- rbinom(1, 1, (1-ext[i]))
      } else {
        #col
        z[i,t+1] <- rbinom(1,1, gam[i])
      }
    }
  }

  #Simulate detection
  ocovs <- data.frame(obs=rep(c('A','B'),N*T))
  Tmax <- 10
  beta_lam <- c(1.6, -0.2, 0.7)
  rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
  rateB <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam - 0.5)
  #Across seasons
  rate2 <- as.numeric(t(cbind(rate, rateB, rate, rateB)))
  ttd <- rexp(N*T*2, 1/rate2)
  ttd <- matrix(ttd, nrow=N, byrow=T)
  ttd[ttd>Tmax] <- Tmax
  ttd[z[,1]==0,1:2] <- Tmax
  ttd[z[,2]==0,3:4] <- Tmax
  

  umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength = Tmax, 
                          siteCovs = scovs, obsCovs=ocovs,
                          yearlySiteCovs=ysc, numPrimary=2) 


  fit <- occuTTD(psiformula=~elev+forest,detformula=~elev+wind+obs,
                 gammaformula=~forest, epsilonformula=~elev, 
                 data=umf,se=T,
                 linkPsi='logit',ttdDist='exp',engine="C")

  truth <- c(beta_psi, c_b0, c_b1, e_b0, e_b1, beta_lam, -0.5)
  checkEqualsNumeric(coef(fit), truth, tol=0.1)
  checkEqualsNumeric(fit@AIC, 45521.18,tol=1e-4)

}

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
  beta_lam <- c(-2, -0.2, 0.7)
  rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
  ttd <- rexp(N, rate)
  ttd[z==0] <- Tmax
  ttd[ttd>Tmax] <- Tmax

  #Build UMF
  umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)

  fitR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="R")

  fitC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="C")

  checkEqualsNumeric(coef(fitR), coef(fitC), tol=1e-5)
  checkEqualsNumeric(coef(fitC), c(-0.6271,0.8157,-0.5982,-2.0588,
                                   -0.4042,1.2328), tol=1e-4)
  checkEqualsNumeric(coef(fitC), c(beta_N, beta_lam), tol=0.3)

  #Check weibull
  fitR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='cloglog',
                  ttdDist='weibull',engine="R")

  fitC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='cloglog', ttdDist='weibull',
                  engine="C")

  checkEqualsNumeric(coef(fitR), coef(fitC))
  checkEqualsNumeric(coef(fitC), c(-0.6846,0.7807,-0.5662,-1.9600,
                                   -0.3779,1.1474,0.06856), tol=1e-4)

  #Check missing value handling
  ttd_na <- ttd; ttd_na[1] <- NA
  umf_na <- unmarkedFrameOccuTTD(y=ttd_na, Tmax, siteCovs=scovs)

  fit_naR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf_na, linkPsi='cloglog',
                  ttdDist='weibull',engine="R")

  fit_naC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf_na, linkPsi='cloglog', ttdDist='weibull',
                  engine="C")

  checkEqualsNumeric(coef(fit_naR), coef(fit_naC), tol=1e-5)
  checkEqualsNumeric(fit_naC@AIC, 1185.15, tol=1e-4)
  checkEqualsNumeric(fit_naC@sitesRemoved, 1)

  set.seed(123)
  p <- parboot(fitC)
  checkEqualsNumeric(p@t.star[1,], 87.90704)

  r <- ranef(fitC)
  checkEqualsNumeric(dim(r@post), c(N,2,1))
  b <- bup(r)
  checkEqualsNumeric(length(b), N)

  #Two observers-----------------------------------
  set.seed(123)
  ocovs <- data.frame(obs=rep(c('A','B'),N))
  Tmax <- 10
  rateB <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam + 0.2)
  rate2 <- as.numeric(t(cbind(rate, rateB)))
  ttd <- rexp(N*2, rate2)
  ttd[ttd>Tmax] <- Tmax
  ttd <- matrix(ttd, nrow=N, byrow=T)
  ttd[z==0,] <- Tmax

  umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax,
                              siteCovs=scovs, obsCovs=ocovs)

  fitR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind+obs,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="R")

  fitC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind+obs,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="C")

  checkEqualsNumeric(coef(fitR), coef(fitC))

  #Check predict
  checkEqualsNumeric(as.numeric(predict(fitC, 'psi')[4,]),
                     c(0.7562,0.0407,0.6385,0.8588), tol=1e-4)
  checkEqualsNumeric(as.numeric(predict(fitC, 'det')[1,]),
                     c(0.16059,0.02349,0.12056,0.2139), tol=1e-4)
  checkException(predict(fitC, 'col'))

  #Check getP
  gp <- getP(fitC)
  checkEqualsNumeric(dim(gp), c(500,2))
  checkEqualsNumeric(gp[1,], c(0.79929,0.8679),tol=1e-4)

  #Check fitted
  ft <- fitted(fitC)
  checkEqualsNumeric(dim(ft),c(500,2))
  checkEqualsNumeric(ft[1,], c(0.17963,0.19505), tol=1e-4)
  #Check residuals
  r <- residuals(fitC)
  checkEqualsNumeric(dim(r),c(500,2))
  checkEqualsNumeric(r[1,], c(0.82036,0.80494), tol=1e-4)
  #Check ranef
  r <- ranef(fitC)
  checkEqualsNumeric(dim(r@post), c(N,2,1))
  b <- bup(r)
  checkEqualsNumeric(length(b), N)

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

  checkEqualsNumeric(coef(fit_naR), coef(fit_naC), tol=1e-5)
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
  beta_lam <- c(-2, -0.2, 0.7)
  rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
  ttd <- rexp(N, rate)
  ttd[z==0] <- Tmax
  ttd[ttd>Tmax] <- Tmax

  umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)

  fitR <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='logit', ttdDist='exp',engine="R")

  fitC <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
                  data=umf, linkPsi='logit', ttdDist='exp',engine="C")

  checkEqualsNumeric(coef(fitR), coef(fitC))
  checkEqualsNumeric(coef(fitC), c(-0.5102,0.67941,-0.88612,-2.1219,
                                   -0.41504,1.3224), tol=1e-4)

  # Check error when random effect in formula
  checkException(occuTTD(~(1|dummy), ~1, data=umf))
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
  beta_lam <- c(-2, -0.2, 0.7)
  rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
  rateB <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam - 0.5)
  #Across seasons
  rate2 <- as.numeric(t(cbind(rate, rateB, rate, rateB)))
  ttd <- rexp(N*T*2, rate2)
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
  checkEqualsNumeric(fit@AIC, 45037.74,tol=1e-4)

  umf_new <- umf[1:100,]

  fit <- occuTTD(psiformula=~elev+forest,detformula=~elev+wind+obs,
                 gammaformula=~forest, epsilonformula=~elev,
                 data=umf_new,se=T,
                 linkPsi='logit',ttdDist='exp',engine="C")

  s <- simulate(fit, nsim=2)
  checkEqualsNumeric(length(s),2)
  checkEqualsNumeric(dim(s[[1]]), c(100,4))
  r <- residuals(fit)
  checkEqualsNumeric(dim(r), c(100,4))

  #Check ranef
  r <- ranef(fit)
  checkEqualsNumeric(dim(r@post), c(100,2,T))
  b <- bup(r)
  checkEqualsNumeric(dim(b), c(100,T))

  # Check T > 2 works
  N <- 100; J <- 1; T <- 3
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
  beta_lam <- c(-2, -0.2, 0.7)
  rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
  rateB <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam - 0.5)
  #Across seasons
  rate2 <- as.numeric(t(cbind(rate, rateB, rate, rateB)))
  ttd <- rexp(N*T*2, rate2)
  ttd <- matrix(ttd, nrow=N, byrow=T)
  ttd[ttd>Tmax] <- Tmax
  ttd[z[,1]==0,1:2] <- Tmax
  ttd[z[,2]==0,3:4] <- Tmax

  umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength = Tmax,
                          siteCovs = scovs, obsCovs=ocovs,
                          yearlySiteCovs=ysc, numPrimary=T)

  fit2 <- occuTTD(psiformula=~elev+forest,detformula=~elev+wind+obs,
                 gammaformula=~forest, epsilonformula=~elev,
                 data=umf,se=T,
                 linkPsi='logit',ttdDist='exp',engine="C")
  checkTrue(inherits(fit2, "unmarkedFitOccuTTD"))

}

test.occuTTD.predict.complexFormulas <- function(){

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
  beta_lam <- c(-2, -0.2, 0.7)
  rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
  ttd <- rexp(N, rate)
  ttd[z==0] <- Tmax
  ttd[ttd>Tmax] <- Tmax

  #Build UMF
  umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)

  fitC <- occuTTD(psiformula=~elev+scale(forest), detformula=~elev+wind,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="C")

  nd1 <- siteCovs(umf)[1:2,]
  nd2 <- siteCovs(umf)[1:5,]
  pr1 <- predict(fitC, 'psi', newdata=nd1)$Predicted
  pr2 <- predict(fitC, 'psi', newdata=nd2)$Predicted[1:2]

  checkEqualsNumeric(pr1,pr2)

  #Check factors
  scovs$fac_cov <- factor(sample(c('a','b','c'), N, replace=T),
                          levels=c('b','a','c'))

  umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)

  fitC <- occuTTD(psiformula=~fac_cov, detformula=~elev+wind,
                  data=umf, linkPsi='cloglog', ttdDist='exp',engine="C")

  pr1 <- predict(fitC, 'psi', newdata=data.frame(fac_cov=c('b','a')))
  pr2 <- predict(fitC, 'psi', newdata=data.frame(fac_cov=c('a','b')))
  pr3 <- predict(fitC, 'psi', newdata=data.frame(fac_cov=factor(c('a','b'))))

  checkEqualsNumeric(as.matrix(pr1),as.matrix(pr2[2:1,]))
  checkEqualsNumeric(as.matrix(pr1),as.matrix(pr3[2:1,]))
  checkException(predict(fitC, 'psi', newdata=data.frame(fac_cov=c('a','d'))))

}

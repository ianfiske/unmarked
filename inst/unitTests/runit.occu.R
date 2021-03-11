test.occu.fit.simple.1 <- function() {

  y <- matrix(rep(1,10),5,2)
  umf <- unmarkedFrameOccu(y = y)
  fm <- occu(~ 1 ~ 1, data = umf)

  occ <- fm['state']
  det <- fm['det']

  occ <- coef(backTransform(occ))
  checkEqualsNumeric(occ,1)

  det <- coef(backTransform(det))
  checkEqualsNumeric(det,1)

  bt <- backTransform(fm, type = 'state')
  checkEqualsNumeric(coef(bt), 1)

  bt <- backTransform(fm, type = 'det')
  checkEqualsNumeric(coef(bt), 1)

  est_obj <- fm@estimates@estimates$state
  checkEquals(est_obj@invlink, "logistic")
  checkEquals(est_obj@invlinkGrad, "logistic.grad")
}

test.occu.fit.simple.0 <- function() {

  y <- matrix(rep(0,10),5,2)
  umf <- unmarkedFrameOccu(y = y)
  fm <- occu(~ 1 ~ 1, data = umf)

  occ <- fm['state']
  det <- fm['det']

  occ <- coef(backTransform(occ))
  checkEqualsNumeric(occ, 0, tolerance = 1e-4)

  det <- coef(backTransform(det))
  checkEqualsNumeric(det,0, tolerance = 1e-4)

  bt <- backTransform(fm, type = 'state')
  checkEqualsNumeric(coef(bt), 0, tolerance = 1e-4)

  bt <- backTransform(fm, type = 'det')
  checkEqualsNumeric(coef(bt), 0, tolerance = 1e-4)


}

test.occu.fit.covs <- function() {

  y <- matrix(rep(0:1,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ x, data = umf)

  occ <- fm['state']
  det <- fm['det']

  checkException(occ <- coef(backTransform(occ)))

  checkEqualsNumeric(coef(occ), c(8.590737, 2.472220), tolerance = 1e-4)
  checkEqualsNumeric(coef(det), c(0.44457, -0.14706, 0.44103), tolerance = 1e-4)

  occ.lc <- linearComb(fm, type = 'state', c(1, 0.5))
  det.lc <- linearComb(fm, type = 'det', c(1, 0.3, -0.3))

  checkEqualsNumeric(coef(occ.lc), 9.826848, tol = 1e-4)
  checkEqualsNumeric(coef(det.lc), 0.2681477, tol = 1e-4)

  checkEqualsNumeric(coef(backTransform(occ.lc)), 1, tol = 1e-4)
  checkEqualsNumeric(coef(backTransform(det.lc)), 0.5666381, tol = 1e-4)

  checkException(backTransform(fm, type = "state"))
  checkException(backTransform(fm, type = "det"))

  fitted <- fitted(fm)
  checkEqualsNumeric(fitted, structure(c(0.5738, 0.5014, 0.4318, 0.38581, 0.50171, 0.53764,
0.46563, 0.40283, 0.39986, 0.79928), .Dim = c(5L, 2L)), tol = 1e-5)

}

test.occu.fit.covs.0 <- function() {

  y <- matrix(rep(0,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  options(warn=2)
  checkException(fm <- occu(~ o1 + o2 ~ x, data = umf))
  options(warn=0)
  fm <- occu(~ o1 + o2 ~ x, data = umf)
  detMat <- fm@estimates@estimates$det@covMat
  stMat <- fm@estimates@estimates$state@covMat
  checkEqualsNumeric(detMat, matrix(rep(NA,9),nrow=3))
  checkEqualsNumeric(stMat, matrix(rep(NA,4),nrow=2))

}

test.occu.fit.NA <- function() {

  y <- matrix(rep(0:1,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  siteCovs[3,1] <- NA
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ x, data = umf)
  checkEquals(fm@sitesRemoved, 3)
  checkEqualsNumeric(coef(fm), c(8.70123, 4.58255, 0.66243, -0.22862, 0.58192), tol = 1e-5)

  obsCovs[10,2] <- NA
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ x, data = umf)
  checkEquals(fm@sitesRemoved, 3)
  checkEqualsNumeric(coef(fm), c(8.91289, 1.89291, -1.42471, 0.67011, -8.44608), tol = 1e-5)

}

## Add some checks here.
test.occu.offest <- function() {

  y <- matrix(rep(0:1,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ offset(x), data = umf)
  checkEqualsNumeric(coef(fm),
                     structure(c(9.74361, 0.44327, -0.14683, 0.44085), .Names = c("psi(Int)",
"p(Int)", "p(o1)", "p(o2)")), tol = 1e-5)
  fm <- occu(~ o1 + offset(o2) ~ offset(x), data = umf)
  checkEqualsNumeric(coef(fm), structure(c(8.59459, 0.97574, -0.3096), .Names = c("psi(Int)",
"p(Int)", "p(o1)")), tol=1e-5)

}

test.occu.cloglog <- function() {

  #Adapted from example by J. Cohen
  set.seed(123)
  M = 500 #sample size
  J = 3 #number of visits

  #standardized covariates
  elev <- runif(n = M, 0,100)
  forest <- runif(n = M, 0,1)
  wind <- array(runif(n = M * J, 0,20), dim = c(M, J))
  elev=as.numeric(scale(elev))
  forest=as.numeric(scale(forest))
  wind[,1] <- as.numeric(scale(wind[,1]))
  wind[,2] <- as.numeric(scale(wind[,2]))
  wind[,3] <- as.numeric(scale(wind[,3]))

  #regression parameters for abundance
  beta0 = -0.69
  beta1 = 0.71
  beta2 = -0.5

  #simulate abundance and derive true occupancy
  lambda <- exp(beta0 + beta1 * elev + beta2 * forest)
  N <- rpois(n = M, lambda = lambda)
  z <- as.numeric(N>0)
  #regression parameters for detection
  alpha0 = -0.84
  alpha1 = 2.
  alpha2 = -1.2

  #simulate detection
  p <- plogis(alpha0 + alpha1 * elev + alpha2 * wind )

  #create vectors of simulation values, for easy comparison to model estimates
  true.beta.p <- c(alpha0,alpha1,alpha2)
  true.beta.occ <- c(beta0,beta1,beta2)

  #generate observed presence
  Obs.pres <- matrix(NA,M,J)
  for (i in 1:M){
    for (j in 1:J){
      Obs.pres[i,j] <- rbinom(1,1,z[i]*p[i,j])
    }
  }
  Obs.ever <- apply(Obs.pres,1,max)

  #create observation-level covariate data frame for unmarked
  sitevec <- rep(1:M,3) #vector of site ID's
  wind.df <- data.frame("wind"=wind)
  colnames(wind.df) <- c("Wind.1","Wind.2","Wind.3")
  wind.vec <- c(wind.df$Wind.1,wind.df$Wind.2,wind.df$Wind.3)
  wind.frame <- data.frame("site"=sitevec,"wind"=wind.vec)
  wind.frame.order <- wind.frame[order(wind.frame$site),]
  wind.for.um <- data.frame(wind.frame.order$wind)
  colnames(wind.for.um)="wind"

  #create unmarked data object
  occ.frame <- unmarkedFrameOccu(Obs.pres,
                                 siteCovs=data.frame("ele"=elev,"forest"=forest),
                                  obsCovs=wind.for.um)

  #create model object
  occ_test <-occu(~ele+wind ~ele+forest, occ.frame, linkPsi="cloglog",
                      se=F)
  truth <- c(true.beta.occ, true.beta.p)
  est <- coef(occ_test)
  checkEqualsNumeric(truth, est, tol=0.1)
  checkEqualsNumeric(est,
    c(-0.7425,0.6600,-0.3333,-0.87547,2.0677,-1.3082), tol=1e-4)

  est_obj <- occ_test@estimates@estimates$state
  checkEquals(est_obj@invlink, "cloglog")
  checkEquals(est_obj@invlinkGrad, "cloglog.grad")

  #Check error if wrong link function
  checkException(occu(~ele+wind ~ele+forest, occ.frame, linkPsi="fake"))
}

test.occu.predict.complexFormulas <- function() {

  y <- matrix(rep(0:1,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ scale(o1) + o2 ~ x, data = umf)

  #Predict values should not depend on means/variance of newdata itself
  nd1 <- obsCovs(umf[1:2,])
  pr1 <- predict(fm, 'det', newdata=nd1)
  nd2 <- obsCovs(umf[1:4,])
  pr2 <- predict(fm, 'det', newdata=nd2)[1:4,]

  checkEqualsNumeric(pr1, pr2)

  #Check factors
  siteCovs$fac_cov <- factor(sample(c('a','b','c'), 5, replace=T),
                             levels=c('b','a','c'))

  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ fac_cov, data = umf)

  pr3 <- predict(fm, 'state', newdata=data.frame(fac_cov=c('a','b')))
  pr4 <- predict(fm, 'state', newdata=data.frame(fac_cov=c('b','a')))

  checkEqualsNumeric(as.matrix(pr3),as.matrix(pr4[2:1,]))
  checkException(predict(fm, 'state', newdata=data.frame(fac_cov=c('a','d'))))

  #Check when original covs contain factor not used in formula
  siteCovs$fac_cov2 <- factor(sample(c('a','b','c'), 5, replace=T),
                             levels=c('b','a','c'))

  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ fac_cov, data = umf)
  #Should error if any warnings appear
  options(warn=2)
  pr <- predict(fm, 'state', newdata=data.frame(fac_cov=c('a','b')))
  options(warn=0)

}

## Add some checks here.
test.occu.offest <- function() {

  y <- matrix(rep(0:1,10),5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ offset(x), data = umf)
  checkEqualsNumeric(coef(fm),
                     structure(c(9.74361, 0.44327, -0.14683, 0.44085), .Names = c("psi(Int)",
"p(Int)", "p(o1)", "p(o2)")), tol = 1e-5)
  fm <- occu(~ o1 + offset(o2) ~ offset(x), data = umf)
  checkEqualsNumeric(coef(fm), structure(c(8.59459, 0.97574, -0.3096), .Names = c("psi(Int)",
"p(Int)", "p(o1)")), tol=1e-5)

}

test.occu.randomeffects <- function(){
  set.seed(123)
  n_sites <- 100
  n_years <- 8
  site_id <- rep(1:n_sites, each=n_years)
  M <- n_sites * n_years
  J <- 5 # number of obs per year
  site_covs <- data.frame(cov1=rnorm(M), site_id=factor(site_id))
  beta <- c(intercept=0.5, cov1=0.3)
  sig <- 1.2
  site_effect <- rnorm(n_sites, 0, sig)
  true_site_means <- plogis(beta[1] + site_effect)

  psi <- rep(NA, M)
  for (i in 1:M){
    #Random group intercept on psi
    psi[i] <- plogis(beta[1] + beta[2]*site_covs$cov1[i]
                      + site_effect[site_id[i]])
  }

  p <- 0.5
  z <- rbinom(M, 1, psi)
  y <- matrix(0, nrow=M, ncol=J)

  for (i in 1:M){
    if(z[i]==1){
      y[i,] <- rbinom(J, 1, p)
    }
  }

  umf <- unmarkedFrameOccu(y=y, siteCovs=site_covs)
  fm <- occu(~1~cov1 + (1|site_id), umf)
  checkEqualsNumeric(coef(fm), c(0.65293, 0.39965, -0.02822), tol=1e-4)
  checkEqualsNumeric(sigma(fm)$sigma, 1.18816, tol=1e-4)

  pr <- predict(fm, "state", newdata=data.frame(cov1=0, site_id=factor(1:100)))
  checkTrue(inherits(pr, "data.frame"))

  ft <- fitted(fm)
  checkEqualsNumeric(dim(ft), c(n_sites*n_years, J))
}

context("gmultmix fitting function")

test_that("unmarkedFrameGMM construction works",{

  y <- matrix(0:3, 5, 4)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  siteCovs[3,1] <- NA
  obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/20)
  yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))

  umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
  expect_is(umf, "unmarkedFrameGMM")

  expect_error(unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="fake", numPrimary=2))


})

test_that("gmultmix removal model works",{
  y <- matrix(0:3, 5, 4)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  siteCovs[3,1] <- NA
  obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/20)
  yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))

  umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
  #fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
  expect_warning(fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C"))

  expect_equal(fm_C@sitesRemoved, 3)
  coef_truth <- c(2.50638554, 0.06226627, 0.21787839, 6.46029769, -1.51885928,
            -0.03409375, 0.43424295)
  #checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
  expect_equivalent(coef(fm_C), coef_truth, tol = 1e-5)

  # NAs in obsCovs
  obsCovs[10,2] <- NA
  umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
  #fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
  expect_warning(fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C"))

  expect_equal(fm_C@sitesRemoved, 3)
  #checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
  expect_equivalent(coef(fm_C), coef_truth, tol = 1e-5)

  # NAs in ysc
  yrSiteCovs[2, 1] <- NA
  umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
  #fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
  expect_warning(fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C"))

  coef_truth <- c(1.17280104, 0.37694710, 2.38249795, 2.87354955, -0.83875134,
            -0.08446507, 1.88056826)
  #checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
  expect_equivalent(coef(fm_C), coef_truth, tol = 1e-5)

  #Negative binomial
  #fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, mixture="NB", K=23, engine="R")
  expect_warning(fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, mixture="NB",
                                  K=23, engine="C"))
  expect_equivalent(coef(fm_C), c(1.1819, 0.3738,2.4571,4.3633,-0.8734,-0.08211,
                                  1.86049,9.38619), tol=1e-4)

  #Check methods
  expect_warning(gp <- getP(fm_C))
  expect_equal(dim(gp), c(4,4)) # missing site dropped

  expect_warning(pr <- predict(fm_C, 'lambda'))
  expect_equal(dim(pr), c(4,4))

  nd <- data.frame(x=c(0,1))
  pr <- predict(fm_C, 'lambda', newdata=nd)
  expect_equal(dim(pr), c(2,4))

  res <- residuals(fm_C)
  expect_equal(dim(res), dim(y))

  expect_warning(r <- ranef(fm_C))
  expect_equal(dim(r@post), c(4,24,1))

  expect_warning(s <- simulate(fm_C, 2))
  expect_equal(length(s), 2)

  expect_warning(pb <- parboot(fm_C, nsim=1))
  expect_is(pb, "parboot")

  expect_error(gmultmix(~(1|dummy),~1,~1,umf))

})

test_that("gmultmix double model works",{
  # Simulate independent double observer data
  nSites <- 50
  lambda <- 10
  p1 <- 0.5
  p2 <- 0.3
  cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
  set.seed(9023)
  N <- rpois(nSites, lambda)
  y <- matrix(NA, nSites, 3)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
  }

  # Fit model
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
  expect_warning(umf <- unmarkedFrameGMM(y=y, obsCovs=list(observer=observer),
         type="double",numPrimary=1))

  fm <- gmultmix(~1,~1,~observer, umf)
  expect_equivalent(coef(fm), c(2.2586,0.17385,-0.7425), tol=1e-4)

  gp <- getP(fm)
  expect_equal(dim(gp), c(nSites, 3))

})

test_that("gmultmix dependent double model works",{
  # Simulate independent double observer data
  nSites <- 50
  lambda <- 10
  p1 <- 0.5
  p2 <- 0.3
  cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
  set.seed(9023)
  N <- rpois(nSites, lambda)
  y <- matrix(NA, nSites, 3)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
  }

  # Fit model
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)

  #expect_error(umf <- unmarkedFrameGMM(y=y, obsCovs=list(observer=observer),
  #       type="depDouble",numPrimary=1))

  expect_warning(umf <- unmarkedFrameGMM(y=y[,1:2], obsCovs=list(observer=observer),
         type="depDouble",numPrimary=1))

  fm <- gmultmix(~1,~1,~observer, umf)
  expect_equivalent(coef(fm), c(1.7762,0.2493,0.2008), tol=1e-4)

  gp <- getP(fm)
  expect_equal(dim(gp), c(nSites, 2))

})

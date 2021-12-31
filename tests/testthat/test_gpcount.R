context("gpcount fitting function")

test_that("gpcount function works", {
  y <- matrix(c(0,0,0, 1,0,1, 2,2,2,
                3,2,3, 2,2,2, 1,1,1,
                NA,0,0, 0,0,0, 0,0,0,
                3,3,3, 3,1,3, 2,2,1,
                0,0,0, 0,0,0, 0,0,0), 5, 9, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,-1,4,-1))
  obsCovs <- list(o1 = matrix(seq(-3, 3, length=length(y)), 5, 9))
  obsCovs$o1[5,4:6] <- NA
  yrSiteCovs <- list(yr=matrix(c('1','2','2'), 5, 3, byrow=TRUE))
  yrSiteCovs$yr[4,2] <- NA

  expect_warning(umf <- unmarkedFrameGPC(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, numPrimary=3))

  expect_warning(fm <- gpcount(~x, ~yr, ~o1, data = umf, K=23))
  expect_equal(fm@sitesRemoved, integer(0))
  expect_equivalent(coef(fm),
        c(1.14754541, 0.44499137, -1.52079283, -0.08881542,
          2.52037155, -0.10950615), tol = 1e-5)

  # Check methods
  expect_warning(gp <- getP(fm))
  expect_equal(dim(gp), dim(y))

  expect_warning(pr <- predict(fm, 'lambda'))
  expect_equal(dim(pr), c(nrow(y), 4))

  nd <- data.frame(x=c(0,1))
  pr <- predict(fm, 'lambda', newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_equal(pr[1,1], c(3.15045), tol=1e-4)

  res <- residuals(fm)
  expect_equal(dim(res), dim(y))

  expect_warning(r <- ranef(fm))
  expect_equal(dim(r@post), c(nrow(y), 24, 1))

  expect_warning(s <- simulate(fm, 2))
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]]), dim(y))

  expect_warning(pb <- parboot(fm, nsim=1))
  expect_is(pb, "parboot")

  # Check error when random effect in formula
  expect_error(gpcount(~(1|dummy),~1,~1,umf))

})

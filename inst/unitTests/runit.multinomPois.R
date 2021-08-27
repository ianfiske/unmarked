

test.removal <- function() {

    y <- matrix(c(
        5, 3, 2,
        3, 3, 1,
        2, 0, 0,
        0, 0, 0,
        0, 0, 0), nrow=5, ncol=3, byrow=TRUE)

    sc <- data.frame(x1 = c(NA, 2, 3, 4, 3))
    oc <- list(x2 = matrix(c(
        1, 1, 1,
        3, NA, 1,
        0, 0, 1,
        NA, NA, NA,
        NA, 1, 0), nrow=5, ncol=3, byrow=TRUE))

    umf1 <- unmarkedFrameMPois(y = y, siteCovs = sc, obsCovs = oc,
        type="removal")

    o2y <- diag(ncol(y))
    o2y[upper.tri(o2y)] <- 1
    checkEquals(obsToY(umf1), o2y)

    m1_R <- multinomPois(~1 ~1, umf1, engine="R")
    m1_C <- multinomPois(~1 ~1, umf1, engine="C")
    checkEqualsNumeric(coef(m1_R), c(1.5257743, -0.2328092), tol=1e-5)
    checkEqualsNumeric(coef(m1_R), coef(m1_C), tol=1e-5)

    m2_R <- multinomPois(~x2 ~1, umf1, engine="R")
    m2_C <- multinomPois(~x2 ~1, umf1, engine="C")
    checkEqualsNumeric(coef(m2_R), c(1.9159845, 0.2248897, -0.1808144), tol=1e-5)
    checkEquals(m2_R@sitesRemoved, 4:5)
    checkEqualsNumeric(coef(m2_R),coef(m2_C), tol=1e-5)

    m3_R <- multinomPois(~x2 ~x1, umf1, engine="R")
    m3_C <- multinomPois(~x2 ~x1, umf1, engine="C")
    checkEqualsNumeric(m3_R@sitesRemoved, c(1, 4:5))
    checkEqualsNumeric(coef(m3_R),
        c(1.9118525, -0.4071202, 8.3569943, 0.3232485), tol=1e-5)
    checkEqualsNumeric(coef(m3_R),coef(m3_C), tol=1e-5)

}



test.double <- function() {
    y <- matrix(c(
        1, 0, 0,
        2, 1, 0,
        1, 0, 1,
        2, 1, 2,
        1, 0, 3,
        1, 1, 1), nrow=6, ncol=3, byrow=TRUE)
    oc <- matrix(c(
        1, 0,
        2, 1,
        1, 1,
        NA, 0,
        1, NA,
        NA, NA), nrow=6, ncol=2, byrow=TRUE)

    umf <- unmarkedFrameMPois(y = y, obsCovs = list(x=oc), type="double")

    m1_R <- multinomPois(~1 ~1, umf, engine="R")
    m1_C <- multinomPois(~1 ~1, umf, engine="C")
    checkEqualsNumeric(coef(m1_R), c(1.3137876, 0.2411609), tol=1e-5)
    checkEqualsNumeric(coef(m1_R),coef(m1_C))

    m2 <- multinomPois(~x ~1, umf, starts=c(1.3, 0, 0.2))
    checkEquals(m2@sitesRemoved, 4:6)
    }

test.ranef.multinomPois <- function(){
  set.seed(9023)
  nSites <- 50
  lambda <- 10
  p1 <- 0.5
  p2 <- 0.3
  cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
  N <- rpois(nSites, lambda)
  y <- matrix(NA, nSites, 3)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
  }

  # Fit model
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
  umf <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
            type="double")
  fm <- multinomPois(~observer-1 ~1, umf)
  checkTrue(inherits(fm, "unmarkedFitMPois"))
  checkTrue(is.null(fm@TMB))
  pr <- predict(fm, "state")
  checkEqualsNumeric(dim(pr), c(50,4))

  set.seed(1)
  nSites <- 100
  lambda <- 5
  sc <- data.frame(ref=sample(letters[1:10], nSites, replace=T),
                   x1=rnorm(nSites))
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)

  ef <- rnorm(10, 0, 0.4)
  names(ef) <- letters[1:10]
  lambda <- exp(log(lambda) + ef[sc$ref])
  N <- rpois(nSites, lambda)

  y <- matrix(NA, nSites, 3)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
  }
  umf2 <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
         type="double", siteCovs=sc)

  fm <- multinomPois(~observer-1 ~x1 + (1|ref), umf2)

  checkTrue(inherits(fm@TMB, "list"))
  checkEqualsNumeric(sigma(fm)$sigma, 0.3655, tol=1e-3)
  checkTrue(inherits(randomTerms(fm), "data.frame"))
  pr <- predict(fm, type='state')
  pr2 <- predict(fm, "state", newdata=umf2@siteCovs[1:5,])
  checkEqualsNumeric(dim(pr), c(100, 4))
  checkEqualsNumeric(dim(pr2), c(5,4))

  umf2@y[1,1] <- NA
  umf2@y[2,] <- NA
  umf2@siteCovs$x1[3] <- NA
  umf2@obsCovs$observer[80] <- NA

  fm_na <- multinomPois(~observer-1 ~x1 + (1|ref), umf2)
  checkTrue(inherits(fm_na, "unmarkedFitMPois"))

  umf3 <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
            piFun="fake", obsToY=umf@obsToY, siteCovs=sc)

  checkException(multinomPois(~observer-1 ~x1 + (1|ref), umf3))
}

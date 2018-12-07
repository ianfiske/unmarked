test.occuMulti.fit.simple.1 <- function() {
  
  y <- list(matrix(rep(1,10),5,2),
            matrix(rep(1,10),5,2))
  umf <- unmarkedFrameOccuMulti(y = y)
  fm <- occuMulti(detformulas=rep("~1",2),
                  stateformulas=rep("~1",3), data = umf, se=FALSE)
  
  occ <- predict(fm,'state')$Predicted[1,1]
  checkEqualsNumeric(occ,1, tolerance = 1e-4)

  detlist <- predict(fm,'det')
  det <- sapply(detlist,function(x) x[1,1])
  checkEqualsNumeric(det, rep(1,length(detlist)), tolerance= 1e-4)

}

test.occuMulti.fit.simple.0 <- function() {

  y <- list(matrix(rep(0,10),5,2),
            matrix(rep(0,10),5,2))
  umf <- unmarkedFrameOccuMulti(y = y)
  fm <- occuMulti(detformulas=rep("~1",2),
                  stateformulas=rep("~1",3), data = umf, se=FALSE)
  
  occ <- predict(fm,'state')$Predicted[1,1]
  checkEqualsNumeric(occ,0, tolerance = 1e-4)

  detlist <- predict(fm,'det')
  det <- sapply(detlist,function(x) x[1,1])
  checkEqualsNumeric(det, rep(0,length(detlist)), tolerance= 1e-4)


}

test.occu.fit.covs <- function() {

  y <- list(matrix(rep(0:1,10),5,2),
            matrix(rep(0:1,10),5,2))
  
  set.seed(123)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('par',1:3,sep='')
  
  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('par',1:2,sep='')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)
  stateformulas <- c('~par1','~par2','~par3')
  detformulas <- c('~par1','~par2')

  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)
 ## 
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
  checkException(fm <- occu(~ o1 + o2 ~ x, data = umf))

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

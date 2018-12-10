test.occuMulti.fit.simple.1 <- function() {
  
  y <- list(matrix(rep(1,10),5,2),
            matrix(rep(1,10),5,2))
  umf <- unmarkedFrameOccuMulti(y = y)
  fm <- occuMulti(detformulas=rep("~1",2),
                  stateformulas=rep("~1",3), data = umf, se=FALSE)
  
  #Probably should not be calling predict here b/c unit test
  #but complicated to get actual occupancy prob otherwise
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

test.occuMulti.fit.covs <- function() {

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
 
  occ <- fm['state']
  det <- fm['det']

  checkEqualsNumeric(coef(occ), c(5.36630,0.79876,5.45492,-0.868451,9.21242,1.14561), 
                     tolerance = 1e-4)
  checkEqualsNumeric(coef(det), c(-0.27586,-0.81837,-0.09537,0.42334), tolerance = 1e-4)

  fit <- fitted(fm)
  checkEqualsNumeric(length(fit),2)
  checkEqualsNumeric(sapply(fit,function(x) x[1,1]),c(0.14954,0.30801), tol = 1e-4)

  res <- residuals(fm)
  checkEqualsNumeric(length(res),2)
  checkEqualsNumeric(sapply(res,function(x) x[1,1]),c(-0.14954,-0.30801), tol= 1e-4)

}

test.occuMulti.fit.NA <- function() {
  
  y <- list(matrix(rep(0:1,10),5,2),
            matrix(rep(0:1,10),5,2))
  
  set.seed(456)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('par',1:3,sep='')
  
  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('par',1:2,sep='')

  stateformulas <- c('~par1','~par2','~par3')
  detformulas <- c('~par1','~par2')

  yna <- y
  yna[[1]][1,1] <- NA
  umf <- unmarkedFrameOccuMulti(y = yna, siteCovs = occ_covs, obsCovs = det_covs)
  
  options(warn=2)
  checkException(occuMulti(detformulas, stateformulas, data=umf, se=FALSE))
  
  options(warn=1)
  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)

  checkEqualsNumeric(coef(fm)[c(1,7)], c(6.63207,0.35323), tol= 1e-4)

  fit <- fitted(fm)
  checkTrue(is.na(fit[[1]][1,1]))

  res <- residuals(fm)
  checkTrue(is.na(res[[1]][1,1]))
  
  yna[[1]][1,] <- NA
  umf <- unmarkedFrameOccuMulti(y = yna, siteCovs = occ_covs, obsCovs = det_covs)
  checkException(occuMulti(detformulas, stateformulas, data=umf, se=FALSE))

}

test.occuMulti.fit.fixed0 <- function(){

  y <- list(matrix(rep(0:1,10),5,2),
            matrix(rep(0:1,10),5,2))
  
  set.seed(123)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('par',1:3,sep='')
  
  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('par',1:2,sep='')

  stateformulas <- c('~par1','~par2','0')
  detformulas <- c('~par1','~par2')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)
  
  occ <- fm['state']
  checkEqualsNumeric(length(coef(occ)),4)
  checkEqualsNumeric(coef(occ),c(12.26043,0.61183,12.41110,0.18764),tol=1e-4)


  stateformulas <- c('~par1','~par2')
  fm2 <- occuMulti(detformulas, stateformulas, data = umf, maxOrder=1,se=FALSE)
  
  occ <- fm2['state']
  checkEqualsNumeric(length(coef(occ)),4)
  checkEqualsNumeric(coef(occ),c(12.26043,0.61183,12.41110,0.18764),tol=1e-4)

}

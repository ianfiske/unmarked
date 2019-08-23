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

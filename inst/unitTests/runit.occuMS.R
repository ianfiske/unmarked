test.unmarkedFrameOccuMS <- function() {
  N <- 100; J <- 3; S <- 3
  psi <- c(0.5,0.3,0.2)

  p11 <- 0.4; p12 <- 0.25; p22 <- 0.3

  z <- sample(0:2, N, replace=T, prob=psi)

  y <- matrix(0,nrow=N,ncol=J)
  for (n in 1:N){
    probs <- switch(z[n]+1,
                  c(0,0,0),
                  c(1-p11,p11,0),
                  c(1-p12-p22,p12,p22))
  
    if(z[n]>0){
      y[n,] <- sample(0:2, J, replace=T, probs)
    }

  }

  site_covs <- as.data.frame(matrix(rnorm(N*2),ncol=2))
  obs_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))

  umf <- unmarkedFrameOccuMS(y=y,siteCovs=site_covs,obsCovs=obs_covs)

  checkEqualsNumeric(umf@numStates,3)
  y[y>1] <- 1
  checkException(unmarkedFrameOccuMS(y=y,siteCovs=site_covs,obsCovs=obs_covs))
}

test.occuMS.na <- function(){
  
  set.seed(123)
  N <- 10; J <- 3; S <- 3
  psi <- c(0.5,0.3,0.2)

  p11 <- 0.4; p12 <- 0.25; p22 <- 0.3

  z <- sample(0:2, N, replace=T, prob=psi)

  y <- matrix(0,nrow=N,ncol=J)
  for (n in 1:N){
    probs <- switch(z[n]+1,
                  c(0,0,0),
                  c(1-p11,p11,0),
                  c(1-p12-p22,p12,p22))
  
    if(z[n]>0){
      y[n,] <- sample(0:2, J, replace=T, probs)
    }

  }

  site_covs <- as.data.frame(matrix(rnorm(N*2),ncol=2))
  obs_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  
  umf <- unmarkedFrameOccuMS(y=y,siteCovs=site_covs,obsCovs=obs_covs)
  fit <- occuMS(rep('~1',3),rep('~1',2),umf,se=F)
  checkEqualsNumeric(fit@AIC,53.19191,tol=1e-4)

  yna <- y
  yna[1,1] <- NA
  umf <- unmarkedFrameOccuMS(y=yna,siteCovs=site_covs,obsCovs=obs_covs)
  fit <- occuMS(rep('~1',3),rep('~1',2),umf,se=F)
  checkEqualsNumeric(fit@AIC,53.06711,tol=1e-4)

  yna <- y
  yna[1,] <- NA
  sc_na <- site_covs
  sc_na[5,1] <- NA
  umf <- unmarkedFrameOccuMS(y=yna,siteCovs=sc_na,obsCovs=obs_covs)
  
  options(warn=2)
  checkException(occuMS(rep('~1',3),rep('~1',2),umf,se=F))
  options(warn=1)
  
  fit <- occuMS(rep('~1',3),rep('~1',2),umf,se=F)
  checkEqualsNumeric(fit@AIC,49.91398,tol=1e-4)
  checkEqualsNumeric(fit@sitesRemoved,c(1,5))

  oc_na <- obs_covs
  oc_na[1,1] <- NA
  umf <- unmarkedFrameOccuMS(y=y,siteCovs=site_covs,obsCovs=oc_na)

  options(warn=2)
  fit <- occuMS(rep('~1',3),rep('~1',2),umf,se=F)
  options(warn=1)

  fit <- occuMS(rep('~1',3),rep('~1',2),umf,se=F)
  checkEqualsNumeric(fit@AIC,53.06711,tol=1e-4)
}



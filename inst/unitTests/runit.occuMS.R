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

test.occuMS.multinom.fit <- function(){
  
  #Simulate data
  set.seed(123)
  N <- 50; J <- 5; S <- 3
  site_covs <- matrix(rnorm(N*2),ncol=2)
  obs_covs <- matrix(rnorm(N*J*2),ncol=2)
  a1 <- -0.5; b1 <- 1; a2 <- -0.6; b2 <- -0.7
  p11 <- -0.4; p12 <- -1.09; p22 <- -0.84
  truth <- c(a1,b1,a2,b2,p11,0,p12,p22)

  lp <- matrix(NA,ncol=S,nrow=N)
  for (n in 1:N){
    lp[n,2] <- exp(a1+b1*site_covs[n,1])
    lp[n,3] <- exp(a2+b2*site_covs[n,2])
    lp[n,1] <- 1  
  }
  psi_mat <- lp/rowSums(lp)

  z <- rep(NA,N)
  for (n in 1:N){
    z[n] <- sample(0:2, 1, replace=T, prob=psi_mat[n,])
  }

  probs_raw <- matrix(c(1,0,0,1,exp(p11),0,1,exp(p12),exp(p22)),nrow=3,byrow=T)
  probs_raw <- probs_raw/rowSums(probs_raw)
  
  y <- matrix(0,nrow=N,ncol=J)
  for (n in 1:N){

  probs <- switch(z[n]+1,
                  probs_raw[1,],
                  probs_raw[2,],
                  probs_raw[3,])
  if(z[n]>0){
    y[n,] <- sample(0:2, J, replace=T, probs)
  }
  }

  umf <- unmarkedFrameOccuMS(y=y,siteCovs=as.data.frame(site_covs),
                           obsCovs=as.data.frame(obs_covs))

  stateformulas <- c('~V1','~V2')
  detformulas <- c('~V1','~1','~1')
  fit_R <- occuMS(detformulas, stateformulas, umf, engine="R")
  fit_C <- occuMS(detformulas, stateformulas, umf, engine="C")
  checkEqualsNumeric(coef(fit_R),coef(fit_C))
  checkEqualsNumeric(coef(fit_C), c(-0.229630097798681, 0.67830519052921, 
                                    -0.0220063419144645,-0.661255952886156, 
                                    -0.554553495521214, 0.510982412286882, 
                                    -1.61783147496373, -1.50645934199995))
  #check state predict
  pr <- predict(fit_C, "state")
  checkEqualsNumeric(length(pr),2)
  checkEqualsNumeric(sapply(pr,function(x) x[1,1]),c(0.22922,0.34897),tol=1e-4)
  checkEquals(names(pr),c('psi[1]','psi[2]'))
  
  #Check bootstrapped error for predict
  checkEqualsNumeric(as.numeric(pr[[1]][1,]), 
                     c(0.2292279,0.1122459,0.07926078,0.5321636), tol=1e-4)

  #det
  pr <- predict(fit_C, "det")
  checkEqualsNumeric(length(pr),3)
  checkEqualsNumeric(sapply(pr,function(x) x[1,1]),
                     c(0.285455,0.13966,0.156119),tol=1e-4)
  checkEquals(names(pr),c('p[11]','p[12]','p[22]'))

  checkEqualsNumeric(as.numeric(pr[[1]][1,]), 
                     c(0.285455,0.069013,0.168485,0.4447024), tol=1e-4)

  #with new data (some missing)
  newdata <- data.frame(V1=rnorm(5),V2=rnorm(5))
  newdata[1,1] <- NA
  pr <- predict(fit_C,"det",newdata=newdata)
  checkTrue(is.na(pr[[1]][1,1]))
  checkEqualsNumeric(as.numeric(pr[[1]][2,]),
                     c(0.343157,0.0703713,0.222039,0.488455),tol=1e-4)

  #check getP
  ps <- getP(fit_C)
  checkEqualsNumeric(length(ps),3)
  checkEqualsNumeric(dim(ps[[1]]),c(numSites(fit_C@data),obsNum(fit_C@data)))
  checkTrue(min(unlist(ps))>=0)
  checkTrue(max(unlist(ps))<=1)
  checkEqualsNumeric(sapply(ps,function(x) x[1,1]),
                     c(0.28545,0.13966,0.156119), tol=1e-4)

  #check simulate
  set.seed(123)
  sim <- simulate(fit_C, 3)
  checkEqualsNumeric(length(sim),3)
  checkTrue(all(unlist(sim)%in%c(0:2)))
  checkEqualsNumeric(mean(fit_C@data@y),0.268)
  checkEqualsNumeric(sapply(sim,mean),c(0.244,0.280,0.288))

  #check fitted
  set.seed(123)
  fitvals <- fitted(fit_C)
  checkEqualsNumeric(dim(fitvals),c(N,J))
  checkEqualsNumeric(fitvals[1,1],0.2231388,tol=1e-4)
}


test.occuMS.condbinom.fit <- function(){
  
  #Simulate data
  set.seed(123)
  N <- 50; J <- 5; S <- 3
  site_covs <- matrix(rnorm(N*2),ncol=2)
  obs_covs <- matrix(rnorm(N*J*2),ncol=2)
  a1 <- -0.5; b1 <- 1; a2 <- -0.6; b2 <- -0.7
  p11 <- -0.4; p12 <- -1.09; p22 <- -0.84
  truth <- c(a1,b1,a2,b2,p11,0,p12,p22)

  psi_mat <- matrix(NA,ncol=S,nrow=N)
  for (n in 1:N){
    psi_mat[n,2] <- plogis(a1+b1*site_covs[n,1])
    psi_mat[n,3] <- plogis(a2+b2*site_covs[n,2])
  }

  psi_bin <- matrix(NA,nrow=nrow(psi_mat),ncol=ncol(psi_mat))
  psi_bin[,1] <- 1-psi_mat[,2]
  psi_bin[,2] <- (1-psi_mat[,3])*psi_mat[,2]
  psi_bin[,3] <- psi_mat[,2]*psi_mat[,3]

  z <- rep(NA,N)
  for (n in 1:N){
    z[n] <- sample(0:2, 1, replace=T, prob=psi_bin[n,])
  }

  p11 <- 0.4
  p12 <- 0.6
  p22 <- 0.8

  y_cb <- matrix(0,nrow=N,ncol=J)
  for (n in 1:N){
  #p11 = p1; p12 = p2; p22 = delta
    probs <- switch(z[n]+1,
                  c(1,0,0),
                  c(1-p11,p11,0),
                  c(1-p12,p12*(1-p22),p12*p22))
  
    if(z[n]>0){
      y_cb[n,] <- sample(0:2, J, replace=T, probs)
    }
  }
  truth_cb <- c(a1,b1,a2,b2,qlogis(p11),0,qlogis(c(p12,p22)))

  umf <- unmarkedFrameOccuMS(y=y_cb,siteCovs=as.data.frame(site_covs),
                           obsCovs=as.data.frame(obs_covs))

  stateformulas <- c('~V1','~V2')
  detformulas <- c('~V1','~1','~1')
  fit_R <- occuMS(detformulas, stateformulas, umf, 
                  parameterization = "condbinom", engine="R")
  fit_C <- occuMS(detformulas, stateformulas, umf, 
                  parameterization = "condbinom", engine="C")
  checkEqualsNumeric(coef(fit_R),coef(fit_C))
  checkEqualsNumeric(coef(fit_C), c(-0.5162987961667, 0.274284662180707, 
                                    -0.272563632366871, -0.85606615784698,
                                    -0.701816583657173, -0.104933853512668, 
                                    -0.21453135304912, 1.35756285443909))

  #check state predict
  pr <- predict(fit_C, "state")
  checkEqualsNumeric(length(pr),2)
  checkEqualsNumeric(as.numeric(pr[[1]][1,]), 
                     c(0.33849,0.08951,0.16304,0.51393), tol=1e-4)
  checkEquals(names(pr),c('psi','R'))
  
  #det
  pr <- predict(fit_C, "det")
  checkEqualsNumeric(length(pr),3)
  checkEqualsNumeric(as.numeric(pr[[1]][1,]), 
                     c(0.34812,0.090899,0.169970,0.526288), tol=1e-4)
  checkEquals(names(pr),c('p[1]','p[2]','delta'))

  #check getP
  ps <- getP(fit_C)
  checkEqualsNumeric(length(ps),3)
  checkEqualsNumeric(dim(ps[[1]]),c(numSites(fit_C@data),obsNum(fit_C@data)))
  checkTrue(min(unlist(ps))>=0)
  checkTrue(max(unlist(ps))<=1)
  checkEqualsNumeric(sapply(ps,function(x) x[1,1]),
                     c(0.34812,0.44657,0.79536), tol=1e-4)

  #check simulate
  set.seed(123)
  sim <- simulate(fit_C, 3)
  checkEqualsNumeric(length(sim),3)
  checkTrue(all(unlist(sim)%in%c(0:2)))
  checkEqualsNumeric(mean(fit_C@data@y),0.2)
  checkEqualsNumeric(sapply(sim,mean),c(0.200,0.156,0.128))
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
  checkException(occuMS(rep('~1',3),rep('~1',2),umf,se=F))
  options(warn=1)

  fit <- occuMS(rep('~1',3),rep('~1',2),umf,se=F)
  checkEqualsNumeric(fit@AIC,53.06711,tol=1e-4)
}



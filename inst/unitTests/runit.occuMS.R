test.unmarkedFrameOccuMS <- function() {
  set.seed(123)
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

  checkEquals(class(umf)[1], "unmarkedFrameOccuMS")
  checkEqualsNumeric(umf@numStates,3)

  umf_sub1 <- umf[1:20,]
  checkEqualsNumeric(numSites(umf_sub1),20)
  checkEquals(class(umf_sub1)[1], "unmarkedFrameOccuMS")

  y[y>1] <- 1
  checkException(unmarkedFrameOccuMS(y=y,siteCovs=site_covs,obsCovs=obs_covs))
}

test.occuMS.multinom.fit <- function(){

  #Simulate data
  set.seed(123)
  N <- 50; J <- 5; S <- 3
  site_covs <- matrix(rnorm(N*2),ncol=2)
  obs_covs <- matrix(rnorm(N*J*2),ncol=2)
  colnames(site_covs) <- paste0("sc",1:2)
  colnames(obs_covs) <- paste0("oc", 1:2)

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

  stateformulas <- c('~sc1','~sc2')
  detformulas <- c('~oc1','~1','~1')
  fit_R <- occuMS(detformulas, stateformulas, data=umf, engine="R")
  fit_C <- occuMS(detformulas, stateformulas, data=umf, engine="C")
  checkEqualsNumeric(coef(fit_R),coef(fit_C))
  checkEqualsNumeric(coef(fit_C), c(-0.229630097798681, 0.67830519052921,
                                    -0.0220063419144645,-0.661255952886156,
                                    -0.554553495521214, 0.510982412286882,
                                    -1.61783147496373, -1.50645934199995))
  #check state predict
  pr <- predict(fit_C, "psi")
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
  newdata <- data.frame(oc1=rnorm(5),oc2=rnorm(5))
  newdata[1,1] <- NA
  pr <- predict(fit_C,"det",newdata=newdata)
  checkTrue(is.na(pr[[1]][1,1]))
  checkEqualsNumeric(nrow(pr[[1]]), nrow(newdata))
  checkEqualsNumeric(as.numeric(pr[[1]][2,]),
                     c(0.343157,0.0703713,0.222039,0.488455),tol=1e-4)

  newdata <- data.frame(sc1=rnorm(5),sc2=rnorm(5))
  newdata[1,1] <- NA
  pr <- predict(fit_C,"psi",newdata=newdata)
  checkTrue(is.na(pr[[1]][1,1]))
  checkEqualsNumeric(nrow(pr[[1]]), nrow(newdata))
  checkEqualsNumeric(pr[[1]][2,1], 0.08791341,tol=1e-4)

  #With site covs in obs covs
  detformulas <- c('~sc1','~oc1','~1')
  fit2 <- occuMS(detformulas, stateformulas, data=umf)

  pr <- predict(fit2, "psi")
  checkEqualsNumeric(nrow(pr[[1]]),N)
  pr <- predict(fit2, "det")
  checkEqualsNumeric(nrow(pr[[1]]), N*J)

  pr_nd <- predict(fit2, "psi", newdata=data.frame(sc1=0, sc2=0))
  checkEqualsNumeric(nrow(pr_nd[[1]]), 1)

  pr_nd <- predict(fit2, "det", newdata=data.frame(sc1=0, oc1=0))
  checkEqualsNumeric(nrow(pr_nd[[1]]), 1)

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

  #Check fitList
  fl <- fitList(fit_C, fit_C)
  checkEquals(class(fl)[1],"unmarkedFitList")
  checkEqualsNumeric(length(fl@fits), 2)
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
  fit_R <- occuMS(detformulas, stateformulas, data=umf,
                  parameterization = "condbinom", engine="R")
  fit_C <- occuMS(detformulas, stateformulas, data=umf,
                  parameterization = "condbinom", engine="C")
  checkEqualsNumeric(coef(fit_R),coef(fit_C))
  checkEqualsNumeric(coef(fit_C), c(-0.5162987961667, 0.274284662180707,
                                    -0.272563632366871, -0.85606615784698,
                                    -0.701816583657173, -0.104933853512668,
                                    -0.21453135304912, 1.35756285443909))

  #check state predict
  pr <- predict(fit_C, "psi")
  checkEqualsNumeric(length(pr),2)
  checkEqualsNumeric(as.numeric(pr[[1]][1,]),
                     c(0.33849,0.08951,0.18945,0.52834), tol=1e-4)
  checkEquals(names(pr),c('psi','R'))

  #det
  pr <- predict(fit_C, "det")
  checkEqualsNumeric(length(pr),3)
  checkEqualsNumeric(as.numeric(pr[[1]][1,]),
                     c(0.34812,0.090899,0.195866,0.53936662), tol=1e-4)
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
  fit <- occuMS(rep('~1',3),rep('~1',2),data=umf,se=F)
  checkEqualsNumeric(fit@AIC,53.19191,tol=1e-4)

  yna <- y
  yna[1,1] <- NA
  umf <- unmarkedFrameOccuMS(y=yna,siteCovs=site_covs,obsCovs=obs_covs)
  fit <- occuMS(rep('~1',3),rep('~1',2),data=umf,se=F)
  checkEqualsNumeric(fit@AIC,53.06711,tol=1e-4)

  yna <- y
  yna[1,] <- NA
  sc_na <- site_covs
  sc_na[5,1] <- NA
  umf <- unmarkedFrameOccuMS(y=yna,siteCovs=sc_na,obsCovs=obs_covs)

  options(warn=2)
  checkException(occuMS(rep('~1',3),rep('~1',2),data=umf,se=F))
  options(warn=1)

  #Check that an NA in a variable not used in formulas doesn't drop
  #entire record
  fit <- occuMS(rep('~1',3),rep('~1',2),data=umf,se=F)
  checkEqualsNumeric(fit@AIC,51.69428,tol=1e-4)
  checkEqualsNumeric(fit@sitesRemoved, 1)

  #Now actually use variable with missing value
  fit <- occuMS(rep('~1',3), c('~V1', '~1'), data=umf,se=F)
  checkEqualsNumeric(fit@sitesRemoved, c(1,5))

  oc_na <- obs_covs
  oc_na[1,1] <- NA
  umf <- unmarkedFrameOccuMS(y=y,siteCovs=site_covs,obsCovs=oc_na)

  #Variable with missing value not used
  fit <- occuMS(rep('~1',3),rep('~1',2), data=umf,se=F)
  checkEqualsNumeric(fit@AIC,53.19191,tol=1e-4)

  #Now actually used
  options(warn=2)
  checkException(occuMS(c('~V1',rep('~1',2)),rep('~1',2), data=umf,se=F))
  options(warn=1)

  fit <- occuMS(c('~V1',rep('~1',2)),rep('~1',2),data=umf,se=F)
  checkEqualsNumeric(fit@AIC,55.03718,tol=1e-4)

  #Check that fitting works when missing site cov and no obs covs
  sc_na <- site_covs
  sc_na[1,1] <- NA
  umf <- unmarkedFrameOccuMS(y=y,siteCovs=sc_na)
  fit <- occuMS(rep('~1',3),rep('~V1',2),data=umf,se=F)
  checkEqualsNumeric(fit@sitesRemoved, 1)
}

test.occuMS.dynamic.multinomial <- function(){

set.seed(123)
N <- 100 #Number of sites
T <- 2 #Number of primary periods
J <- 3 #Number of secondary periods
S <- 3 #Number of occupancy states (0,1,2)

#Generate covariates
site_covs <- as.data.frame(matrix(rnorm(N*2),ncol=2))
yearly_site_covs <- as.data.frame(matrix(rnorm(N*T*2),ncol=2))
obs_covs <- as.data.frame(matrix(rnorm(N*J*T*2),ncol=2))

#True parameter values
b <- c(
  #Occupancy parameters
  a1=-0.5, b1=1, a2=-0.6, b2=-0.7,
  #Transition prob (phi) parameters
  phi01=0.7, phi01_cov=-0.5, phi02=-0.5, phi10=1.2,
  phi12=0.3, phi12_cov=1.1, phi20=-0.3, phi21=1.4, phi21_cov=0,
  #Detection prob parameters
  p11=-0.4, p11_cov=0, p12=-1.09, p22=-0.84
)

#Generate occupancy probs (multinomial parameterization)
lp <- matrix(1, ncol=S, nrow=N)
lp[,2] <- exp(b[1]+b[2]*site_covs[,1])
lp[,3] <- exp(b[3]+b[4]*site_covs[,2])
psi <- lp/rowSums(lp)

#True occupancy state matrix
z <- matrix(NA, nrow=N, ncol=T)

#Initial occupancy
for (n in 1:N){
  z[n,1] <- sample(0:(S-1), 1, prob=psi[n,])
}

#Raw phi probs
phi_raw <- matrix(NA, nrow=N*T, ncol=S^2-S)
phi_raw[,1] <- exp(b[5]+b[6]*yearly_site_covs[,1]) #p[0->1]
phi_raw[,2] <- exp(b[7]) #p[0->2]
phi_raw[,3] <- exp(b[8]) #p[1->0]
phi_raw[,4] <- exp(b[9]+b[10]*yearly_site_covs[,2]) #p[1->2]
phi_raw[,5] <- exp(b[11]) #p[2->0]
phi_raw[,6] <- exp(b[12]+b[13]*yearly_site_covs[,1])

#Generate states in times 2..T
px <- 1
for (n in 1:N){
  for (t in 2:T){
    phi_mat <- matrix(c(1, phi_raw[px,1], phi_raw[px,2],  # phi|z=0
                        phi_raw[px,3], 1, phi_raw[px,4],  # phi|z=1
                        phi_raw[px,5], phi_raw[px,6], 1), # phi|z=2
                      nrow=S, byrow=T)
    phi_mat <- phi_mat/rowSums(phi_mat)
    z[n, t] <- sample(0:(S-1), 1, prob=phi_mat[z[n,(t-1)]+1,])
    px <- px + 1
    if(t==T) px <- px + 1 #skip last datapoint for each site
  }
}

#Raw p probs
p_mat <- matrix(c(1, 0, 0, #p|z=0
                  1, exp(b[14]), 0, #p|z=1
                  1, exp(b[16]), exp(b[17])), #p|z=2
                nrow=S, byrow=T)
p_mat <- p_mat/rowSums(p_mat)

#Simulate observation data
y <- matrix(0, nrow=N, ncol=J*T)
for (n in 1:N){
  yx <- 1
  for (t in 1:T){
    if(z[n,t]==0){
      yx <- yx + J
      next
    }
    for (j in 1:J){
      y[n, yx] <- sample(0:(S-1), 1, prob=p_mat[z[n,t]+1,])
      yx <- yx+1
    }
  }
}

umf <- unmarkedFrameOccuMS(y=y, siteCovs=site_covs,
                           obsCovs=obs_covs,
                           yearlySiteCovs=yearly_site_covs,
                           numPrimary=T)

psiformulas <- c('~V1','~V2') #on psi[1] and psi[2]
phiformulas <- c('~V1','~1','~1','~V2','~1','~V1')
detformulas <- c('~V1','~1','~1') #on p[1|1], p[1|2], p[2|2]

fitC <- occuMS(detformulas=detformulas, psiformulas=psiformulas,
              phiformulas=phiformulas, data=umf)

fitR <- occuMS(detformulas=detformulas, psiformulas=psiformulas,
              phiformulas=phiformulas, data=umf,engine="R")

checkEqualsNumeric(fitC@AIC,799.1723,tol=1e-4)
checkEqualsNumeric(fitC@AIC,fitR@AIC,tol=1e-4)
checkEqualsNumeric(length(coef(fitC)),17)

phiformulas_new <- rep('~1',6)
fit_new <- update(fitC,phiformulas=phiformulas_new)
checkEqualsNumeric(fit_new@AIC,800.8553,tol=1e-4)
checkEqualsNumeric(length(coef(fit_new)),14)

set.seed(123)
fit_sim <- simulate(fitC,nsim=2)
checkEqualsNumeric(fit_sim[[1]][2,],c(0,2,1,0,0,2))

pr_phi <- predict(fitC,'phi')
pr_phi <- sapply(pr_phi, function(x) x$Predicted[1])
checkEqualsNumeric(pr_phi,
                   c(0.055117,0.57195,0.931102,0.02733675,
                     0.2192255,0.1819882),tol=1e-4)

#Check predicting phi with newdata works
nd <- data.frame(V1=c(1,2), V2=1)
pr_nd <- predict(fitC, "phi", newdata=nd)
checkEqualsNumeric(length(pr_nd), 6)
checkEqualsNumeric(sapply(pr_nd, nrow), rep(2,6))
#No differing covariates on either phi value so estimates should be the same
checkEqualsNumeric(pr_nd[[3]][1,1], pr_nd[[3]][2,1])

umf_new <- umf
umf_new@y[1,1:3] <- NA

options(warn=2)
checkException(
  occuMS(detformulas=detformulas, psiformulas=psiformulas,
              phiformulas=phiformulas, data=umf_new)
)
options(warn=1)

umf_miss <- umf
umf_miss@yearlySiteCovs[1,1] <- NA
checkException(
  occuMS(detformulas=detformulas, psiformulas=psiformulas,
              phiformulas=phiformulas, data=umf_miss)
)

}

test.occuMS.dynamic.condbinom <- function(){

set.seed(123)
N <- 100
J <- 3
S <- 3
T <- 2

site_covs <- matrix(rnorm(N*2),ncol=2)
obs_covs <- matrix(rnorm(N*J*T*2),ncol=2)
yearly_site_covs <- matrix(rnorm(N*T*2),ncol=2)

a1 <- -0.5; b1 <- 1; a2 <- -0.6; b2 <- -0.7
phi0 <- 0.7; phi1 <- -0.5; phi2 <- 1.2
R0 <- 0.3; R1 <- -0.3; R2 <- 0.5
p11 <- 0.4; p12 <- 0.6; p22 <- 0.8

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

y_cb <- matrix(0,nrow=N,ncol=J*T)
z_new <- rep(NA,N)
for (n in 1:N){

  phi <- matrix(NA,nrow=S,ncol=S)
  phi[1,] <- c(1-plogis(phi0),plogis(phi0)*(1-plogis(R0)),plogis(phi0)*plogis(R0))
  phi[2,] <- c(1-plogis(phi1),plogis(phi1)*(1-plogis(R1)),plogis(phi1)*plogis(R1))
  phi[3,] <- c(1-plogis(phi2),plogis(phi2)*(1-plogis(R2)),plogis(phi2)*plogis(R2))


  #p11 = p1; p12 = p2; p22 = delta
  probs <- switch(z[n]+1,
                  c(1,0,0),
                  c(1-p11,p11,0),
                  c(1-p12,p12*(1-p22),p12*p22))

  if(z[n]>0){
    y_cb[n,1:J] <- sample(0:2, J, replace=T, probs)
  }

  for (i in 2:T){
    phi_t <- switch(z[n]+1,
                    phi[1,],phi[2,],phi[3,])
    z_new[n] <- sample(0:2, 1, prob=phi_t)
    if(z_new[n]>0){
      probs <- switch(z_new[n]+1,
                      c(1,0,0),
                      c(1-p11,p11,0),
                      c(1-p12,p12*(1-p22),p12*p22))
      y_cb[n,((T-1)*J+1):(T*J)] <- sample(0:2, J, replace=T, probs)
    }
  }

}

umf <- unmarkedFrameOccuMS(y=y_cb,siteCovs=as.data.frame(site_covs),
                           obsCovs=as.data.frame(obs_covs),
                           #yearlySiteCovs=as.data.frame(yearly_site_covs),
                           numPrimary=2)

stateformulas <- c('~V1','~V2')
detformulas <- c('~V1','~1','~1')
phiformulas<- rep('~1',6)

fit_cbC <- occuMS(detformulas=detformulas, psiformulas=stateformulas,
              phiformulas=phiformulas,
              parameterization='condbinom',
              data=umf, se=T,engine="C")

checkEqualsNumeric(length(coef(fit_cbC)),14)
checkEqualsNumeric(fit_cbC@AIC,820.0645,tol=1e-4)

fit_cbR <- occuMS(detformulas=detformulas, psiformulas=stateformulas,
              phiformulas=phiformulas,
              parameterization='condbinom',
              data=umf, se=T,engine="R")
checkEqualsNumeric(fit_cbC@AIC,fit_cbR@AIC,tol=1e-4)

set.seed(123)
fit_sim <- simulate(fit_cbC,nsim=2)
checkEqualsNumeric(fit_sim[[1]][1,],c(0,0,0,2,1,0))

pr_phi <- predict(fit_cbC,'phi')
pr_phi <- sapply(pr_phi, function(x) x$Predicted[1])
checkEqualsNumeric(pr_phi,
                   c(0.72966,0.54682,0.728597,0.3856194,
                     0.999950,0.5841778),tol=1e-4)

}

test.occuMS.predit.complexFormulas <- function(){

  #Simulate data
  set.seed(123)
  N <- 50; J <- 5; S <- 3
  site_covs <- matrix(rnorm(N*2, mean=2),ncol=2)
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

  stateformulas <- c('~scale(V1)','~V2')
  detformulas <- c('~V1','~1','~1')
  fit_C <- occuMS(detformulas, stateformulas, data=umf, engine="C")

  #Check with newdata; contents of newdata should not
  #effect resulting predictions (scale should be based on
  #original data)
  nd <- siteCovs(umf)[1:5,]
  pr_nd <- predict(fit_C, type='psi', newdata=nd, se=F)$Predicted
  nd <- siteCovs(umf)[1:2,]
  pr_nd2 <- predict(fit_C, type='psi', newdata=nd, se=F)$Predicted
  nd <- siteCovs(umf)[c(1,1),]
  pr_nd3 <- predict(fit_C, type='psi', newdata=nd, se=F)$Predicted

  checkEqualsNumeric(pr_nd[1:2,], pr_nd2)
  checkEqualsNumeric(pr_nd[c(1,1),], pr_nd3)

  #Check for factor level handling
  site_covs2 <- as.data.frame(site_covs)
  site_covs2$occ_fac <- factor(sample(c('a','b','c'),N,replace=T),
                              levels=c('b','a','c'))

  umf <- unmarkedFrameOccuMS(y=y,siteCovs=site_covs2,
                           obsCovs=as.data.frame(obs_covs))
  stateformulas <- c('~occ_fac','~1')
  fm <- occuMS(detformulas, stateformulas, data = umf)

  nd <- siteCovs(umf)[1:2,]
  pr_nd <- predict(fm, type='psi', newdata=nd, se=F)$Predicted

  nd2 <- data.frame(occ_fac=factor(c('a','b'),levels=c('b','a','c')))
  pr_nd2 <- predict(fm, type='psi', newdata=nd2, se=F)$Predicted

  checkEqualsNumeric(pr_nd, pr_nd2[c(2,1),])

  nd3 <- data.frame(occ_fac=c('a','b'))
  pr_nd3 <- predict(fm, type='psi', newdata=nd3, se=F)$Predicted

  checkEqualsNumeric(pr_nd, pr_nd3[c(2,1),])

  nd4 <- data.frame(occ_fac=c('a','d'))
  checkException(predict(fm, type='psi', newdata=nd4, se=F))

}

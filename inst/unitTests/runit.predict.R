
test.occu <- function() {

    if(!require(raster))
        stop("raster package required")
    set.seed(55)
    R <- 20
    J <- 4
    x1 <- rnorm(R)
    x2 <- factor(c(rep('A', R/2), rep('B', R/2)))
    x3 <- matrix(rnorm(R*J), R, J)
    z <- rbinom(R, 1, 0.5)
    y <- matrix(rbinom(R*J, 1, z*0.6), R, J)
    x1[1] <- NA
    x3[2,1] <- NA
    x3[3,] <- NA
    umf1 <- unmarkedFrameOccu(y=y, siteCovs=data.frame(x1=x1, x2=x2),
                              obsCovs=list(x3=x3))
    fm1 <- occu(~x3 ~x1+x2, umf1)
    E1.1 <- predict(fm1, type="state")
    E1.2 <- predict(fm1, type="det")

    nd1.1 <- data.frame(x1=0, x2=factor('A', levels=c('A','B')))
    nd1.2 <- data.frame(x3=0)
    E1.3 <- predict(fm1, type="state", newdata=nd1.1, appendData=TRUE)
    E1.4 <- predict(fm1, type="det", newdata=nd1.2)

    r1 <- raster(matrix(rnorm(100), 10))
    checkException(predict(fm1, type="state", newdata=r1))
    s1 <- stack(r1)
    checkException(predict(fm1, type="state", newdata=s1))
    names(s1) <- c("x3")
    E1.5 <- predict(fm1, type="det", newdata=s1)
    E1.5 <- predict(fm1, type="det", newdata=s1, appendData=TRUE)

    E1.6 <- predict(fm1, type="state", level=0.9)
    checkEquals(as.numeric(E1.6[1,3:4]), c(0.01881844, 0.8538048))

}



test.pcount <- function() {

    set.seed(55)
    R <- 20
    J <- 4
    N <- rpois(R, 2)
    y <- matrix(rbinom(R*J, N, 0.7), R, J)
    umf1 <- unmarkedFramePCount(y=y)

    fm1 <- pcount(~1 ~1, umf1, K=40)
    E1.1 <- predict(fm1, type="state")
    E1.2 <- predict(fm1, type="det")

    fm2 <- pcount(~1 ~1, umf1, K=40, mixture="NB")
    E2.1 <- predict(fm2, type="state")
    checkException(predict(fm2, type="alpha"))

    fm3 <- pcount(~1 ~1, umf1, K=40, mixture="ZIP")
    E3.1 <- predict(fm3, type="state")
    checkException(predict(fm3, type="psi"))
    checkEquals(E3.1[1,1], 1.818512, tol=1e-6)

}





test.pcountOpen <- function() {

    set.seed(55)
    R <- 10
    J <- 4
    T <- 3
    N <- matrix(NA, R, T)
    N[,1] <- rpois(R, 4)
    N[,2] <- rbinom(R, N[,1], 0.8) + rpois(R, 1)
    N[,3] <- rbinom(R, N[,2], 0.8) + rpois(R, 1)
    y1 <- matrix(rbinom(R*J, N[,1], 0.7), R, J)
    y2 <- matrix(rbinom(R*J, N[,2], 0.7), R, J)
    y3 <- matrix(rbinom(R*J, N[,3], 0.7), R, J)
    umf1 <- unmarkedFramePCO(y=cbind(y1,y2,y3), numPrimary=T)

#    fm1 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=30)
#    E1.1 <- predict(fm1, type="lambda")
#    E1.2 <- predict(fm1, type="det")

    fm2 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=40, mixture="NB")
    checkException(predict(fm2, type="alpha"))

    fm3 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=40, mixture="ZIP")
    E3.1 <- predict(fm3, type="lambda")
    checkException(predict(fm3, type="psi"))
    checkEquals(E3.1[1,1], 2.472029, tol=1e-6)

    #With newdata
    siteCovs(umf1) <- data.frame(x=rnorm(10))
    fm4 <- pcountOpen(~x, ~1, ~1, ~1, umf1, K=40, mixture="NB")
    nd <- data.frame(x=1)
    checkEqualsNumeric(predict(fm4, "lambda", newdata=nd)$Predicted,
                       2.427701, tol=1e-6)

}

test.occuMulti <- function() {

  set.seed(123)
  N <- 10
  nspecies <- 2
  J <- 5

  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- list()
  for (i in 1:nspecies){
    det_covs[[i]] <- matrix(rnorm(N*J),nrow=N)
  }
  names(det_covs) <- paste('det_cov',1:nspecies,sep='')

  #True vals
  beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
  f1 <- beta[1] + beta[2]*occ_covs$occ_cov1
  f2 <- beta[3] + beta[4]*occ_covs$occ_cov2
  f3 <- beta[5] + beta[6]*occ_covs$occ_cov3
  f <- cbind(f1,f2,f3)
  z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
  colnames(z) <- paste('sp',1:nspecies,sep='')
  dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)

  psi <- exp(f %*% t(dm))
  psi <- psi/rowSums(psi)

  #True state
  ztruth <- matrix(NA,nrow=N,ncol=nspecies)
  for (i in 1:N){
    ztruth[i,] <- as.matrix(z[sample(4,1,prob=psi[i,]),])
  }

  p_true <- c(0.6,0.7)

  # fake y data
  y <- list()

  for (i in 1:nspecies){
    y[[i]] <- matrix(NA,N,J)
    for (j in 1:N){
      for (k in 1:J){
        y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p_true[i])
      }
    }
  }
  names(y) <- c('coyote','tiger')

  #Create the unmarked data object
  data = unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)

  occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov3')
  detFormulas <- c('~1','~1')

  fit <- occuMulti(detFormulas,occFormulas,data)

  pr_state <- predict(fit,'state')
  checkEqualsNumeric(pr_state$Predicted[1,],
                     c(0.34936,0.18146,0.21590,0.25327),
                     tol=1e-4)
  checkEqualsNumeric(pr_state$SE[1,],
                     c(0.2023365,0.1334475,0.2009201,0.1551536),
                     tol=1e-4)
  pr_det <- predict(fit,'det')
  checkEqualsNumeric(length(pr_det),nspecies)
  checkEqualsNumeric(sapply(pr_det,function(x) x[1,1]),c(0.59429,0.64731),tol=1e-4)

  #marginal occupancy
  pr_marg <- predict(fit,'state',species=2)
  checkEqualsNumeric(as.numeric(pr_marg[1,1:4]),
                     c(0.56527,0.1937941,0.2380479,0.9207503),tol=1e-4)

  #conditional occupancy
  pr_cond <- predict(fit,'state',species=1,cond=2)
  checkEqualsNumeric(as.numeric(pr_cond[1,1:4]),
                     c(0.61805,0.25368,0.089551,0.96615),tol=1e-4)

  #check newdata
  newdata <- data.frame(occ_cov1=rnorm(1),occ_cov2=rnorm(1),occ_cov3=rnorm(1))
  pr_new <- predict(fit,'state',newdata=newdata)
  checkEqualsNumeric(sapply(pr_new,function(x) x[1,1]),
                     c(0.307815,0.221529,0.0264599,0.7406483),
                     tol=1e-4)
}


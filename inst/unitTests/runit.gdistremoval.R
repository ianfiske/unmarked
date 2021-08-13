# Simulation function

simData <- function(lambda=1, sigma=40, scale=NULL, remP=0.4, remJ=4,
                    M=100, J=4,type="point", keyfun="halfnorm", T=1, phi=1)
{
    y <- array(NA, c(M, T*J))
    yRem <- matrix(NA, M, T*remJ)
    N <- rep(NA, M)

piRem <- rep(NA,remJ+1)
piRem[1] <- remP
for (j in 2:remJ){
  piRem[j] <- piRem[j-1] / remP * (1-remP) * remP
}
piRem[remJ+1] <- 1-sum(piRem, na.rm=T)

    db <- c(0, 25, 50, 75, 100)
    if(length(db)-1 != J)
      stop("hey, what")
if(keyfun=="halfnorm"){
    if(type=="point")
       g <- function(x, sig) exp(-x^2/(2*sig^2))*x
    if(type=="line")
       g <- function(x, sig) exp(-x^2/(2*sig^2))
}
if(keyfun=="exp"){
    if(type=="point")
       g <- function(x, sig) exp(-x/sig)*x
    if(type=="line")
       g <- function(x, sig) exp(-x/sig)
}
if(keyfun=="hazard"){
    shape<- sigma
    if(type=="point")
       g <- function(x, shape, scale) (1 - exp(-(r/shape)^-scale)) * x
    if(type=="line")
       g <- function(x, sig) (1 - exp(-(x/shape)^-scale) )
}
if(keyfun=="uniform"){
     if(type=="point")
       g <-  function(x,sig) x
    if(type=="line")
       g <- function(x,sig) 1
}

    cp <- u <- a <- numeric(J)
if(keyfun!="uniform"){
if(type=="point"){
   a[1] <- pi*db[2]^2
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value * 2 * pi
    for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value * 2*pi
    }
    }
if(type=="line"){
   L <-  1
   a[1] <- L*db[2]
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value
    for(j in 2:J) {
      a[j] <-  db[j+1]  - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value
    }
    }
    u <- a / sum(a)
    cp <- cp / a * u
    cp[j+1] <- 1-sum(cp)

}
if(keyfun=="uniform"){
if(type=="point"){
   a[1] <- pi*db[2]^2
   for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])

    }
    }
if(type=="line"){
   L <-  1
   a[1] <- L*db[2]

    for(j in 2:J) {
      a[j] <-  db[j+1]  - sum(a[1:j])

    }
    }
    u <- a / sum(a)
    cp <- a
    cp[j+1] <- 0

}
    for(i in 1:M) {

      N[i] <- rpois(1, lambda)

      yRem_sub <- y_sub <- c()

      for (t in 1:T){

        yRemT <- rep(0, remJ)
        yT <- rep(0, J)

        dist_class <- rem_class <- rep(NA, N[i])

        prob_in <- sum(piRem[1:remJ])*sum(cp[1:J])*phi

        if(N[i] > 0){
        for (n in 1:N[i]){

          keep <- rbinom(1,1,prob_in)
          if(keep==1){
            rem_class[n] <- sample(1:remJ, 1, prob=piRem[1:remJ]/sum(piRem[1:remJ]))
            dist_class[n] <- sample(1:J, 1, prob=cp[1:J]/sum(cp[1:J]))

            yT[dist_class[n]] <- yT[dist_class[n]] + 1
            yRemT[rem_class[n]] <- yRemT[rem_class[n]] + 1
          }
        }
        }
        yRem_sub <- c(yRem_sub, yRemT)
        y_sub <- c(y_sub, yT)
      }
      yRem[i,] <- yRem_sub
      y[i,] <- y_sub
    }

    return(list(y=y,N=N,yRem=yRem))
}

simDataRand <- function(lambda=1, lam_sd=1, groups=10,
                        sigma=40, remP=0.4, remJ=4,
                        M=100, J=4, T=1, phi=1)
{
    y <- array(NA, c(M, T*J))
    yRem <- matrix(NA, M, T*remJ)
    N <- rep(NA, M)

group_eff <- rnorm(groups, 0, lam_sd)
names(group_eff) <- letters[1:groups]
group <- sample(letters[1:groups], M, replace=T)
rand_ef <- group_eff[group]
lambda <- exp(log(lambda) + rand_ef)

piRem <- rep(NA,remJ+1)
piRem[1] <- remP
for (j in 2:remJ){
  piRem[j] <- piRem[j-1] / remP * (1-remP) * remP
}
piRem[remJ+1] <- 1-sum(piRem, na.rm=T)

    db <- c(0, 25, 50, 75, 100)
    if(length(db)-1 != J)
      stop("hey, what")
g <- function(x, sig) exp(-x^2/(2*sig^2))*x

cp <- u <- a <- numeric(J)

  a[1] <- pi*db[2]^2
  cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value * 2 * pi
    for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value * 2*pi
    }
    u <- a / sum(a)
    cp <- cp / a * u
    cp[j+1] <- 1-sum(cp)

    for(i in 1:M) {

      N[i] <- rpois(1, lambda[i])

      yRem_sub <- y_sub <- c()

      for (t in 1:T){

        yRemT <- rep(0, remJ)
        yT <- rep(0, J)

        dist_class <- rem_class <- rep(NA, N[i])

        prob_in <- sum(piRem[1:remJ])*sum(cp[1:J])*phi

        if(N[i] > 0){
        for (n in 1:N[i]){

          keep <- rbinom(1,1,prob_in)
          if(keep==1){
            rem_class[n] <- sample(1:remJ, 1, prob=piRem[1:remJ]/sum(piRem[1:remJ]))
            dist_class[n] <- sample(1:J, 1, prob=cp[1:J]/sum(cp[1:J]))

            yT[dist_class[n]] <- yT[dist_class[n]] + 1
            yRemT[rem_class[n]] <- yRemT[rem_class[n]] + 1
          }
        }
        }
        yRem_sub <- c(yRem_sub, yRemT)
        y_sub <- c(y_sub, yT)
      }
      yRem[i,] <- yRem_sub
      y[i,] <- y_sub
    }

    return(list(y=y,N=N,yRem=yRem,group=group,group_eff=group_eff))
}

test.gdistremoval.frame <- function(){
  set.seed(123)

  # Single primary period
  simdat <- simData(remJ=5)
  sc <- data.frame(sc1=rnorm(nrow(simdat$y)))
  oc <- data.frame(oc1=rnorm(length(simdat$yRem)))

  umf <- unmarkedFrameGDR(simdat$y, simdat$yRem, siteCovs=sc,
                          obsCovs=oc, dist.breaks=c(0,10,20,30,40),
                          unitsIn='m')
  checkTrue(inherits(umf, "unmarkedFrameGDR"))

  # Check subsetting
  umf_sub <- umf[1:3,]
  checkTrue(inherits(umf_sub, "unmarkedFrameGDR"))
  checkEqualsNumeric(numSites(umf_sub), 3)
  checkException(umf[,1:2])

  # Input mistake handling

  # Wrong number of dist.breaks
  checkException(unmarkedFrameGDR(simdat$y, simdat$yRem, siteCovs=sc,
                          obsCovs=oc, dist.breaks=c(0,10,20,30),
                          unitsIn='m'))

  # Wrong number of period.lengths
  checkException(unmarkedFrameGDR(simdat$y, simdat$yRem, siteCovs=sc,
                          obsCovs=oc, dist.breaks=c(0,10,20,30,40),
                          unitsIn='m', period.lengths=c(1,1,1)))

  # row sums of yDistance and yRemoval don't match
  yRem_bad <- simdat$yRem
  yRem_bad[which(yRem_bad>0)[1]] <- 0
  checkException(unmarkedFrameGDR(simdat$y, yRem_bad, siteCovs=sc,
                          obsCovs=oc, dist.breaks=c(0,10,20,30,40),
                          unitsIn='m'))

  # Multi primary period
  simdat <- simData(T=3, remJ=5)
  sc <- data.frame(sc1=rnorm(nrow(simdat$y)))
  oc <- data.frame(oc1=rnorm(length(simdat$yRem)))
  ysc <- data.frame(ysc1=rnorm(100*3))

  umf2 <- unmarkedFrameGDR(simdat$y, simdat$yRem, siteCovs=sc,
                           obsCovs=oc, dist.breaks=c(0,10,20,30,40),
                           yearlySiteCovs=ysc, unitsIn='m', numPrimary=3)

  # Check subsetting
  umf_sub <- umf2[1:3,]
  checkTrue(inherits(umf_sub, "unmarkedFrameGDR"))
  checkEqualsNumeric(numSites(umf_sub), 3)
  checkEqualsNumeric(nrow(siteCovs(umf_sub)), 3)
  checkEqualsNumeric(nrow(obsCovs(umf_sub)), 3*15)
  checkEqualsNumeric(nrow(yearlySiteCovs(umf_sub)), 3*3)
  umf_sub <- umf2[,1:2]
  checkEqualsNumeric(ncol(umf_sub@yRemoval), 10)
  checkEqualsNumeric(ncol(umf_sub@yDistance), 8)
  checkEqualsNumeric(nrow(umf_sub@obsCovs), 10*100)
  checkEqualsNumeric(nrow(umf_sub@yearlySiteCovs), 2*100)
  checkEqualsNumeric(umf_sub@numPrimary, 2)

}

test.gdistremoval <- function(){
  set.seed(123)

  sc <- data.frame(sc1=rnorm(300))
  oc <- data.frame(oc1=rnorm(5*300))

  # Half-normal
  dat <- simData(lambda=5, sigma=50, M=300, J=4, remP=0.2, remJ=5)
  umf <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')

  fit <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf)
  checkTrue(inherits(fit, "unmarkedFitGDR"))
  checkEqualsNumeric(coef(fit), c(1.7759, 0.0491, 3.8030, -1.4560, 0.0731), tol=1e-3)

  # With unequal period lengths
  umf2 <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m',
                         period.lengths=c(1,5,1,1,1))

  fit2 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf2)
  checkEqualsNumeric(coef(fit2), c(3.24070,0.04638,3.8045,-4.00957,0.09266), tol=1e-3)

  # With negative binomial
  fit3 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf,
                       mixture="NB")
  checkEqualsNumeric(coef(fit3), c(1.7771,0.0491,6.061,3.8032,-1.4585,0.0730), tol=1e-3)

  # With exponential
  set.seed(123)
  dat <- simData(lambda=5, sigma=50, M=300, J=4, remP=0.2, remJ=5, keyfun="exp")
  umf4 <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')
  fit4 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf4,
                       keyfun="exp")
  checkEqualsNumeric(coef(fit4), c(1.5876,-0.0194,3.9263,-1.3335,-0.03879), tol=1e-3)

  # With hazard
  set.seed(123)
  dat <- simData(lambda=5, sigma=50, M=300, J=4, remP=0.2, remJ=5)
  umf5 <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')
  fit5 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf5,
                       keyfun="hazard")
  checkEqualsNumeric(coef(fit5), c(1.3710,0.04585,4.0143,0.9900,-1.2050,-0.06878), tol=1e-3)

  # With uniform
  set.seed(123)
  dat <- simData(lambda=5, sigma=50, M=300, J=4, remP=0.2, remJ=5)
  umf6 <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')
  fit6 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf6,
                       keyfun="uniform")
  checkEqualsNumeric(coef(fit6), c(0.6838,0.0459,-1.207,-0.0687), tol=1e-3)

  # Methods
  gp <- getP(fit)
  checkEqualsNumeric(dim(gp$dist), c(300,4,1))
  checkEqualsNumeric(dim(gp$rem), c(300,5,1))

  s <- simulate(fit, 2)
  checkEqualsNumeric(length(s), 2)
  checkEqualsNumeric(dim(s[[1]]$yDistance), dim(fit@data@yDistance))
  checkEqualsNumeric(dim(s[[1]]$yRemoval), dim(fit@data@yRemoval))

  r <- ranef(fit)
  checkEqualsNumeric(length(bup(r)), 300)

  pb <- parboot(fit, nsim=2)
  checkTrue(inherits(pb, "parboot"))

}

test.gdistremoval.predict <- function(){
 set.seed(123)

  sc <- data.frame(sc1=rnorm(300), sc2=sample(letters[1:5],300,replace=T))
  oc <- data.frame(oc1=rnorm(5*300))

  # Half-normal
  dat <- simData(lambda=5, sigma=50, M=300, J=4, remP=0.2, remJ=5)
  umf <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')

  fit <- gdistremoval(~sc1+sc2,removalformula=~oc1+sc1,distanceformula=~1, data=umf)

  pr <- predict(fit, 'lambda')
  checkEqualsNumeric(dim(pr), c(300,4))
  nd <- data.frame(sc1=c(0,1), sc2='a')
  pr <- predict(fit, 'lambda', newdata=nd)
  checkEqualsNumeric(dim(pr), c(2,4))

  nd <- data.frame(sc1=c(0,1), sc2=letters[6])
  checkException(predict(fit, 'lambda', newdata=nd))

  pr <- predict(fit, 'rem')
  checkEqualsNumeric(dim(pr), c(5*300,4))
  nd <- data.frame(oc1=c(0,1))
  checkException(predict(fit, 'lambda', newdata=nd))
  nd <- data.frame(oc1=c(0,1), sc1=c(0,1))
  pr <- predict(fit, 'rem', newdata=nd)
  checkEqualsNumeric(dim(pr), c(2,4))

  pr <- predict(fit, 'dist')
  checkEqualsNumeric(dim(pr), c(300,4))

  checkException(predict(fit, 'test'))
}

test.gdistremoval.na <- function(){
  set.seed(123)

  sc <- data.frame(sc1=rnorm(300), sc2=sample(letters[1:5],300,replace=T))
  oc <- data.frame(oc1=rnorm(5*300))

  # Half-normal
  dat <- simData(lambda=5, sigma=50, M=300, J=4, remP=0.2, remJ=5)

  yDist <- dat$y
  yDist[1,1] <- NA
  yDist[2,] <- NA

  yRem <- dat$yRem
  yRem[3,1] <- NA
  yRem[2,] <- NA

  umf <- unmarkedFrameGDR(yDist, yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')

  fit <- gdistremoval(~1,removalformula=~1,distanceformula=~1, data=umf)
  checkEqualsNumeric(coef(fit), c(1.7084,3.81973,-1.4819), tol=1e-3)

  umf2 <- umf
  umf2@siteCovs$sc1[1] <- NA
  checkException(gdistremoval(~sc1,removalformula=~1,distanceformula=~1, data=umf2))
}

test.gdistremoval.multiperiod <- function(){
  set.seed(123)

  sc <- data.frame(sc1=rnorm(300))
  oc <- data.frame(oc1=rnorm(5*300*10))
  ysc <- data.frame(ysc1=rnorm(300*10))

  dat <- simData(lambda=5, sigma=30, M=300, T=10, J=4, remP=0.2, remJ=5, phi=0.5)
  umf <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc, yearlySiteCovs=ysc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m', numPrimary=10)

  fit <- gdistremoval(~sc1,phiformula=~ysc1, removalformula=~oc1,distanceformula=~1, data=umf)

  checkEqualsNumeric(coef(fit),
                     c(1.6884,-0.0156,0.05479,0.00021,3.39992,-1.5703,0.01918),
                     tol=1e-3)

  # Predict
  pr <- predict(fit, 'phi')
  checkEqualsNumeric(dim(pr), c(300*10,4))
  checkEqualsNumeric(as.numeric(pr[1,]), c(0.5145,0.1225,0.2884,0.7348), tol=1e-3)

  # getP
  gp <- getP(fit)
  checkEqualsNumeric(lapply(gp, dim)$phi, c(300,10))

  # ranef
  r <- ranef(fit)
  checkEqualsNumeric(dim(bup(r)), c(300,10))

}

test.gdistremoval.random <- function(){

  set.seed(123)
  dat <- simDataRand(lambda=5, lam_sd=0.4, groups=10, sigma=50, M=500, J=4, remP=0.2, remJ=5) #
  sc <- data.frame(sc1=rnorm(500), group=dat$group)
  oc <- data.frame(oc1=rnorm(5*500))
  umf <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')
  fit <- gdistremoval(~sc1+(1|group),removalformula=~oc1,distanceformula=~1, data=umf)

  checkEqualsNumeric(coef(fit), c(1.5854, -0.0452, 3.8999, -1.2795, 0.08595),
                     tol=1e-4)
  checkEqualsNumeric(sigma(fit)$sigma, 0.4080, tol=1e-4)

  pr <- predict(fit, "lambda")
  checkTrue(inherits(pr, "data.frame"))
  pr <- predict(fit, "lambda", newdata=umf@siteCovs[1:2,])
  checkTrue(inherits(pr, "data.frame"))

  s <- simulate(fit, 2)
  checkTrue(inherits(s, "list"))

  pb <- parboot(fit, nsim=2)
  checkTrue(inherits(pb, "parboot"))

}

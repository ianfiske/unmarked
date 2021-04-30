simData <- function(lambda=1, gamma=0.5, omega=0.8, sigma=40, scale=NULL,
                    M=100, T=5, J=4,type="line", keyfun="halfnorm")
{
    y <- array(NA, c(M, J, T))
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
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
    N[i,1] <- rpois(1, lambda)

    y[i,1:J,1] <- rmultinom(1, N[i,1], cp)[1:J]

    for(t in 1:(T-1)) {
        S[i,t] <- rbinom(1, N[i,t], omega)
        G[i,t] <- rpois(1, gamma)
        N[i,t+1] <- S[i,t] + G[i,t]
        y[i,1:J,t+1] <- rmultinom(1, N[i,t+1], cp)[1:J]
        }
    }
    cp <- array(cp, c(J, M, T))
    cp <- matrix(aperm(cp, c(2,1,3)), M)
    return(list(y=matrix(y, M),N=N))
}

test.unmarkedFrameDSO <- function(){
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=15,type="line",
            keyfun="halfnorm")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=15,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 100))
  checkException(unmarkedFrameDSO(y=y,numPrimary=15,
                  dist.breaks=c(0,25,50,75,100), survey="line", tlength=c(1,1)))

  checkException(unmarkedFrameDSO(y=y,numPrimary=15,
                  dist.breaks=c(0,25,50,75,100), survey="point",
                  tlength=rep(1,100)))
  checkException(unmarkedFrameDSO(y=y, numPrimary=15, dist.breaks=c(25,50,75,100),
                                  survey='line', unitsIn='m', tlength=rep(1,100)))
}

test.distsampOpen.halfnormal <- function()
{
  set.seed(456)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=15,type="line",
            keyfun="halfnorm")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=15,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 100))

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=30,keyfun="halfnorm")

  checkEqualsNumeric(coef(fm), c(1.38017,0.69961,0.053022,3.17838,0.0043299),
                     tol=1e-5)

  pr <- predict(fm, type='lambda')
  checkEqualsNumeric(as.numeric(pr[1,]),
                     c(3.9756,0.3474,3.3497,4.7183), tol=1e-4)

  pval <- getP(fm)
  checkEqualsNumeric(dim(pval), dim(y))
  checkEqualsNumeric(pval[1,1:4], c(0.211395,0.078579,0.0107615,0.0005353),
                     tol=1e-5)

  r <- residuals(fm)
  checkEqualsNumeric(dim(r), dim(y))
  checkEqualsNumeric(r[1,1:4], c(-0.84042,-0.31240,-0.042783,-0.0021283),tol=1e-4)

  ran <- ranef(fm)
  checkEqualsNumeric(bup(ran)[1,1], 2.777855, tol=1e-5)

  set.seed(123)
  sim <- simulate(fm, nsim=2)
  checkEqualsNumeric(length(sim), 2)
  checkEqualsNumeric(sim[[1]][1,1:3], c(1,0,0))

  fm2 <- update(fm, pformula=~1)
  checkEqualsNumeric(length(coef(fm2)), 4)

  #Point
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=15,type="point",
            keyfun="halfnorm")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=15,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="point",
                          unitsIn="m")

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=30,keyfun="halfnorm")
  checkEqualsNumeric(coef(fm), c(1.43259,0.82993,-0.295220,3.205348,-0.000132),
                                 tol=1e-4)

  #Check error with random effects in formula
  checkException(distsampOpen(~(1|dummy), ~1, ~1, ~1, data=umf, K=30))
}

test.distsampOpen.NA <- function(){
  set.seed(456)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=30, T=5,type="line",
            keyfun="halfnorm")$y

  y[5,1:4] <- NA
  y[2,1] <- NA
  sc = data.frame(x1=rnorm(30))
  sc$x1[3] <- NA
  ysc = data.frame(x2=rnorm(30*5))
  ysc$x2[46] <- NA

  umf <- unmarkedFrameDSO(y = y, numPrimary=5,
            siteCovs=sc,
            yearlySiteCovs=ysc,
            dist.breaks = c(0, 25, 50, 75, 100), survey="line",
            unitsIn="m",tlength=rep(1, 30))

  fm <- distsampOpen(~x1, ~x2, ~1, ~1, data=umf, K=25, keyfun="halfnorm")
  checkEqualsNumeric(coef(fm), c(1.497405,-0.0826876,-0.662144,
                                 0.651976,2.054032,3.1728838), tol=1e-4)

}

test.distsampOpen.exp <- function(){
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=10,type="line",
            keyfun="exp")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=10,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 100))

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=25,keyfun="exp")
  checkEqualsNumeric(coef(fm), c(1.34009,0.69997,-0.11887,3.17950,0.029042),
                     tol=1e-4)

  #Point
  set.seed(456)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=10,type="point",
            keyfun="exp")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=10,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="point",
                          unitsIn="m")

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=25,keyfun="exp")
  checkEqualsNumeric(coef(fm), c(1.39598,0.64463,0.053240,3.23198,0.012271),
                     tol=1e-4)

}

test.distsampOpen.unif <- function(){
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=10,type="line",
            keyfun="uniform")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=10,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 100))

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25,keyfun="unif")
  checkEqualsNumeric(coef(fm), c(1.47853,0.7475,-0.115096),
                     tol=1e-4)

  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=10,type="point",
            keyfun="uniform")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=10,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="point",
                          unitsIn="m")

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25,keyfun="unif")
  checkEqualsNumeric(coef(fm), c(1.36098,0.69191,-0.03537),
                     tol=1e-4)
}

test.distsampOpen.hazard <- function(){
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, scale=1,
               M=100, T=10,type="line", keyfun="hazard")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=10,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 100))

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=25,keyfun="hazard")
  checkEqualsNumeric(coef(fm), c(1.40843,0.64105,-0.010841,3.297099,-0.02168,0.07719),
                     tol=1e-4)
}

test.distsampOpen.NB <- function(){

  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=10,type="line",
            keyfun="exp")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=10,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 100))

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=25,keyfun="exp",
                 mixture="NB")
  checkEqualsNumeric(coef(fm), c(1.34009,0.699979,-0.118878,3.179589,
                                 0.0290938,9.238497), tol=1e-4)

}

test.distsampOpen.dynamics <- function(){
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=10,type="line",
            keyfun="uniform")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=10,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 100))

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25,keyfun="unif",
                    dynamics="notrend")
  checkEqualsNumeric(coef(fm), c(1.4080889, -0.1006024), tol=1e-5)

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25, keyfun="unif",
                     dynamics="trend")
  checkEqualsNumeric(coef(fm), c(1.518695, -0.0143889), tol=1e-5)

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25, keyfun="unif",
                     dynamics="autoreg")
  checkEqualsNumeric(coef(fm), c(1.518686, -0.018026, -5.628779), tol=1e-5)

  #Sketchy estimates
  #Maybe just because data were simulated using a different process?
  #Leaving these in for now just to make sure they run without errors
  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25, keyfun="unif",
                     dynamics="gompertz")

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25, keyfun="unif",
                     dynamics="ricker")

}

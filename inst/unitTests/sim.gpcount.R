library(unmarked)


sim1 <- function(R=50, J=3, K=3, lambda=5, phi=0.6, p=0.4) {
    M <- rpois(R, lambda) # super-population size
    N <- matrix(NA, R, J) # Population available
    y <- array(NA, c(R, K, J)) # Detected
    for(i in 1:R) {
        for(j in 1:J) {
            N[i,j] <- rbinom(1, M[i], phi)
            y[i,,j] <- rbinom(K, N[i,j], p)
        }
    }
    y <- matrix(y, R)
    return(list(y=y, N=N))
}

set.seed(348)
y1 <- sim1()$y

y1[1,] <- NA
y1[2, 1:3] <- NA
y1[3, 4:6] <- NA
umf <- unmarkedFrameGPC(y=y1, numPrimary=3)

fm1.1 <- gpcount(~1, ~1, ~1, umf, K=40, control=list(trace=TRUE, REPORT=1))
fm1.2 <- gpcount(~1, ~1, ~1, umf, K=40, mixture="NB",
                 control=list(trace=TRUE, REPORT=1))
fm1.3 <- gpcount(~1, ~1, ~1, umf, K=40, mixture="ZIP",
                 control=list(trace=TRUE, REPORT=1))


nsim1 <- 50
simout1 <- matrix(NA, nsim1, 3)
lam1 <- 5
phi1 <- 0.5
p1 <- 0.4
nPrimary1 <- 3
for(i in 1:nsim1) {
    if(i %% 10 == 0) cat("doing", i, "\n")
    sim1.i <- sim1(lambda=lam1, phi=phi1, p=p1, J=nPrimary1)
    umf1.i <- unmarkedFrameGPC(y=sim1.i, numPrimary=nPrimary1)
    fm1.i <- gpcount(~1, ~1, ~1, umf1.i, K=50, engine="C")
    mle1.i <- coef(fm1.i)
    simout1[i,] <- c(exp(mle1.i[1]), plogis(mle1.i[2:3]))
}

hist(simout1[,1]); abline(v=lam1, lwd=2, col=4)
hist(simout1[,2]); abline(v=phi1, lwd=2, col=4)
hist(simout1[,3]); abline(v=p1, lwd=2, col=4)














# Covariates

set.seed(568)
R <- 50
J <- 4
K <- 3
x1 <- rnorm(R)
x2 <- matrix(rnorm(R*J), R, J)
x3 <- matrix(rnorm(R*K*J), R, K*J)

sim2 <- function(x1, x2, x3,
                 lam0=0, lam1=1, phi0=1, phi1=1, p0=0, p1=1) {
    R <- length(x1)
    J <- ncol(x2)
    K <- ncol(x3)/J
    lambda <- exp(lam0 + lam1*x1)
    phi <- plogis(phi0 + phi1*x2)
    p <- plogis(p0 + p1*x3)
    p <- array(p, c(R, K, J))
    M <- rpois(R, lambda) # super-population size
    N <- matrix(NA, R, J) # Population available
    y <- array(NA, c(R, K, J)) # Detected
    for(i in 1:R) {
        for(j in 1:J) {
            N[i,j] <- rbinom(1, M[i], phi[i,j])
            y[i,,j] <- rbinom(K, N[i,j], p[i,,j])
        }
    }
    y <- matrix(y, R)
    return(list(y=y, N=N))
}

y2 <- sim2(x1, x2, x3)$y
umf2 <- unmarkedFrameGPC(y=y2,
                         siteCovs=data.frame(x1),
                         yearlySiteCovs=list(x2=x2),
                         obsCovs = list(x3=x3), numPrimary=J)
summary(umf2)

fm2.1 <- gpcount(~x1, ~x2, ~x3, umf2, K=40,
                 control=list(trace=TRUE, REPORT=1))

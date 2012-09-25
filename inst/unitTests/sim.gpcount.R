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
for(i in 1:nsim1) {
    if(i %% 10 == 0) cat("doing", i, "\n")
    lam1 <- 5
    phi1 <- 0.5
    p1 <- 0.4
    nPrimary1 <- 3
    sim1.i <- sim1(lambda=lam1, phi=phi1, p=p1, J=nPrimary1)
    umf1.i <- unmarkedFrameGPC(y=y1, numPrimary=nPrimary1)
    fm1.i <- gpcount(~1, ~1, ~1, umf1.i, K=50, engine="C")
    mle1.i <- coef(fm1.i)
    simout1[i,] <- c(exp(mle1.i[1]), plogis(mle1.i[2:3]))
}


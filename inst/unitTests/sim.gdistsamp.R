

library(unmarked)

sim <- function(lambda=5, phi=0.5, shape=20, scale=10, R=100, T=3,
    breaks=seq(0, 50, by=10), survey="pt", detfun="hn")
{
    J <- length(breaks)-1
    maxDist <- max(breaks)
    tlength <- 1000
    switch(survey,
           pt = A <- (2*maxDist)^2 / 10000,   # Area (ha) of sqauare
           line = A <- maxDist*2*100 / 10000 # Area (ha) 100m transect
           )
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda * A) # Individuals within the rectangle
        N <- rbinom(T, M, phi)    # Individuals available at time t
        for(t in 1:T) {
            switch(survey,
                pt = {
                    X <- runif(N[t], -maxDist, maxDist)
                    Y <- runif(N[t], -maxDist, maxDist)
                    d <- sqrt(X^2+Y^2)
                    },
                line = {
                    d <- runif(N[t], 0, maxDist)
                    })

            # Detection process
            if(length(d) > 0) {
                switch(detfun,
                    hn   = p <- exp(-d^2 / (2 * shape^2)),
                    exp  = p <- exp(-d/shape),
                    haz  = p <- 1-exp(-(d/shape)^-scale),
                    unif = p <- 1
                    )
                d1 <- d[rbinom(length(d), 1, p) == 1]
                y[i,,t] <- table(cut(d1, breaks, include.lowest=TRUE))
                }
            }
        }
    y <- matrix(y, nrow=R)
    return(y)
}








library(unmarked)

breaks <- seq(0, 50, by=10)
T <- 3
set.seed(3)
umf <- unmarkedFrameGDS(y = sim(T=T, breaks=breaks), survey="point",
    unitsIn="m", dist.breaks=breaks, numPrimary=T)
summary(umf)

system.time(m <- gdistsamp(~1, ~1, ~1, umf)) # 28s

backTransform(m, type="lambda")
backTransform(m, type="phi")
backTransform(m, type="det")



# Point-transect, half-normal
nsim <- 5
simout <- matrix(NA, nsim, 3)
colnames(simout) <- c('lambda', 'phi', 'sigma')
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y <- sim(lambda=20, phi=0.7, R=100, T=T, breaks=breaks)
    umf <- unmarkedFrameGDS(y = y, survey="point",
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf, output="density")
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3]))
    }

par(mfrow=c(3, 1))
hist(simout[,1], xlab=expression(lambda), main="")
abline(v=20, col=4)
hist(simout[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout[,3], xlab=expression(sigma), main=""); abline(v=20, col=4)



# Point-transect, neg exp
nsim <- 5
simout <- matrix(NA, nsim, 3)
colnames(simout) <- c('lambda', 'phi', 'rate')
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y <- sim(lambda=20, phi=0.6, R=100, T=T, breaks=breaks, detfun="exp")
    umf <- unmarkedFrameGDS(y = y, survey="point",
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf, keyfun="exp", output="density", K=100)
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3]))
    }

par(mfrow=c(3, 1))
hist(simout[,1], xlab=expression(lambda), main="")
abline(v=20, col=4, lwd=2)
hist(simout[,2], xlab=expression(phi), main="")
abline(v=0.7, col=4, lwd=2)
hist(simout[,3], xlab=expression(rate), main="")
abline(v=20, col=4, lwd=2)




# Point-transect, hazard
nsim <- 100
simout <- matrix(NA, nsim, 4)
colnames(simout) <- c('lambda', 'phi', 'shape', 'scale')
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y <- sim(lambda=20, phi=0.7, shape=20, scale=1, R=100, T=T,
             breaks=breaks, detfun="haz")
    umf <- unmarkedFrameGDS(y = y, survey="point",
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf, keyfun="hazard", output="density",
                   K=100, starts=c(3, 0, 2, 0))
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3:4]))
    }

par(mfrow=c(3, 1))
hist(simout[,1], xlab=expression(lambda), main="")
abline(v=20*(pi*50^2/10000), col=4)
hist(simout[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout[,3], xlab=expression(sigma), main=""); abline(v=20, col=4)




# Point-transect, uniform
nsim <- 5
simout <- matrix(NA, nsim, 2)
colnames(simout) <- c('lambda', 'phi')
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y <- sim(lambda=20, phi=0.6, R=100, T=T, breaks=breaks, detfun="unif",
             survey="pt")
    umf <- unmarkedFrameGDS(y = y, survey="point",
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf, keyfun="uniform", output="density",
                   K=100, unitsOut="ha")
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2]))
    }

par(mfrow=c(3, 1))
hist(simout[,1], xlab=expression(lambda), main="")
abline(v=20, col=4)
hist(simout[,2], xlab=expression(phi), main=""); abline(v=0.5, col=4)










# Line-transect, half-normal
nsim <- 100
simout <- matrix(NA, nsim, 3)
colnames(simout) <- c('lambda', 'phi', 'sigma')
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y <- sim(lambda=20, phi=0.7, R=100, T=T, breaks=breaks, survey="line")
    umf <- unmarkedFrameGDS(y = y, survey="line", tlength=rep(100,100),
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf,
        starts=c(1.5, 0.5, 3), output="density", K=100)
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3]))
    }

par(mfrow=c(3, 1))
hist(simout[,1], xlab=expression(lambda), main="")
abline(v=20, col=4)
hist(simout[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout[,3], xlab=expression(sigma), main=""); abline(v=20, col=4)



# Line-transect, neg exp
nsim <- 5
simout <- matrix(NA, nsim, 3)
colnames(simout) <- c('lambda', 'phi', 'rate')
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y <- sim(lambda=20, phi=0.7, R=100, T=T, breaks=breaks, detfun="exp",
             survey="line")
    umf <- unmarkedFrameGDS(y = y, survey="line", tlength=rep(1000,100),
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf, rel.tol=0.01,
        starts=c(1.5, 0.5, 3), keyfun="exp", output="density")
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3]))
    }

par(mfrow=c(3, 1))
hist(simout[,1], xlab=expression(lambda), main="")
abline(v=20*(pi*50^2/10000), col=4)
hist(simout[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout[,3], xlab=expression(sigma), main=""); abline(v=20, col=4)




# Line-transect, hazard
nsim <- 5
simout <- matrix(NA, nsim, 4)
colnames(simout) <- c('lambda', 'phi', 'shape', 'scale')
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y <- sim(lambda=20, phi=0.7, shape=20, scale=1, R=100, T=T,
             breaks=breaks, detfun="haz", survey="line")
    umf <- unmarkedFrameGDS(y = y, survey="line", tlength=rep(100,100),
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf, rel.tol=0.01, keyfun="hazard",
        starts=c(1.5, 0.5, 3, 0), output="density")
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2]), exp(e[3:4]))
    }

par(mfrow=c(3, 1))
hist(simout[,1], xlab=expression(lambda), main="")
abline(v=20*(pi*50^2/10000), col=4)
hist(simout[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)
hist(simout[,3], xlab=expression(sigma), main=""); abline(v=20, col=4)




# Line-transect, uniform
nsim <- 10
simout <- matrix(NA, nsim, 2)
colnames(simout) <- c('lambda', 'phi')        )
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    breaks <- seq(0, 50, by=10)
    T <- 5
    y <- sim(lambda=20, phi=0.7, R=100, T=T, breaks=breaks, detfun="unif",
             survey="line")
    umf <- unmarkedFrameGDS(y = y, survey="line", tlength=rep(100,100),
        unitsIn="m", dist.breaks=breaks, numPrimary=T)
    m <- gdistsamp(~1, ~1, ~1, umf, keyfun="uniform", output="density",
                   K=100, starts=c(3, 0))
    e <- coef(m)
    simout[i,] <- c(exp(e[1]), plogis(e[2]))
    }

par(mfrow=c(3, 1))
hist(simout[,1], xlab=expression(lambda), main="")
abline(v=20, col=4)
hist(simout[,2], xlab=expression(phi), main=""); abline(v=0.7, col=4)








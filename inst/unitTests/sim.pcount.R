


sim1 <- function(R=100, J=5, lambda=2, p=0.3, mix="P", disp=1) {
    y <- matrix(NA, R, J)
    switch(mix, 
        P = N <- rpois(R, lambda),
        NB = N <- rnbinom(R, mu=lambda, size=disp)
        )
    for(i in 1:R) {
        y[i,] <- rbinom(J, N[i], p)
        }
    return(y)
    }

set.seed(7)

umf1 <- unmarkedFramePCount(y=sim1(mix="NB", disp=0.7))
m1 <- pcount(~1 ~1, umf1, K=10)
coef(m1)

m2 <- pcount(~1 ~1, umf1, K=100)
coef(m2)


m3 <- pcount(~1 ~1, umf1, K=10, mixture="NB")
coef(m3)

m4 <- pcount(~1 ~1, umf1, K=100, mixture="NB")
coef(m4)

m5 <- pcount(~1 ~1, umf1, K=500, mixture="NB")
coef(m5)




set.seed(11)
nsims <- 50
simout1 <- matrix(NA, nsims, 2)
lam <- 20
sig <- 30
for(i in 1:nsims) {
    cat("sim", i, "\n"); flush.console()
    umf <- unmarkedFrameDS(y = simPt(lambda=lam, sigma=sig), survey="point", 
        dist.breaks=seq(0, 50, by=10), unitsIn="m")
    m <- distsamp(~1 ~1, umf, starts=c(log(lam), log(sig)), rel.tol=1e-3, 
        output="abund")
    simout1[i,] <- exp(coef(m))
    }
hist(simout1[,1]); abline(v=lam*pi*50^2/10000, lwd=2, col=3)
hist(simout1[,2]); abline(v=sig, lwd=2, col=3)





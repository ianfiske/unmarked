
sim7 <- function(lambda=1, gamma=0.5, omega=0.8, sigma=40, scale=NULL, M=100, T=5,
                 J=4, type="line", keyfun="halfnorm") {
    keyfun <- "halfnorm" ## FIXME!!!
    y <- array(NA, c(M, J, T))
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    db <- c(0, 25, 50, 75, 100)
    w <- diff(db) ## rbc added
    if(length(db)-1 != J)
        stop("hey, what")
    if(keyfun=="halfnorm"){
        if(type=="point")
            g <- function(x, sig) exp(-x^2/(2*sig^2))*x
        if(type=="line")
            g <- function(x, sig) exp(-x^2/(2*sig^2))
    }

    cp <- u <- a <- numeric(J)
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
        cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value #/ w[1]
        for(j in 2:J) {
            a[j] <-  db[j+1]  - sum(a[1:j])
            cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value #/ w[j]
        }
    }
    u <- a / sum(a)
    cp <- cp / a * u
    cp[j+1] <- 1-sum(cp)

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
    cat("max(N) =", max(N), "\n")
    return(list(y=matrix(y, M),N=N, cp=cp))
}

library(unmarked)
set.seed(711)
lambda <- 4
gamma <- 2
omega <- 0.5
sigma <- 30
T <- 10


out7 <- sim7(lambda, gamma, omega, sigma=50, M=200, T=T,type="line", keyfun="halfnorm")

cbind(out7$y[,1:4], out7$N[,1])

cbind(rowSums(out7$y[,1:4]), out7$N[,1])

cbind(out7$y[,5:8], out7$N[,2])

cbind(rowSums(out7$y[,5:8]), out7$N[,2])

colSums(out7$y)

y.sim7 <- sim7(lambda, gamma, omega, sigma=sigma, M=50, T=T,type="line",
               keyfun="halfnorm")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=120,keyfun="half",
##                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 starts=c(0, 0, 0, 3), method="BFGS",
                 se=FALSE, nintervals=5, control=list(trace=TRUE, REPORT=1, maxit=100))

exp(coef(fm)[c(1,2,4)])
plogis(coef(fm)[3])



y.sim7 <- sim7(lambda, gamma, omega, sigma=sigma, M=50, T=T,type="point",
               keyfun="halfnorm")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="point", unitsIn="m")
fm <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=120,keyfun="half",
##                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 starts=c(0, 0, 0, 3),
                 se=FALSE, nintervals=5, control=list(trace=TRUE, REPORT=1))





y.sim7 <- sim7(lambda, gamma, omega, sigma=25, M=200, T=T,type="line", keyfun="exp")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm2 <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=60,keyfun="exp",
                 starts=rnorm(4,0,.1)+c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 se=FALSE, nintervals=8)


y.sim7 <- sim7(lambda, gamma, omega, sigma=NULL, M=200, T=T,type="line", keyfun="uniform", scale=1)$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm3 <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=60,keyfun="uniform",
                 starts=c(log(c(lambda, gamma)),plogis(omega), log(sigma)),
                 se=FALSE, nintervals=8)



y.sim7 <- sim7(lambda, gamma, omega, sigma=20, scale = 1, M=200, T=T,type="line", keyfun="hazard")$y
umf7b <- unmarkedFrameDSO(y = y.sim7, numPrimary=T,
    dist.breaks = c(0, 25, 50, 75, 100), survey="line", unitsIn="m",tlength=rep(1, 200))
fm4 <- distsampOpen(~1, ~1, ~1, ~1, data = umf7b, K=50,keyfun="hazard", method="Nelder-Mead",
   control=list(trace=2),
  starts=rnorm(5,
         0.8*c(log(c(lambda, gamma)),plogis(omega), log(sigma), 0),0),
                 se=FALSE, nintervals=8)

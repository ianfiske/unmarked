
library(unmarked)

## Simulate no covariates, constant sampling period intervals	

sim1 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=50, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}
                            
# Prepare data                               
set.seed(3)
umf <- unmarkedFramePCO(y = sim1())

summary(umf)

# Fit model and backtransform
system.time(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10))  # 14s on 64bit.

backTransform(m1, "lambda") # 0.71 initial abundance
backTransform(m1, "gamma")  # 0.45 recruitment rate
backTransform(m1, "omega")  # 0.75 survival rate
backTransform(m1, "det")    # 0.72 detection probability


set.seed(3223)
nsim1 <- 50
simout1 <- matrix(NA, nsim1, 4)
colnames(simout1) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim1) {
    cat("sim", i, "\n"); flush.console()
    lambda <- 1
    gamma <- 0.5
    omega <- 0.8
    p <- 0.7
    y.sim1 <- sim1(lambda, gamma, omega, p, M=100)
    umf1 <- unmarkedFramePCO(y = y.sim1)
    m1 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=15, 
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m1)
    simout1[i, 1:2] <- exp(e[1:2])
    simout1[i, 3:4] <- plogis(e[3:4])
    }

png("pcountOpenSim1.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout1[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout1[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout1[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout1[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)    
dev.off()    





## Simulate covariate model with constant intervals



sim2 <- function(lam=c(0,1), gam=c(-1,-1), om=c(2,-1), p=c(-1,1), M=50, T=5)
{
    y <- gamma <- omega <- det <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    veght <- rnorm(M)
    isolation <- matrix(rnorm(M*T), M, T)
    time <- matrix(rnorm(M*T, 1), M, T)
    lambda <- exp(lam[1] + lam[2]*veght)
    gamma[] <- exp(gam[1] + gam[2]*isolation)
    omega[] <- plogis(om[1] + om[2]*isolation)
    det[] <- plogis(p[1] + p[2]*time)

    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega[,t])
        G[,t] <- rpois(M, gamma[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, det)
    return(list(y=y, covs=data.frame(veght=veght, 
        isolation=isolation, time=time)))
}
                            
sim2()



nsim2 <- 50
simout2 <- matrix(NA, nsim2, 8)
colnames(simout2) <- c('lam0', 'lam1', 'gam0', 'gam1', 'om0', 'om1', 'p0', 'p1')
for(i in 1:nsim2) {
    cat("sim", i, "\n"); flush.console()
    lam <- c(-2, 1)
    gam <- c(-1, -1)
    om <- c(0, -1)
    p <- c(-1, 1)
    sim2out <- sim2(lam, gam, om, p)
    y.sim2 <- sim2out$y
    covs <- sim2out$covs
    cn <- colnames(covs)
    siteCovs <- covs[,grep("veght", cn), drop=FALSE]
    obsCovs <- list(time = covs[,grep("time", cn)], 
        isolation = covs[,grep("isolation", cn)])     
    umf2 <- unmarkedFramePCO(y = y.sim2, siteCovs=siteCovs, obsCovs=obsCovs)
    m2 <- pcountOpen(~veght, ~isolation, ~isolation, ~time, umf2, K=40, se=F, 
        starts=c(lam, gam, om, p))
    e <- coef(m2)
    simout2[i, ] <- e
    }

png("pcountOpenSim2.png", width=6, height=8, units="in", res=360)
par(mfrow=c(4,2))
hist(simout2[,1], xlab=expression(lambda)); abline(v=lam[1], lwd=2, col=4)
hist(simout2[,2], xlab=expression(lambda)); abline(v=lam[2], lwd=2, col=4)
hist(simout2[,3], xlab=expression(gamma)); abline(v=gam[1], lwd=2, col=4)
hist(simout2[,4], xlab=expression(gamma)); abline(v=gam[2], lwd=2, col=4)
hist(simout2[,5], xlab=expression(omega)); abline(v=om[1], lwd=2, col=4)
hist(simout2[,6], xlab=expression(omega)); abline(v=om[2], lwd=2, col=4)
hist(simout2[,7], xlab=expression(p)); abline(v=p[1], lwd=2, col=4)    
hist(simout2[,8], xlab=expression(p)); abline(v=p[2], lwd=2, col=4)     
dev.off()













## Simulate uneven sampling period intervals
set.seed(333)

sim3 <- function(lambda=4, gamma=0.1, omega=0.8, p=0.7, M=50, T=5)
{
    y <- N <- date <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    date[,1] <- 1
    
    for(i in 1:M) {
    for(t in 2:T) {
        delta <- max(rpois(1, 3), 1)
        date[i, t] <- date[i, t-1] + delta
        S[i, t-1] <- rbinom(1, N[i, t-1], omega)
        G[i, t-1] <- rpois(1, gamma)
        N[i, t] <- S[i, t-1] + G[i, t-1]
        if(delta > 1) {
            for(d in 2:delta) {
                S[i, t-1] <- rbinom(1, N[i, t], omega)
                G[i, t-1] <- rpois(1, gamma)
                N[i, t] <- S[i, t-1] + G[i, t-1]
                }
            }
        }}
    y[] <- rbinom(M*T, N, p)
    mode(date) <- "integer"
    return(list(y=y, dates=date))
}

yd <- sim3()
y <- yd$y
dates <- yd$dates

formatDelta <- unmarked:::formatDelta


# y[1, 1] <- y[2,1] <- y[3,2] <- NA

head(y)
head(dates)
head(formatDelta(dates, y))
colSums(y)



# Prepare data
umfO <- unmarkedFramePCO(y = y, dates = dates)
umfO
max(y)

# Fit model
system.time(m3 <- pcountOpen(~1, ~1, ~1, ~1, umfO, K=15, se=TRUE, 
    starts=c(1.4, -0.7, 0.6, 0.6),
    control=list(maxit=50, trace=T, REPORT=1))) # 30 min using pure R
                                                #  1 min using C++
backTransform(m3, "lambda")
backTransform(m3, "gamma")
backTransform(m3, "omega")
backTransform(m3, "det")









set.seed(3223)
nsim3 <- 50
simout3 <- matrix(NA, nsim3, 4)
colnames(simout3) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim3) {
    cat("sim3", i, "\n"); flush.console()
    lambda <- 1
    gamma <- 0.5
    omega <- 0.8
    p <- 0.7
    yd <- sim3(lambda, gamma, omega, p, M=100)
    y.sim3 <- yd$y
    dates3 <- yd$dates
    umf3 <- unmarkedFramePCO(y = y.sim3, dates=dates3)
    m3 <- pcountOpen(~1, ~1, ~1, ~1, umf3, K=15, 
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m3)
    simout3[i, 1:2] <- exp(e[1:2])
    simout3[i, 3:4] <- plogis(e[3:4])
    }



png("pcountOpenSim3.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout3[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout3[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout3[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout3[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)    
dev.off()    














# Auto-regressive model




sim4 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=50, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma*N[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}
                            
# Prepare data
set.seed(434)
umf4 <- unmarkedFramePCO(y = sim4())

summary(umf4)

# Fit model and backtransform
system.time(m4 <- pcountOpen(~1, ~1, ~1, ~1, umf4, K=20,
    dynamics="autoreg")) # 52s on 64bit.

backTransform(m4, "lambda") # 1.1 initial abundance
backTransform(m4, "gamma") # 0.47 recruitment rate
backTransform(m4, "omega") # 0.84 survival rate
backTransform(m4, "det") # 0.76 detection probability


set.seed(3223)
nsim4 <- 5
simout4 <- matrix(NA, nsim4, 4)
colnames(simout4) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim4) {
    cat("sim4", i, "\n"); flush.console()
    lambda <- 1
    gamma <- 0.5
    omega <- 0.7
    p <- 0.7
    y.sim4 <- sim4(lambda, gamma, omega, p, M=100)
    umf4 <- unmarkedFramePCO(y = y.sim4)
    m4 <- pcountOpen(~1, ~1, ~1, ~1, umf4, K=30, dynamics="autoreg",
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m4)
    simout4[i, 1:2] <- exp(e[1:2])
    simout4[i, 3:4] <- plogis(e[3:4])
    }

png("pcountOpenSim4.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout4[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout4[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout4[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout4[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
dev.off()










# No trend model




sim5 <- function(lambda=1, omega=0.8, p=0.7, M=50, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    gamma <- (1-omega)*lambda
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma*N[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}
                            
# Prepare data                               
set.seed(434)
umf4 <- unmarkedFramePCO(y = sim4())

summary(umf4)

# Fit model and backtransform
system.time(m4 <- pcountOpen(~1, ~1, ~1, ~1, umf4, K=20, 
    dynamics="autoreg"))  # 52s on 64bit.

backTransform(m4, "lambda") # 1.1 initial abundance
backTransform(m4, "gamma")  # 0.47 recruitment rate
backTransform(m4, "omega")  # 0.84 survival rate
backTransform(m4, "det")    # 0.76 detection probability


set.seed(3223)
nsim4 <- 5
simout4 <- matrix(NA, nsim4, 4)
colnames(simout4) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim4) {
    cat("sim4", i, "\n"); flush.console()
    lambda <- 1
    gamma <- 0.5
    omega <- 0.7
    p <- 0.7
    y.sim4 <- sim4(lambda, gamma, omega, p, M=100)
    umf4 <- unmarkedFramePCO(y = y.sim4)
    m4 <- pcountOpen(~1, ~1, ~1, ~1, umf4, K=30, dynamics="autoreg",
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m4)
    simout4[i, 1:2] <- exp(e[1:2])
    simout4[i, 3:4] <- plogis(e[3:4])
    }

png("pcountOpenSim4.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout4[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout4[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout4[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout4[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)    
dev.off()    








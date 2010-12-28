
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
    y.sim1 <- sim1(lambda, gamma, omega, p)
    umf1 <- unmarkedFramePCO(y = y.sim1)
    m1 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=10, 
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
    om <- c(1, -1)
    p <- c(-1, 1)
    sim2out <- sim2(lam, gam, om, p)
    y.sim2 <- sim2out$y
    covs <- sim2out$covs
    cn <- colnames(covs)
    siteCovs <- covs[,grep("veght", cn), drop=FALSE]
    obsCovs <- list(time = covs[,grep("time", cn)], 
        isolation = covs[,grep("isolation", cn)])     
    umf2 <- unmarkedFramePCO(y = y.sim2, siteCovs=siteCovs, obsCovs=obsCovs)
    m2 <- pcountOpen(~veght, ~isolation, ~isolation, ~time, umf2, K=10, se=F, 
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
M <- 25
T <- 5
date <- matrix(NA, M, T)
date[,1] <- pmax(rpois(M, 2), 1)
for(t in 2:T) {
    date[,t] <- date[,t-1] + pmax(rpois(M, 10), 1)
    }
datediff <- t(apply(date, 1, diff))    
lambda <- 4
gamma <- 0.05   # Daily survival rate
omega <- 0.95   # Daily recruitment rate
p <- 0.7
y <- N <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], omega^datediff[,t])   
	G[,t] <- rpois(M, gamma*datediff[,t])           
	N[,t+1] <- S[,t] + G[,t]
	}
y[] <- rbinom(M*T, N, p)


formatDelta <- unmarked:::formatDelta


y[1, 1] <- y[2,1] <- y[3,2] <- NA

head(y)
head(date)
head(formatDelta(date, y))



# Prepare data
umfO <- unmarkedFramePCO(y = y, delta = date)
umfO

# Fit model
(m3 <- pcountOpen(~1, ~1, ~1, ~1, umfO, K=15, se=TRUE, 
    starts=c(1.3, -3, 5, 0.5),
    control=list(maxit=30, trace=T, REPORT=1)))
backTransform(m3, "lambda")
backTransform(m3, "gamma")
backTransform(m3, "omega")
backTransform(m3, "det")





# Auto-regressive model


set.seed(3)
M <- 50
T <- 5
lambda <- 1
gamma <- 0.5
omega <- 0.8
p <- 0.7
y <- N <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], omega)
	G[,t] <- rpois(M, gamma * N[,t])
	N[,t+1] <- S[,t] + G[,t]
	}
y[] <- rbinom(M*T, N, p)


# Prepare data
umf4 <- unmarkedFramePCO(y = y)
umf4

# Fit model
(m4 <- pcountOpen(~1, ~1, ~1, ~1, umf4, K=20, dynamics="autoreg", se=TRUE, 
    starts=c(0, 0.5, 0.5, 0.6), 
    control=list(maxit=50, trace=T, REPORT=1)))
backTransform(m4, "lambda") # 0.79 initial abundance
backTransform(m4, "gamma")  # 0.43 recruitment rate
backTransform(m4, "omega")  # 0.88 survival rate
backTransform(m4, "det")    # 0.60 detection probability


       
